%% dataImport %%

%% Import Data
% Uncomment the csvReadIsodat line in this section when running the script 
% for the first time in a session. It's easier and faster to then comment
% these lines back out and not clear the importedData variable as the files
% take a long time to load in.

% set(0,'defaultFigureVisible','off'); disp('Run dataImport: Turning Figures Off');
clearvars -EXCEPT importedData;
cluk; clc;

filesToImport = ["XP-2018(excelExportIntensityJDM).csv"; ...
    "XP-2017(excelExportIntensityJDM).csv"; ...
    "XP-2016(excelExportIntensityJDM).csv"];

%importedData = csvReadIsodat(filesToImport);

%% Calculate Delta Values and Compile Useful Metadata
% NOTE: The current method of interpolating onto the ~isRef indices
% excludes the CO2 checks from the blocks I consider as they are recorded
% in isRef as 'true'. This is probably the right way to do it but I have to
% add them back in below, including their metadata.

[cycle_deltas,cycle_metadata] = calcCycles(importedData,'IsRef__');

delta_cols = string(cycle_deltas.Properties.VariableNames);
delta_labels = string(cycle_deltas.Properties.VariableDescriptions);
cycle_deltas = table2array(cycle_deltas);

metadata_fields = fieldnames(cycle_metadata)'; % Transpose it to a row vector so that it works as a loop index

%% Reject Erroneous Cycles
% Get rid of some data that is of little use for evaluating the health of
% the mass spec and which causes problems with PIS etc.
%
% First, discard all data with a gas config different to 'Air+'
% This is mostly Jeff's neon experiments.
    % N.B. Note that there is a difference between Gas Config and Gas Name!
    % N.B. This does not reject the data from the CO2 check blocks as they
    % are already misisng from the cycle deltas and metadata.

iAirPlus = cycle_metadata.gasConfig == 'Air+';
cycle_deltas(~iAirPlus,:) = [];
for ii = 1:numel(metadata_fields)
    cycle_metadata.(metadata_fields{ii})(~iAirPlus,:) = [];
end

% Could do something similar here with folder names e.g. only use data from
% the "SPICE" Folders:
% importedData = importedData(contains(importedData.FileHeader_Filename,'SPICE'),:);


% Second, reject cycles where the mass 28 beam is saturated on the cup
% (>9000 mV) or where it is less than 10 mV, presumably indicating some 
% problem with the source.
    % N.B. This also happens sometimes when measuring in an unusual gas
    % configuration that doesn't focus the 28 beam into a cup but collects 
    % the signal anyway. If the 28 beam isn't included in the gas config 
    % then its value is recorded as NaN.

iBeamReject = cycle_metadata.int28SA > 9000 | cycle_metadata.int28SA < 10 | cycle_metadata.int28ST > 9000 | cycle_metadata.int28ST < 10;

cycle_deltas = cycle_deltas(~iBeamReject,:);
for ii = 1:numel(metadata_fields)
    cycle_metadata.(metadata_fields{ii}) = cycle_metadata.(metadata_fields{ii})(~iBeamReject);
end

%% Reshape into Isotope Ratio-x-Block-x-Cycle
% Identify the different blocks by the unique filenames
[~,idx_blocksAll,~] = unique(cycle_metadata.filename,'stable');
blockLengthsAll = diff(idx_blocksAll); blockLengthsAll(end+1)=length(cycle_deltas)-(idx_blocksAll(end)-1);

% TAKE ONLY THE BLOCKS WITH 16 CYCLES! - This causes me to lose 156 blocks
% (17672 -> 17516, <1%), mostly with <5 cycles in them.
idx_blocks = idx_blocksAll(blockLengthsAll==16);
blockLengths = blockLengthsAll(blockLengthsAll==16);

numberOfBlocks = length(idx_blocks);
longestBlock = max(blockLengths);

% Fill the Array of Delta Values
block_deltas = nan(numel(delta_cols),numberOfBlocks,longestBlock);

for ii = 1:length(idx_blocks)
    block_deltas(:,ii,1:blockLengths(ii)) = cycle_deltas(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:)';
end

for ii = 1:length(idx_blocks)
    for jj = 1:numel(metadata_fields)
        block_metadata.(metadata_fields{jj})(ii,1:blockLengths(ii)) = cycle_metadata.(metadata_fields{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end

%% Reshape into Cycles-x-Isotope Ratio-x-Block-x-Sample Aliquot
% Define Methods that Are Run at the Start of Each New Sample
sampleStartMethods = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];

idx_working = false(size(block_metadata.msDatetime,1),1);

for ii=1:length(sampleStartMethods)
    idx_working = idx_working | block_metadata.method(:,1)==sampleStartMethods(ii);
end

% Find new aliquots by finding above methods or by finding the start of a
% new sequence (sequence row decreases from N to 1).
idx_SampleAliquotsAll = find(idx_working | ([0; diff(block_metadata.sequenceRow(:,1))]<0));
aliquotLengthsAll = diff(idx_SampleAliquotsAll); aliquotLengthsAll(end+1)=length(block_deltas)-(idx_SampleAliquotsAll(end)-1);

% TAKE ONLY THE ALIQUOTS WITH 4 OR 5 BLOCKS! - This causes me to lose 56
% aliquots (4360 -> 4304, ~1%), mostly with only 1 or 2 blocks in them.
idx_SampleAliquots = idx_SampleAliquotsAll(aliquotLengthsAll>3 & aliquotLengthsAll<6);
aliquotLengths = aliquotLengthsAll(aliquotLengthsAll>3 & aliquotLengthsAll<6);

numberOfAliquots = length(idx_SampleAliquots);
longestAliquot = max(aliquotLengths);

% Reshape
aliquot_deltas = nan(numberOfAliquots,numel(delta_cols),longestAliquot,longestBlock);

for ii = 1:length(idx_SampleAliquots)
    aliquot_deltas(ii,:,1:aliquotLengths(ii),:) = block_deltas(:,idx_SampleAliquots(ii):idx_SampleAliquots(ii)+aliquotLengths(ii)-1,:);
    for jj=1:numel(metadata_fields)
        aliquot_metadata.(metadata_fields{jj})(ii,1:aliquotLengths(ii),:) = block_metadata.(metadata_fields{jj})(idx_SampleAliquots(ii):idx_SampleAliquots(ii)+aliquotLengths(ii)-1,:);
    end
end

%% Add the CO2 Check Data Back In
% Find the Timestamp for the CO2 Check Blocks
% I do this by interpolating to find the index/cycle number (y) of the
% block nearest to the CO2 blocks in datetime (x). I then use this index to
% find the datetime of these blocks.
% I have to use the index as the interpolant as the 'nearest' neighbour is
% calculated in the x variable and I want to interpolate to the nearest
% neightbour in datetime, not in cycle number.

iCO2 = importedData.Method=="CO2nonB" | importedData.Method=="CO2_non_B"; % Identify the CO2 blocks using the method name
rCO2N2_SA = importedData.x1_CycleInt_Samp_40(iCO2)./importedData.x1_CycleInt_Samp_29(iCO2);
rCO2N2_ST = importedData.x1_CycleInt_Ref_40(iCO2)./importedData.x1_CycleInt_Ref_29(iCO2);

x = importedData.datenum + cumsum(ones(size(importedData.datenum))).*importedData.datenum*eps; % Generate an x-variable with unique points by adding a quasi-infinitesimal increment (eps) to each row
y = 1:length(x);

idxCO2NearestBlocks = interp1(x(~iCO2),y(~iCO2),x(iCO2),'nearest'); % Interpolate to find the index of the nearest blocks
datetimeCO2NearestBlocks = importedData.datetime(idxCO2NearestBlocks); % Use the index to find the datetime of those blocks

[iA,idxB]=ismember(datetimeCO2NearestBlocks,aliquot_metadata.msDatetime); % Find the linear index of these blocks in the aliquot_metadata structure, not all will be present as some have <16 cycles or 6 < blocks < 4
[idxCO2Aliquots,idx2,idx3] = ind2sub(size(aliquot_metadata.msDatetime),idxB(iA)); % Convert to a subscript so that I know which aliquot they correspond to

aliquot_metadata.rCO2N2_SA = nan(size(aliquot_metadata.msDatetime,1),1);
aliquot_metadata.rCO2N2_ST = nan(size(aliquot_metadata.msDatetime,1),1);
aliquot_metadata.dCO2N2 = nan(size(aliquot_metadata.msDatetime,1),1);

aliquot_metadata.rCO2N2_SA(idxCO2Aliquots) = rCO2N2_SA(iA);
aliquot_metadata.rCO2N2_ST(idxCO2Aliquots) = rCO2N2_ST(iA);
aliquot_metadata.dCO2N2(idxCO2Aliquots) = (rCO2N2_SA(iA)./rCO2N2_ST(iA)-1)*1000;


%% Do Some Staistics on the Blocks
% There are some weird looking blocks here that plot off the top of the
% y-axis. I should take a closer look. Set a threshold for inclusion?

figure;
subplot(211); hold on;
plot(cycle_metadata.msDatetime(idx_blocksAll),blockLengthsAll,'-','Color',lineCol(9),'LineWidth',1)
plot(cycle_metadata.msDatetime(idx_blocks),blockLengths,'-','Color',lineCol(10),'LineWidth',3)
ylabel('Number of Cycles in Block');
title('Block Length (by method)');
ylim([0 20])

subplot(212); hold on;
plot(block_metadata.msDatetime(idx_SampleAliquotsAll,1),aliquotLengthsAll,'-','Color',lineCol(9),'LineWidth',1);
plot(block_metadata.msDatetime(idx_SampleAliquots,1),aliquotLengths,'-','Color',lineCol(10),'LineWidth',3);
ylabel('Number of Blocks in Aliquot');
title('Aliquot Lengths (by method)')
ylim([0 10])

figure
subplot(121); hold on;
histogram(blockLengthsAll(blockLengthsAll~=16))
text(15,75,{['No. Total: ' num2str(length(blockLengthsAll))]; ...
    ['No. Included: ' num2str(numberOfBlocks)]; ...
    ['No. Rejected: ' num2str(length(blockLengthsAll(blockLengthsAll~=16)))]; ...
    ['Rejected: ' num2str(length(blockLengthsAll(blockLengthsAll~=16))/length(blockLengthsAll)*100) '%']}, ...
    'HorizontalAlignment','right','VerticalAlignment','Top');
xlabel('Block Length (# cycles)');
ylabel('Counts')

subplot(122); hold on;
histogram(aliquotLengthsAll(aliquotLengthsAll < 4 | aliquotLengthsAll > 5))
text(9,17,{['No. Total: ' num2str(length(aliquotLengthsAll))]; ...
    ['No. Included: ' num2str(numberOfAliquots)]; ...
    ['No. Rejected: ' num2str(length(aliquotLengthsAll(aliquotLengthsAll < 4 | aliquotLengthsAll > 5)))]; ...
    ['Rejected: ' num2str(length(aliquotLengthsAll(aliquotLengthsAll < 4 | aliquotLengthsAll > 5))/length(aliquotLengthsAll)*100) '%']}, ...
    'HorizontalAlignment','right','VerticalAlignment','Top');
xlabel('Aliquot Length (# blocks)');
ylabel('Counts');

suptitle('Properties of Rejected Blocks and Aliquots')

figure
subplot(211)
semilogy(block_metadata.msDatetime(:,1),std(block_deltas(:,:,:),0,3),'.');
ylabel('Std Dev of Cycles in a Block [per mil]');
ylim([0 150]);
legend(delta_cols,'Location','N','Orientation','Horizontal');

subplot(212)
semilogy(aliquot_metadata.msDatetime(:,1,1),nanstd(mean(aliquot_deltas(:,:,:,:),3),0,4),'.');
ylabel('Std Dev of Blocks in an Aliquot [per mil]')
ylim([0 1500])


%% Calculate the PIS for each PIS experiment
% Identifies the aliquots containing PIS experiments and calculates the
% pressure imbalance sensitivity of each delta vlaue for each experiment.

iPIS = all(~isnan(mean(aliquot_deltas,4)),3); % identify aliquots where the mean of ANY of the delta values is not NaN for ALL the blocks
iPIS = iPIS & aliquot_metadata.ID1(:,5,1)=='PIS'; % limit the selection to just those identified as a PIS experiment by their Sample ID1

if sum(any(iPIS,2)) ~= sum(all(iPIS,2)) % Check that all delta values identify each PIS experiment
    warning('Warning: Some delta values are misisng a PIS block for one or more experiments')
    iPIS = any(iPIS,2); % If some delta values are missing a PIS block somehow, calculate the PIS for the other delta values anyway
else
    iPIS = iPIS(:,1); % Otherwise, just make iPis a vector (from a matrix) by taking the first column.
end


% Calculate the PIS for each experiment for each of the delta values
calcPis = nan(size(aliquot_deltas,[1 2]));
calcPisRsq = nan(size(aliquot_deltas,[1 2]));
calcPisPval = nan(size(aliquot_deltas,[1 2]));
calcPisImbal = nan(size(aliquot_deltas,[1 2]));

for ii=find(iPIS)' % find the indices of the PIS aliquots and loop through them
    for jj=1:numel(delta_cols) % loop all through delta values, skip the first three columns as these are voltages and pressure imbalance
        
        d = squeeze(nanmean(aliquot_deltas(ii,jj,:,:),4)); % response variable = the looped delta value from the looped aliquot
        G = [ones(size(d)) nanmean(aliquot_metadata.pressureImbal(ii,:,:),3)']; % predictor variable = the pressure imbalance (col 3) from the looped variable
        m = (G'*G)\G'*d; % Calculate the PIS
        
        [R_corr,pVal] = corrcoef(G(:,2),d); % Find the r-squared correlation coefficient for the PIS test
        [pImbal, idx] = max(abs(G(:,2))); % Find the block with the max P Imbalance
        
        if idx ~= 5
            warning(['For the PIS experiment on ' datestr(aliquot_metadata.msDatenum(ii,1,1),'dd-mmm-yyyy HH:MM') ' the block with the largest imbalance is block ' num2str(idx) ', not block 5.'])
        end
        
        calcPis(ii,jj)=m(2);
        calcPisRsq(ii,jj) = R_corr(1,2).^2;
        calcPisPval(ii,jj) = pVal(1,2);
        calcPisImbal(ii,jj) = pImbal * sign(G(idx,2));
    end
end


%% Identify Anomalous PIS values

% Reject all PIS values where the P Imbalance is smaller than 100 mV
iPisRejections = abs(calcPisImbal)<100;
calcPis(iPisRejections) = nan;
calcPisRsq(iPisRejections) = nan;
calcPisPval(iPisRejections) = nan;
calcPisImbal(iPisRejections) = nan;

%pisRejectionMethodTesting; % Run to test different rejection criteria

% Identify PIS values that differ significantly from the running median for a given delta value
numRej = zeros(1,numel(delta_cols));
movWindow=49; % Moving median window = 7 weeks


for ii = 1:numel(delta_cols)
    x = aliquot_metadata.msDatenum(iPIS,1,1);
    y = calcPis(iPIS,ii,1,1);

    cen = movmedian(y,movWindow,'omitnan','SamplePoints',x); % Calculate moving median for given window width
    dev = y-cen; % Detrend time-series by calculating deviation of each point from moving median

    [CDF,edges] = histcounts(dev,'BinMethod','fd','Normalization','cdf'); % Calculate CDF of the deviations, using Freedman-Diaconis rule for bin widths

    low = cen + edges(find(CDF>0.01,1,'first')); % Lower Bound = Median - First Percentile Deviation
    upp = cen + edges(find(CDF>0.99,1,'first')); % Upper Bound = Median + Ninety Ninth Percentile Deviation
    iRej = (y > upp) | (y < low);

    iPisRejections(iPIS,ii) = iPisRejections(iPIS,ii) | iRej; % assign the rejections back to the full-size variable
    numRej(ii) = sum(iRej); % tally the number of rejections for each combination of delta value and window size

end


%% Plot a time-series of the PIS and related parameters
% This is useful to identify aliquots where the PIS block did no run
% correctly or where the r-squared is low, suggesting a poor determination
% of the PIS. These cases can be filtered out below.

stackedFig(3,'RelSize',[0.4 1.7 0.9],'Overlap',[-10 -10]);
stackedFigAx
xlim(datenum(['01 Jan 2016'; '31 Dec 2018']))

% Plot R-Squared of each PIS Experiment
stackedFigAx(1)
for ii = 1:numel(delta_cols)
    set(gca,'ColorOrderIndex',ii)
    plot(aliquot_metadata.msDatenum(iPIS,1,1),calcPisRsq(iPIS,ii),'-ok','MarkerIndices',find(iPisRejections(iPIS,ii)),'MarkerFaceColor',lineCol(ii));
    legH(ii)=plot(aliquot_metadata.msDatenum(~iPisRejections(:,ii) & iPIS,1,1),calcPisRsq(~iPisRejections(:,ii) & iPIS,ii),'s-','MarkerFaceColor',lineCol(ii));
end
legend(legH,delta_cols,'Orientation','Horizontal','Location','South')
ylabel('r^2');
ylim([0.7 1]);

% Plot PIS Value for each PIS Experiment
stackedFigAx(2)
for ii = 1:numel(delta_cols)
    set(gca,'ColorOrderIndex',ii)
    plot(aliquot_metadata.msDatenum(iPIS,1,1),calcPis(iPIS,ii),'-ok','MarkerIndices',find(iPisRejections(iPIS,ii)),'MarkerFaceColor',lineCol(ii));
    legH(ii)=plot(aliquot_metadata.msDatenum(~iPisRejections(:,ii) & iPIS,1,1),calcPis(~iPisRejections(:,ii) & iPIS,ii),'s-','MarkerFaceColor',lineCol(ii));
end
ylabel('PIS [per mil/per mil]');
ylim([-0.01 0.005]);

% Plot Pressure Imbalance for each PIS Experiment
stackedFigAx(3)
plot(aliquot_metadata.msDatenum(iPIS & any(~iPisRejections,2),1,1),calcPisImbal(iPIS & any(~iPisRejections,2),:),'-^','Color',lineCol(9))
text(aliquot_metadata.msDatenum(iPIS & any(~iPisRejections,2),1,1),calcPisImbal(iPIS & any(~iPisRejections,2)),aliquot_metadata.ID1(iPIS & any(~iPisRejections,2),5,1))
ylabel('Pressure Imbalance [per mil]')
ylim([-600 0])

stackedFigAx();
datetick('x');
xlim(datenum(['01 Jan 2016'; '31 Dec 2018']))

stackedFigReset


%% Remove PIS Blocks
% Now remove the fifth blocks from the array of aliquot delta values so
% they don't get used in with further analysis
aliquot_deltasPisExp = aliquot_deltas(:,:,5,:);
aliquot_deltasPisExp(~iPIS,:,:,:)=nan;
aliquot_deltas(:,:,5,:) = [];

for ii = 1:numel(metadata_fields)
    aliquot_metadataPisExp.(metadata_fields{ii})(iPIS,:,:) = aliquot_metadata.(metadata_fields{ii})(iPIS,5,:);
    aliquot_metadata.(metadata_fields{ii})(:,5,:) = [];
end


%% Make PIS Correction
PIS = calcPis;
PIS(iPisRejections) = nan;
PIS = repmat(PIS,[1 1 size(aliquot_deltas,[3 4])]);
PIS = fillmissing(PIS,'previous',1);

stackedFig(numel(delta_cols))
for ii=1:numel(delta_cols)
    stackedFigAx(ii)
    plot(aliquot_metadata.msDatenum(~iPisRejections(:,ii),1,1),calcPis(~iPisRejections(:,ii),ii),'o','Color','none','MarkerFaceColor',lineCol(ii))
    plot(aliquot_metadata.msDatenum(:,1,1),PIS(:,ii,1,1),'.','Color',lineCol(ii)*0.5)
    ylabel(delta_cols{ii})
end
stackedFigAx
title('PIS Values used for Correction')
xlabel('Date')
xlim(datenum(["01-Jan-2016" "01-Jan-2019"]))
datetick('x','keeplimits')
stackedFigReset

aliquot_deltas_pisCorr = aliquot_deltas - permute(aliquot_metadata.pressureImbal,[1 4 2 3]).*PIS;

%% Calculate the Chemical Slopes
% Still need to:
%   1) Weed out the questionable blocks/aliquots that mess up some of the
%      chem slope experiments
%   2) Figure out how I'm going to do the CS Correction for Ar isotopes

iCS_15N = contains(squeeze(aliquot_metadata.ID1(:,1,1)),'CS') ...
    & (contains(squeeze(aliquot_metadata.ID1(:,1,1)),'15') ...
    | contains(squeeze(aliquot_metadata.ID1(:,1,1)),'N'));
iCS_18O = contains(squeeze(aliquot_metadata.ID1(:,1,1)),'CS') ...
    & (contains(squeeze(aliquot_metadata.ID1(:,1,1)),'18') ...
    | contains(squeeze(aliquot_metadata.ID1(:,1,1)),'O'));

iCS_Ar = iCS_15N | iCS_18O; % Create joint CS index for d40/36Ar CS experiments

dates = {aliquot_metadata.msDatetime(iCS_15N,1,1) aliquot_metadata.msDatetime(iCS_18O,1,1)}; % just for reference

% Find the different CS experiments by finding the CS aliquots separated by more than 10 hours
diffsCS_15N = duration(nan(3223,3));
diffsCS_18O = diffsCS_15N;
diffsCS_Ar = diffsCS_15N;

diffsCS_15N(iCS_15N) = [diff(aliquot_metadata.msDatetime(iCS_15N,1,1)); hours(999)+minutes(59)+seconds(59)];
diffsCS_18O(iCS_18O) = [diff(aliquot_metadata.msDatetime(iCS_18O,1,1)); hours(999)+minutes(59)+seconds(59)];
diffsCS_Ar(iCS_Ar) = [diff(aliquot_metadata.msDatetime(iCS_Ar,1,1)); hours(999)+minutes(59)+seconds(59)];

% N.B. It's important to not use too big a number here. Using 24 hours
% fails to resolve a re-do of the 18O CS in Feb-2016 as it was run the
% morning after the previous attempt was run in the afternoon w/ diff=15 hr
endCS_15N = diffsCS_15N > 10/24;
endCS_18O = diffsCS_18O > 10/24;
endCS_Ar = diffsCS_Ar > 24; % Use 24 hours for Ar to make sure the re-do is lumped together with the right 15N CS experiment

% Increment an index by one each time its a new CS experiment
idxCS_15NEnd = nan(size(aliquot_deltas_pisCorr,1),1);
idxCS_15NEnd(endCS_15N) = cumsum(endCS_15N(endCS_15N));

idxCS_18OEnd = nan(size(aliquot_deltas_pisCorr,1),1);
idxCS_18OEnd(endCS_18O) = cumsum(endCS_18O(endCS_18O));

idxCS_ArEnd = nan(size(aliquot_deltas_pisCorr,1),1);
idxCS_ArEnd(endCS_Ar) = cumsum(endCS_Ar(endCS_Ar));

% Assign the same index to all the aliquots from each CS experiment
idxCS_15N = zeros(size(aliquot_deltas_pisCorr,1),1); % Create a vector of zeros
idxCS_15N(iCS_15N)=nan; % Assign NaN to the values I want to replace
idxCS_15N(endCS_15N) = idxCS_15NEnd(endCS_15N); % Fill in the indices of the end of each CS experiment
idxCS_15N(iCS_15N) = fillmissing(idxCS_15N(iCS_15N),'next'); % Replace the nans for each aliquot with the next index - must fill only the (iCS_15N) subset otherwise it fills some indices with interspersed zeros that are then reassigned to NaN in the next line
idxCS_15N = standardizeMissing(idxCS_15N,0); % Change the zeros to nans to properly represent the fact that they are missing values

idxCS_18O = zeros(size(aliquot_deltas_pisCorr,1),1); % Create a vector of zeros
idxCS_18O(iCS_18O)=nan; % Assign NaN to the values I want to replace
idxCS_18O(endCS_18O) = idxCS_18OEnd(endCS_18O); % Fill in the indices of the end of each CS experiment
idxCS_18O(iCS_18O) = fillmissing(idxCS_18O(iCS_18O),'next'); % Replace the nans for each aliquot with the next index - must fill only the (iCS_15N) subset otherwise it fills some indices with interspersed zeros that are then reassigned to NaN in the next line
idxCS_18O = standardizeMissing(idxCS_18O,0); % Change the zeros to nans to properly represent the fact that they are missing values

idxCS_Ar = zeros(size(aliquot_deltas_pisCorr,1),1); % Create a vector of zeros
idxCS_Ar(iCS_Ar)=nan; % Assign NaN to the values I want to replace
idxCS_Ar(endCS_Ar) = idxCS_ArEnd(endCS_Ar); % Fill in the indices of the end of each CS experiment
idxCS_Ar(iCS_Ar) = fillmissing(idxCS_Ar(iCS_Ar),'next'); % Replace the nans for each aliquot with the next index - must fill only the (iCS_15N) subset otherwise it fills some indices with interspersed zeros that are then reassigned to NaN in the next line
idxCS_Ar = standardizeMissing(idxCS_Ar,0); % Change the zeros to nans to properly represent the fact that they are missing values


fitCS = @(x,y)(([x ones(length(x),1)]\y)'); % Transpose output to a row vector so that split apply can concatanate the outputs from each group
rsqCS = @(x,y){(corrcoef(x,y)).^2};

% dO2/N2 Effect on d15N
x = squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:),4),3));
y = squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='d15N',:,:),4),3));
P_15N = splitapply(fitCS,x,y,idxCS_15N); % Calculate the slope with the anonymous function
rsq_15N = splitapply(rsqCS,x,y,idxCS_15N);
CS_15N = nan(size(aliquot_deltas_pisCorr(:,1,:,:)));
CS_15N(endCS_15N,:,:,:) = repmat(P_15N(:,1), [1 1 size(CS_15N,[3 4])]); % Assign the slope values to the end aliquots in the time-series of CS values

figure
for ii=1:sum(endCS_15N)
    subplot(1,sum(endCS_15N),ii); hold on
    plot(x(idxCS_15N==ii),y(idxCS_15N==ii),'xk') % Plot the individual aliquots
    plot(x(idxCS_15N==ii),polyval(P_15N(ii,:),x(idxCS_15N==ii)),'-r') % Plot the fitted line
    text(max(x(idxCS_15N==ii)),min(y(idxCS_15N==ii)),compose('CS = %.2f per meg/per mil\nr^2 = %.4f',P_15N(ii,1)*1000,rsq_15N{ii}(1,2)),'HorizontalAlignment','Right','VerticalAlignment','bottom')
    axis([-10 350 -0.01 0.25]);
    xlabel('\deltaO_2/N_2 [per mil]');
    ylabel('\delta^{15}N [per mil]');
    title(['\delta^{15}N CS: ' datestr(aliquot_metadata.msDatetime(idxCS_15NEnd==ii,1,1),'yyyy-mmm-dd')])
end

% dO2/N2 Effect on dArN2
x = squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:),4),3));
y = squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dArN2',:,:),4),3));
P_ArN2 = splitapply(fitCS,x,y,idxCS_15N); % Calculate the slope with the anonymous function
rsq_ArN2 = splitapply(rsqCS,x,y,idxCS_15N);
CS_ArN2 = nan(size(aliquot_deltas_pisCorr(:,1,:,:)));
CS_ArN2(endCS_15N,:,:,:) = repmat(P_ArN2(:,1), [1 1 size(CS_ArN2,[3 4])]); % Assign the slope values to the end aliquots in the time-series of CS values

figure
for ii=1:sum(endCS_15N)
    subplot(1,sum(endCS_15N),ii); hold on
    plot(x(idxCS_15N==ii),y(idxCS_15N==ii),'xk') % Plot the individual aliquots
    plot(x(idxCS_15N==ii),polyval(P_ArN2(ii,:),x(idxCS_15N==ii)),'-r') % Plot the fitted line
    text(max(x(idxCS_15N==ii)),min(y(idxCS_15N==ii)),compose('CS = %.2f per meg/per mil\nr^2 = %.4f',P_ArN2(ii,1)*1000,rsq_ArN2{ii}(1,2)),'HorizontalAlignment','Right','VerticalAlignment','bottom')
    axis([-10 350 -0.1 1.2]);
    xlabel('\deltaO_2/N_2 [per mil]');
    ylabel('\deltaAr/N_2 [per mil]');
    title(['\deltaAr/N_2 CS: ' datestr(aliquot_metadata.msDatetime(idxCS_15NEnd==ii,1,1),'yyyy-mmm-dd')])
end

% dN2/O2 Effect on d18O
x = ((squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:),4),3))/1000+1).^-1-1)*1000;
y = squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='d18O',:,:),4),3));
P_18O = splitapply(fitCS,x,y,idxCS_18O); % Calculate the slope with the anonymous function
rsq_18O = splitapply(rsqCS,x,y,idxCS_18O);
CS_18O = nan(size(aliquot_deltas_pisCorr(:,1,:,:)));
CS_18O(endCS_18O,:,:,:) = repmat(P_18O(:,1), [1 1 size(CS_18O,[3 4])]); % Assign the slope values to the end aliquots in the time-series of CS values

figure
for ii=1:sum(endCS_18O)
    subplot(1,sum(endCS_18O),ii); hold on
    plot(x(idxCS_18O==ii),y(idxCS_18O==ii),'xk') % Plot the individual aliquots
    plot(x(idxCS_18O==ii),polyval(P_18O(ii,:),x(idxCS_18O==ii)),'-r') % Plot the fitted line
    text(max(x(idxCS_18O==ii)),min(y(idxCS_18O==ii)),compose('CS = %.2f per meg/per mil\nr^2 = %.4f',P_18O(ii,1)*1000,rsq_18O{ii}(1,2)),'HorizontalAlignment','Right','VerticalAlignment','bottom')
    axis([-10 350 -0.1 0.1]);
    xlabel('\deltaN_2/O_2 [per mil]');
    ylabel('\delta^{18}O [per mil]');
    title(['\delta^{18}O CS: ' datestr(aliquot_metadata.msDatetime(idxCS_18OEnd==ii,1,1),'yyyy-mmm-dd')])
end

% dN2/O2 Effect on d17O
x = ((squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:),4),3))/1000+1).^-1-1)*1000;
y = squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='d17O',:,:),4),3));
P_17O = splitapply(fitCS,x,y,idxCS_18O); % Calculate the slope with the anonymous function
rsq_17O = splitapply(rsqCS,x,y,idxCS_18O);
CS_17O = nan(size(aliquot_deltas_pisCorr(:,1,:,:)));
CS_17O(endCS_18O,:,:,:) = repmat(P_17O(:,1), [1 1 size(CS_17O,[3 4])]); % Assign the slope values to the end aliquots in the time-series of CS values

figure
for ii=1:sum(endCS_18O)
    subplot(1,sum(endCS_18O),ii); hold on
    plot(x(idxCS_18O==ii),y(idxCS_18O==ii),'xk') % Plot the individual aliquots
    plot(x(idxCS_18O==ii),polyval(P_17O(ii,:),x(idxCS_18O==ii)),'-r') % Plot the fitted line
    text(max(x(idxCS_18O==ii)),min(y(idxCS_18O==ii)),compose('CS = %.2f per meg/per mil\nr^2 = %.4f',P_17O(ii,1)*1000,rsq_17O{ii}(1,2)),'HorizontalAlignment','Right','VerticalAlignment','bottom')
    axis([-10 350 -0.1 1]);
    xlabel('\deltaN_2/O_2 [per mil]');
    ylabel('\delta^{17}O [per mil]');
    title(['\delta^{17}O CS: ' datestr(aliquot_metadata.msDatetime(idxCS_18OEnd==ii,1,1),'yyyy-mmm-dd')])
end

% d4036Ar
% N.B. This one is a little different as there are two chemical slopes, one
% for the effect of the N2/Ar ratio and one for the isobaric interference
% of 18O18O, i.e. the O2/Ar ratio.
x = [((squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dArN2',:,:),4),3))./1000+1).^-1-1)*1000 ((squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:),4),3))/1000+1)./((squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dArN2',:,:),4),3))./1000+1))-1)*1000]; % predictor variables = dN2/Ar AND dO2Ar (= [q_o2n2/q_arn2 -1]*1000)
y = squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='d4036Ar',:,:),4),3));
P_4036Ar = splitapply(fitCS,x,y,idxCS_Ar); % Calculate the slope with the anonymous function
CS_4036Ar = nan(size(aliquot_deltas_pisCorr(:,1:2,:,:)));
CS_4036Ar(endCS_Ar,:,:,:) = repmat(P_4036Ar(:,1:2), [1 1 size(CS_4036Ar,[3 4])]); % Assign the slope values to the end aliquots in the time-series of CS values

figure
for ii=1:sum(endCS_Ar)
    subplot(1,sum(endCS_Ar),ii); hold on
    plot3(x(idxCS_Ar==ii,1),x(idxCS_Ar==ii,2),y(idxCS_Ar==ii),'xk'); % Plot the individual aliquots
    [X1,X2]=meshgrid(linspace(min(x(idxCS_Ar==ii,1)),max(x(idxCS_Ar==ii,1)),20),linspace(min(x(idxCS_Ar==ii,2)),max(x(idxCS_Ar==ii,2)),20)); % Create a regularly spaced grid
    surf(X1,X2,X1.*P_4036Ar(ii,1)+X2*P_4036Ar(ii,2)+P_4036Ar(ii,1));  % Plot the fitted surface on the grid
    
    text(max(x(idxCS_Ar==ii,1)),min(x(idxCS_Ar==ii,2)),min(y(idxCS_Ar==ii)),compose('CS N_2/Ar = %.2f per meg/per mil\nCS O_2/Ar = %.2f per meg/per mil',P_4036Ar(ii,1)*1000,P_4036Ar(ii,2)*1000),'HorizontalAlignment','Right','VerticalAlignment','bottom');
    
    xlabel('\deltaN_2/Ar [per mil]')
    ylabel('\deltaO_2/Ar [per mil]')
    zlabel('\delta^{40}/_{36}Ar [per mil]')
    title(['\delta^{40}/_{36}Ar CS: ' datestr(aliquot_metadata.msDatetime(idxCS_ArEnd==ii,1,1),'yyyy-mmm-dd')])
    
    colormap(cbrewer('seq','Greens',20));
    axis([-20 350 -20 350 -8 1]);
    caxis([-8 0]);
    view([40 45]);
end
pos=get(gca,'Position');
colorbar;
set(gca,'Position',pos);

% d4038Ar
% N.B. This one is a little different as there are two chemical slopes, one
% for the effect of the N2/Ar ratio and one for the isobaric interference
% of 18O18O, i.e. the O2/Ar ratio.
x = [((squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dArN2',:,:),4),3))./1000+1).^-1-1)*1000 ((squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:),4),3))/1000+1)./((squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='dArN2',:,:),4),3))./1000+1))-1)*1000]; % predictor variables = dN2/Ar AND dO2Ar (= [q_o2n2/q_arn2 -1]*1000)
y = squeeze(nanmean(mean(aliquot_deltas_pisCorr(:,delta_cols=='d4038Ar',:,:),4),3));
P_4038Ar = splitapply(fitCS,x,y,idxCS_Ar); % Calculate the slope with the anonymous function
CS_4038Ar = nan(size(aliquot_deltas_pisCorr(:,1:2,:,:)));
CS_4038Ar(endCS_Ar,:,:,:) = repmat(P_4036Ar(:,1:2), [1 1 size(CS_4038Ar,[3 4])]); % Assign the slope values to the end aliquots in the time-series of CS values

figure
for ii=1:sum(endCS_Ar)
    subplot(1,sum(endCS_Ar),ii); hold on
    plot3(x(idxCS_Ar==ii,1),x(idxCS_Ar==ii,2),y(idxCS_Ar==ii),'xk'); % Plot the individual aliquots
    [X1,X2]=meshgrid(linspace(min(x(idxCS_Ar==ii,1)),max(x(idxCS_Ar==ii,1)),20),linspace(min(x(idxCS_Ar==ii,2)),max(x(idxCS_Ar==ii,2)),20)); % Create a regularly spaced grid
    surf(X1,X2,X1.*P_4038Ar(ii,1)+X2*P_4038Ar(ii,2)+P_4038Ar(ii,1));  % Plot the fitted surface on the grid
    
    text(max(x(idxCS_Ar==ii,1)),min(x(idxCS_Ar==ii,2)),min(y(idxCS_Ar==ii)),compose('CS N_2/Ar = %.2f per meg/per mil\nCS O_2/Ar = %.2f per meg/per mil',P_4038Ar(ii,1)*1000,P_4038Ar(ii,2)*1000),'HorizontalAlignment','Right','VerticalAlignment','bottom');
    
    xlabel('\deltaN_2/Ar [per mil]')
    ylabel('\deltaO_2/Ar [per mil]')
    zlabel('\delta^{40}/_{38}Ar [per mil]')
    title(['\delta^{40}/_{38}Ar CS: ' datestr(aliquot_metadata.msDatetime(idxCS_ArEnd==ii,1,1),'yyyy-mmm-dd')])
    
    colormap(cbrewer('seq','Greens',20));
    axis([-50 350 -50 350 -16 1]);
    caxis([-8 0]);
    view([40 45]);
end
pos=get(gca,'Position');
colorbar;
set(gca,'Position',pos);


%% Make the CS Corrections

CS = [CS_15N CS_18O CS_17O zeros(size(CS_15N)) zeros(size(CS_15N)) zeros(size(CS_15N)) CS_ArN2];
CS = fillmissing(CS,'previous',1);

CS_predictors = [aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:) ... % O2N2 CS on d15N
    ((aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:)/1000+1).^-1-1)*1000 ... % N2O2 CS on d18O
    ((aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:)/1000+1).^-1-1)*1000 ... % N2O2 CS on d17O
    zeros(size(aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:))) ... % d4036Ar CS Below
    zeros(size(aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:))) ... % d4038Ar CS Below
    zeros(size(aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:))) ... % No CS Corr for dO2N2
    aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:)]; % O2N2 CS on dArN2

aliquot_deltas_pisCorr_csCorr = aliquot_deltas_pisCorr - CS.*CS_predictors;

CS_4036Ar = fillmissing(CS_4036Ar,'previous',1);
aliquot_deltas_pisCorr_csCorr(:,delta_cols=='d4036Ar',:,:) = aliquot_deltas_pisCorr_csCorr(:,delta_cols=='d4036Ar',:,:) - (CS_4036Ar(:,1,:,:).*((aliquot_deltas_pisCorr(:,delta_cols=='dArN2',:,:)/1000+1).^-1-1)*1000) - (CS_4036Ar(:,2,:,:).*((aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:)/1000+1)./(aliquot_deltas_pisCorr(:,delta_cols=='dArN2',:,:)/1000+1)-1)*1000);

CS_4038Ar = fillmissing(CS_4038Ar,'previous',1);
aliquot_deltas_pisCorr_csCorr(:,delta_cols=='d4038Ar',:,:) = aliquot_deltas_pisCorr_csCorr(:,delta_cols=='d4038Ar',:,:) - (CS_4038Ar(:,1,:,:).*((aliquot_deltas_pisCorr(:,delta_cols=='dArN2',:,:)/1000+1).^-1-1)*1000) - (CS_4038Ar(:,2,:,:).*((aliquot_deltas_pisCorr(:,delta_cols=='dO2N2',:,:)/1000+1)./(aliquot_deltas_pisCorr(:,delta_cols=='dArN2',:,:)/1000+1)-1)*1000);


%% Calculate the LJA Normalization Values
% Identified the aliquots that correspond to the different batches of LJA
% measurements. Plots the distribution (box plots) of aliquot means for
% each batch and identifies and removes any outliers. Calculates the mean
% of aliquot means for each batch to use for LJA normalization.

% Identify the LJA aliquots
iLja = contains(aliquot_metadata.ID1(:,1,1),'LJA');

% Identify each batch of LJA measurements as separated in time by >5 days
diffsLJA = duration(nan(size(aliquot_deltas_pisCorr_csCorr,1),3));
diffsLJA(iLja) = [diff(aliquot_metadata.msDatetime(iLja,1,1)); hours(999)+minutes(59)+seconds(59)];

iLjaEnd = diffsLJA > 5;

% Increment an index by one each time its a new batch of LJA measurements
idxLjaAliquotsEnd = nan(size(aliquot_deltas_pisCorr_csCorr,1),1);
idxLjaAliquotsEnd(iLjaEnd) = cumsum(iLjaEnd(iLjaEnd));

% Assign the same index to all the aliquots from each given batch
idxLjaAliquots = zeros(size(aliquot_deltas_pisCorr_csCorr,1),1); % Create a vector of zeros
idxLjaAliquots(iLja)=nan; % Assign NaN to the values I want to replace
idxLjaAliquots(iLjaEnd) = idxLjaAliquotsEnd(iLjaEnd); % Fill in the indices of the end of each LJA aliquot
idxLjaAliquots(iLja) = fillmissing(idxLjaAliquots(iLja),'next'); % Replace the nans for each aliquot with the next index - must fill only the (iLJA) subset otherwise it fills some indices with interspersed zeros that are then reassigned to NaN in the next line
idxLjaAliquots = standardizeMissing(idxLjaAliquots,0); % Change the zeros to nans to properly represent the fact that they are missing values

% Plot the distribution of all LJA aliquot means
figure
for ii = 1:numel(delta_cols)
    subplot(1,numel(delta_cols),ii)
    histogram(nanmean(mean(aliquot_deltas_pisCorr_csCorr(iLja,ii,:,:),4),3))
    axis('square');
    xlabel('\delta [per mil]'); ylabel('Counts')
    title(delta_cols(ii))
end

% Calculate stats for each batch and make box plots
N = nan(sum(iLjaEnd),1);
stdev = nan(sum(iLjaEnd),numel(delta_cols));
for ii=1:sum(iLjaEnd)
    N(ii) = sum(idxLjaAliquots==ii); % number of aliquots measured for each batch of LJA measurements
    stdev(ii,:) = nanstd(nanmean(mean(aliquot_deltas_pisCorr_csCorr(idxLjaAliquots==ii,:,:,:),4),3)); % standard deviation of aliquot means for each batch of LJA measurements
end
SEM = stdev./sqrt(N); % standard error of aliquot means for each batch of LJA measurements

% Make Box Plots
labels = datestr(aliquot_metadata.msDatenum(iLjaEnd,1,1),'dd-mmm-yyyy');
for ii = 1:numel(delta_cols)
    figure; hold on;
    hBoxPlot=boxplot(nanmean(mean(aliquot_deltas_pisCorr_csCorr(iLja,ii,:,:),4),3),idxLjaAliquots(iLja),'Labels',labels,'Notch','off');
    plot(idxLjaAliquots(iLja),nanmean(mean(aliquot_deltas_pisCorr_csCorr(iLja,ii,:,:),4),3),'.k');
    text(1:sum(iLjaEnd),repmat(min(nanmean(mean(aliquot_deltas_pisCorr_csCorr(iLja,ii,:,:),4),3))*1.05,1,sum(iLjaEnd)),compose(['N = %d\n\\sigma = %.3f' char(8240) '\nSEM = %.3f' char(8240)],N,stdev(:,ii),SEM(:,ii)));
    
    ylim('auto');
    ylabel('\delta_{LJA} [per mil]')
    title(['LJA ' delta_cols(ii)])
end

% Reject Outliers
isoutlierAnonFn = @(x){(isoutlier(x,'quartiles'))};
tf=splitapply(isoutlierAnonFn,nanmean(mean(aliquot_deltas_pisCorr_csCorr(iLja,:,:,:),4),3),idxLjaAliquots(iLja));

iLjaAfterRej = repmat(iLja,[1 numel(delta_cols)]);
for ii=1:length(tf)
    iLjaAfterRej(idxLjaAliquots==ii,:) = ~tf{ii};
end

% Calculate LJA Values for Normalization
ljaValues = nan(size(aliquot_deltas_pisCorr_csCorr));
ljaLoCI = nan(size(aliquot_deltas_pisCorr_csCorr));
ljaHiCI = nan(size(aliquot_deltas_pisCorr_csCorr));

for ii=1:7
    [batchMean, batchMeanCI] = grpstats(nanmean(mean(aliquot_deltas_pisCorr_csCorr(iLjaAfterRej(:,ii),ii,:,:),4),3),idxLjaAliquots(iLjaAfterRej(:,ii)),{'mean','meanci'});
    ljaValues(iLjaEnd,ii,:,:) = repmat(batchMean,[1 1 size(aliquot_deltas_pisCorr_csCorr,[3 4])]);
    ljaLoCI(iLjaEnd,ii,:,:) = repmat(batchMeanCI(:,1),[1 1 size(aliquot_deltas_pisCorr_csCorr,[3 4])]);
    ljaHiCI(iLjaEnd,ii,:,:) = repmat(batchMeanCI(:,2),[1 1 size(aliquot_deltas_pisCorr_csCorr,[3 4])]);
end


%% Make the LJA Correction

LJA = ljaValues;
LJA = fillmissing(LJA,'previous',1);

aliquot_deltas_pisCorr_csCorr_ljaCorr = ((aliquot_deltas_pisCorr_csCorr/1000+1)./(LJA/1000+1)-1)*1000;


%%
set(0,'defaultFigureVisible','on');
disp('>> Script Complete: Turning Figures Back On');