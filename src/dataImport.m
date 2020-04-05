%% dataImport %%

%% Import Data and Fix Variable Types
% Uncomment this section when running the script for the first time in a
% session. It's easiest to then comment these lines back out and not clear
% the xp2018 etc. variables as they take a long time to load in.

% clearvars; clc;
% disp({'Loading file 1: Working...'})
% xp2018 = readtable('XP-2018(excelExportIntensityJDM).csv','Delimiter',',');
% disp({'Loading file 1: Complete'}); disp({'Loading file 2: Working...'})
% xp2017 = readtable('XP-2017(excelExportIntensityJDM).csv','Delimiter',',');
% disp({'Loading file 2: Complete'}); disp({'Loading file 3: Working...'})
% xp2016 = readtable('XP-2016(excelExportIntensityJDM).csv','Delimiter',',');
% disp({'Loading file 3: Complete'}); disp({'Loading file 4: Working...'})
% % xp2015 = readtable('XP-2015(excelExportIntensityJDM).csv','Delimiter',',');
% % disp({'Loading file 4: Complete'})

%%
xp2018.MeasurmentErrors = num2cell(xp2018.MeasurmentErrors);
xp2017.MeasurmentErrors = num2cell(xp2017.MeasurmentErrors);
xp2016.MeasurmentErrors = num2cell(xp2016.MeasurmentErrors);
%xp2015.MeasurmentErrors = num2cell(xp2015.MeasurmentErrors);
importedData = [xp2018; xp2017; xp2016]; % ; xp2015];

importedData.TimeCode = datetime(datenum(importedData.TimeCode),'ConvertFrom','Datenum'); %This is probably redundant...
importedData.Date = datestr(datenum(importedData.Date) + datenum('31 Dec 1999')); %Correct for two-character month '0018'
importedData.IsRef__ = logical(importedData.IsRef__);
importedData.Method = string(importedData.Method);
importedData.Identifier1 = string(importedData.Identifier1);

idxTimeCode = find(cellfun(@(varName) strcmp('TimeCode',varName),importedData.Properties.VariableNames)==1);
importedData = sortrows(importedData,idxTimeCode);

%% Check the Gas Configuration and Gas Names
% Do something here to highlight the different unique values of the Gas
% Config or the Gas Name (CAREFUL! These can be different) and ask the user
% to select which ones they want to include? For example, discard
% everything that is not 'Air+' or discard everything that is 'Ne'.
%
% Could do something similar here with folder names etc.

% 1. Only Use Data from the "SPICE" Folders
% importedData = importedData(contains(importedData.FileHeader_Filename,'SPICE'),:);

% 2. Remove Cycles with No ID1 (Sample ID/Bottom Depth)
% disp(['Removing ' num2str(sum(ismissing(importedData.Identifier1))) ' cycles with no ID1']);
% importedData = importedData(~ismissing(importedData.Identifier1),:);

%% Interpolate ST Voltages and Calculate Delta Values
% NOTE: The current method of interpolating onto the ~isRef indices
% excludes the CO2 checks from the blocks I consider as they are recorded 
% in isRef as 'true'. This is probably the right way to do it for now but I
% must remember to add them back in, including their metadata, which is
% also excluded when I subsample the metadata columns for the ~isRef rows.

intSA = importedData(~importedData.IsRef__,1:9);
intST = array2table(interp1(find(importedData.IsRef__), ...
                    importedData{importedData.IsRef__,1:9}, ...
                    find(~importedData.IsRef__)));
intST.Properties.VariableNames = intSA.Properties.VariableNames;

cycle_deltas = table();
cycle_deltas.int28SA = intSA.rIntensity28;
cycle_deltas.intST28 = intST.rIntensity28;
%cycle_deltas.pressure = (intSA.rIntensity28./intST.rIntensity28 - 1)*1000; % N.B. - Ross typically calculates this just as a raw SA - ST rather than a delta value. Does this make a difference?
cycle_deltas.pressure = (intSA.rIntensity28 - intST.rIntensity28);

cycle_deltas.d15N = ((intSA.rIntensity29./intSA.rIntensity28)./(intST.rIntensity29./intST.rIntensity28) - 1)*1000;
cycle_deltas.d18O = ((intSA.rIntensity34./intSA.rIntensity32)./(intST.rIntensity34./intST.rIntensity32) - 1)*1000;
cycle_deltas.d17O = ((intSA.rIntensity33./intSA.rIntensity32)./(intST.rIntensity33./intST.rIntensity32) - 1)*1000;
cycle_deltas.d4036Ar = ((intSA.rIntensity40./intSA.rIntensity36)./(intST.rIntensity40./intST.rIntensity36) - 1)*1000;
cycle_deltas.d4038Ar = ((intSA.rIntensity40./intSA.rIntensity38)./(intST.rIntensity40./intST.rIntensity38) - 1)*1000;
cycle_deltas.dO2N2 = ((intSA.rIntensity32./intSA.rIntensity28)./(intST.rIntensity32./intST.rIntensity28) - 1)*1000;
cycle_deltas.dArN2 = ((intSA.rIntensity40./intSA.rIntensity28)./(intST.rIntensity40./intST.rIntensity28) - 1)*1000;

delta_cols = cycle_deltas.Properties.VariableNames;
cycle_deltas = table2array(cycle_deltas);

%% Compile Useful Metadata
% NOTE: The current approach of subsampling the colums for the ~isRef rows
% excludes the CO2 check rows. I must remember to add these back in later.
% Also, the cycles all have exactly the same metadata so this step is
% somewhat pointless.

cycle_metadata.msDatetime = importedData.TimeCode(~importedData.IsRef__);
cycle_metadata.filename = importedData.FileHeader_Filename(~importedData.IsRef__);
cycle_metadata.sequenceRow = importedData.Row(~importedData.IsRef__);
cycle_metadata.ASInlet = importedData.AS_SIOInlet(~importedData.IsRef__);
cycle_metadata.ID1 = importedData.Identifier1(~importedData.IsRef__);
cycle_metadata.method = importedData.Method(~importedData.IsRef__);
cycle_metadata.scriptName = importedData.ScriptName(~importedData.IsRef__);
cycle_metadata.gasConfig = importedData.GasConfiguration(~importedData.IsRef__);
cycle_metadata.gasName = importedData.GasName(~importedData.IsRef__);

metadata_fields = fieldnames(cycle_metadata)'; % Transpose it to a row vector so that it works as a loop index

%% Reshape into Cycles-x-Isotope Ratio-x-Block
% Identify the different blocks by the unique filenames
[~,idx_blocks,~] = unique(cycle_metadata.filename,'stable');
blockLengths = diff(idx_blocks); blockLengths(end+1)=length(cycle_deltas)-(idx_blocks(end)-1);

% TAKE ONLY THE BLOCKS WITH 16 CYCLES! - This causes me to lose 156 blocks
% (17672 -> 17516, <1%), mostly with <5 cycles in them.
idx_blocks = idx_blocks(blockLengths==16); 
blockLengths = blockLengths(blockLengths==16);

numberOfBlocks = length(idx_blocks);
longestBlock = max(blockLengths);

% Fill the Array of Delta Values
block_deltas = nan(longestBlock,size(cycle_deltas,2),numberOfBlocks);

for ii = 1:length(idx_blocks)
    block_deltas(1:blockLengths(ii),:,ii) = cycle_deltas(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    for jj = 1:numel(metadata_fields)
        block_metadata.(metadata_fields{jj})(1:blockLengths(ii),ii) = cycle_metadata.(metadata_fields{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end


%% Reshape into Cycles-x-Isotope Ratio-x-Block-x-Sample Aliquot
% Define Methods that Are Run at the Start of Each New Sample
sampleStartMethods = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];

idx_working = false(size(block_metadata.msDatetime,2),1);

for ii=1:length(sampleStartMethods)
    idx_working = idx_working | block_metadata.method(1,:)'==sampleStartMethods(ii);
end

% Find new aliquots by finding above methods or by finding the start of a
% new sequence (sequence row decreases from N to 1).
idx_SampleAliquots = find(idx_working | ([0; diff(block_metadata.sequenceRow(1,:)')]<0));
aliquotLengths = diff(idx_SampleAliquots); aliquotLengths(end+1)=length(block_deltas)-(idx_SampleAliquots(end)-1);

% TAKE ONLY THE ALIQUOTS WITH 4 OR 5 BLOCKS! - This causes me to lose 56
% aliquots (4360 -> 4304, ~1%), mostly with only 1 or 2 blocks in them.
idx_SampleAliquots = idx_SampleAliquots(aliquotLengths>3 & aliquotLengths<6);
aliquotLengths = aliquotLengths(aliquotLengths>3 & aliquotLengths<6);

numberOfAliquots = length(idx_SampleAliquots);
longestAliquot = max(aliquotLengths);


% Reshape
aliquot_deltas = nan(longestBlock,size(cycle_deltas,2),longestAliquot,numberOfAliquots);

for ii = 1:length(idx_SampleAliquots)
    aliquot_deltas(:,:,1:aliquotLengths(ii),ii) = block_deltas(:,:,idx_SampleAliquots(ii):idx_SampleAliquots(ii)+aliquotLengths(ii)-1);
    for jj=1:numel(metadata_fields)
        aliquot_metadata.(metadata_fields{jj})(:,1:aliquotLengths(ii),ii) = block_metadata.(metadata_fields{jj})(:,idx_SampleAliquots(ii):idx_SampleAliquots(ii)+aliquotLengths(ii)-1);
    end
end


%% Do Some Staistics on the Blocks
% There are some weird looking blocks here that plot off the top of the
% y-axis. I should take a closer look. Set a threshold for inclusion?

figure
subplot(211)
plot(block_metadata.msDatetime(1,:),blockLengths)
xlabel('Block Number'); ylabel('Number of Cycles in Block'); 
title('Block Length (by method)');
ylim([0 20])

subplot(212)
plot(squeeze(aliquot_metadata.msDatetime(1,1,:)),aliquotLengths)
xlabel('Aliquot Number'); ylabel('Number of Blocks in Aliquot'); 
title('Aliquot Lengths (by method)')
ylim([0 10])


figure
subplot(211)
semilogy(block_metadata.msDatetime(1,:),squeeze(std(block_deltas)),'.');
ylabel('Std Dev of Cycles in a Block [per mil]');
ylim([0 150]);
legend(delta_cols,'Location','N','Orientation','Horizontal');

subplot(212)
semilogy(squeeze(aliquot_metadata.msDatetime(1,1,:)),squeeze(nanstd(mean(aliquot_deltas,1),0,3)),'.');
ylabel('Std Dev of Blocks in an Aliquot [per mil]')
ylim([0 1500])


%% Calculate the PIS for each PIS experiment
block_means = nanmean(aliquot_deltas,1);

iPIS = ~isnan(block_means(:,:,:,:)); % PIS aliquots are the aliquots with no nans for any of the blocks

% First, check that there is a delta value for each measured ratio
allDeltasPisCheck = sum(squeeze(iPIS),3); % sums the number of non-nan blocks for each delta value
if length(unique(allDeltasPisCheck(:,5))) > 1 % if there is more than one unique value in the fifth column then some deltas are missing a PIS value
    warning('CAUTION: Some delta values are missing a PIS value')
end

 % Now restrict PIS blocks to just those with PIS as an identifier to weed
 % out unusual aliquots where a fifth, non-PIS block was run
iPIS = squeeze(iPIS(1,1,5,:)) & squeeze(aliquot_metadata.ID1(1,5,:))=='PIS';

calcPis = nan(size(block_means,1),size(block_means,2)-3,5,size(block_means,4));
calcPisRsq = nan(size(block_means,1),size(block_means,2)-3,5,size(block_means,4));
calcPisImbal = nan(size(block_means,4),1);

for ii=find(iPIS)' % find the indices of the PIS aliquots and loop through them, just look at the first cycle, delta value, and the fifth block for each aliquot
    for jj=4:size(block_means,2) % loop all through delta values, skip the first three columns as these are voltages and pressure imbalance

        d = squeeze(block_means(:,jj,:,ii)); % response variable = the looped delta value from the looped aliquot
        G = [ones(size(d)) squeeze(block_means(:,3,:,ii))]; % predictor variable = the pressure imbalance (col 3) from the looped variable
        m = (G'*G)\G'*d; % Calculate the PIS
        
        r_sq = corrcoef(G(:,2),d).^2; % Find the r-squared correlation coefficient for the PIS test
        [pImbal, idx] = max(abs(G(:,2))); % Find the block with the max P Imbalance
        
        calcPis(1,jj-3,1:5,ii)=m(2);
        calcPisRsq(1,jj-3,1:5,ii) = r_sq(1,2);
        calcPisImbal(ii) = pImbal * sign(G(idx,2));
    end
end


%% Filter the PIS values
% Some of the PIS values are likely to be erroneous, here they get weeded
% out.

% Remove those with a P Imbalance smaller than 100 mV
iSmallImbal = abs(calcPisImbal)<100;
calcPis(1,:,1:5,iSmallImbal) = nan;
calcPisRsq(1,:,1:5,iSmallImbal) = nan;
calcPisImbal(iSmallImbal) = nan;

% Remove those with an r-squared of less than 0.7
iBadRsq = calcPisRsq < 0.7;
calcPis(iBadRsq) = nan;
calcPisRsq(iBadRsq) = nan;
calcPisImbal(iBadRsq(1,1,1,:)) = nan;

% Manual Removal
% Remove two spurious looking d4038 values where the sign of the PIS
% changes back and forth and the magnitude jumps by two orders.
toRemove = find(aliquot_metadata.msDatetime(1,1,:)==datetime(2017,12,12,10,22,33));
calcPis(:,5,:,toRemove)=nan;

toRemove = find(aliquot_metadata.msDatetime(1,1,:)==datetime(2016,04,26,02,38,05));
calcPis(:,5,:,toRemove)=nan;


%% Plot a time-series of the PIS and related parameters
% This is useful to identify aliquots where the PIS block did no run
% correctly or where the r-squared is low, suggesting a poor determination
% of the PIS. These cases can be filtered out below.

stackedFig(3,'RelSize',[0.4 1.7 0.9],'Overlap',[-10 -10]);

stackedFigAx(1)
plot(squeeze(aliquot_metadata.msDatetime(1,1,:)),squeeze(calcPisRsq(1,:,1,:)),'s')
set(gca,'ColorOrderIndex',1)
plot(squeeze(aliquot_metadata.msDatetime(1,1,~isnan(calcPisImbal))),squeeze(calcPisRsq(1,:,1,~isnan(calcPisImbal))),'-');
ylabel('r^2');
ylim([0.7 1]);

stackedFigAx(2)
plot(squeeze(aliquot_metadata.msDatetime(1,1,:)),squeeze(calcPis(1,:,1,:)),'o')
set(gca,'ColorOrderIndex',1)
plot(squeeze(aliquot_metadata.msDatetime(1,1,~isnan(calcPisImbal))),squeeze(calcPis(1,:,1,~isnan(calcPisImbal))),'-')
legend(delta_cols(4:end),'Orientation','Horizontal','Location','North')
ylabel('PIS [per mil/per mil]');
ylim([-0.01 0.005]);

stackedFigAx(3)
plot(squeeze(aliquot_metadata.msDatetime(1,1,:)),calcPisImbal,'^')
set(gca,'ColorOrderIndex',1)
plot(squeeze(aliquot_metadata.msDatetime(1,1,~isnan(calcPisImbal))),calcPisImbal(~isnan(calcPisImbal)),'-','Color',lineCol(1));
text(squeeze(aliquot_metadata.msDatetime(1,1,:)),calcPisImbal,aliquot_metadata.ID1(1,5,:))
ylabel('Pressure Imbalance [per mil]')
ylim([-600 600])

stackedFigAx();
% Have to do some workaround stuff here because stackedFig doesn;t play
% nicely with datetime axes rulers. It won't automatically change the
% masterAx ruler to datetime
set(stackedFigAx(),'XRuler',matlab.graphics.axis.decorator.DatetimeRuler);
set(stackedFigAx(),'XLim',get(stackedFigAx(1),'XLim'));

stackedFigReset


%% Make PIS Correction
PIS = calcPis;
PIS = fillmissing(PIS,'previous',4);

figure; hold on;
plot(squeeze(aliquot_metadata.msDatetime(1,1,:)),squeeze(calcPis(1,:,1,:)),'o')
set(gca,'ColorOrderIndex',1);
plot(squeeze(aliquot_metadata.msDatetime(1,1,:)),squeeze(PIS(1,:,1,:)),'.')


block_means_pisCorr = block_means;
block_means_pisCorr(1,4:end,:,:) = block_means(1,4:end,:,:) - block_means(1,3,:,:).*PIS;


%% Calculate the Chemical Slopes

aliquot_means = nanmean(block_means_pisCorr,3);

iCS_15N = contains(squeeze(aliquot_metadata.ID1(1,1,:)),'CS') ...
          & (contains(squeeze(aliquot_metadata.ID1(1,1,:)),'15') ...
          | contains(squeeze(aliquot_metadata.ID1(1,1,:)),'N'));
iCS_18O = contains(squeeze(aliquot_metadata.ID1(1,1,:)),'CS') ...
          & (contains(squeeze(aliquot_metadata.ID1(1,1,:)),'18') ...
          | contains(squeeze(aliquot_metadata.ID1(1,1,:)),'O'));

dates = {squeeze(aliquot_metadata.msDatetime(1,1,iCS_15N)) squeeze(aliquot_metadata.msDatetime(1,1,iCS_18O))}; % just for reference

% Find the different CS experiments by finding the CS aliquots separated by more than 10 hours
diffsCS_15N = duration(nan(3223,3)); diffsCS_18O = duration(nan(3223,3));
diffsCS_15N(iCS_15N) = [diff(squeeze(aliquot_metadata.msDatetime(1,1,iCS_15N))); hours(999)+minutes(59)+seconds(59)];
diffsCS_18O(iCS_18O) = [diff(squeeze(aliquot_metadata.msDatetime(1,1,iCS_18O))); hours(999)+minutes(59)+seconds(59)];

% N.B. It's important to not use too big a number here. Using 24 hours
% fails to resolve a re-do of the 18O CS in Feb-2016 as it was run the
% morning after the previous attempt was run in the afternoon w/ diff=15 hr
endCS_15N = diffsCS_15N > 10/24; endCS_18O = diffsCS_18O > 10/24;


% dO2/N2 Effect on d15N
% Create logical indices to update within the loop as I "check off" the
% different sets of aliquots that make up the different CS experiments
iCS_15Nloop = iCS_15N;
endCS_15Nloop = endCS_15N;

figure;
CS_15N = nan(size(aliquot_means,4),1);
for ii=1:sum(endCS_15N)
    idxFinalAliquot = find(endCS_15Nloop,1);
    iAliquotsToUse = iCS_15Nloop & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    d = squeeze(aliquot_means(1,4,1,iAliquotsToUse)); % response variable = d15N
    G = [ones(size(d)) squeeze(aliquot_means(1,9,1,iAliquotsToUse))]; % predictor variable = dO2N2
    
    m = (G'*G)\G'*d; % Calculate the 15N CS
    r_sq = corrcoef(G(:,2),d).^2;
    
    CS_15N(idxFinalAliquot) = m(2);
    
    subplot(1,sum(endCS_15N),ii); hold on;
    plot(G(:,2),d,'xk')
    plot(G(:,2),G*m,'-r')
    text(50,0.015,['CS = ' num2str(m(2)*1000) ' per meg/per mil'])
    text(50,0.005,['r^2 = ' num2str(r_sq(2,1))])
    xlabel('\deltaO_2/N_2 [per mil]')
    ylabel('\delta^{15}N [per mil]')
    axis([-10 350 -0.01 0.25]);
    
    title(['\delta^{15}N CS: ' datestr(aliquot_metadata.msDatetime(1,1,idxFinalAliquot),'yyyy-mmm-dd')])    
    
    iCS_15Nloop(1:idxFinalAliquot)=false;
    endCS_15Nloop(idxFinalAliquot)=false;
end

% dO2/N2 effect on dArN2
% Reset logical indices to update within the loop as I "check off" the
% different sets of aliquots that make up the different CS experiments
iCS_15Nloop = iCS_15N;
endCS_15Nloop = endCS_15N;

figure;
CS_ArN2 = nan(size(aliquot_means,4),1);
for ii=1:sum(endCS_15N)
    idxFinalAliquot = find(endCS_15Nloop,1);
    iAliquotsToUse = iCS_15Nloop & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    d = squeeze(aliquot_means(1,10,1,iAliquotsToUse)); % response variable = dArN2
    G = [ones(size(d)) squeeze(aliquot_means(1,9,1,iAliquotsToUse))]; % predictor variable = dO2N2
    
    m = (G'*G)\G'*d; % Calculate the ArN2 CS
    r_sq = corrcoef(G(:,2),d).^2;
    
    CS_ArN2(idxFinalAliquot) = m(2);
    
    subplot(1,sum(endCS_15N),ii); hold on;
    plot(G(:,2),d,'xk')
    plot(G(:,2),G*m,'-r')
    text(50,0,['CS = ' num2str(m(2)*1000) ' per meg/per mil'])
    text(50,-0.05,['r^2 = ' num2str(r_sq(2,1))])
    xlabel('\deltaO_2/N_2 [per mil]')
    ylabel('\deltaAr/N_2 [per mil]')
    axis([-10 350 -0.1 1.2]);
    
    title(['\deltaAr/N_2 CS: ' datestr(aliquot_metadata.msDatetime(1,1,idxFinalAliquot),'yyyy-mmm-dd')])    
    
    iCS_15Nloop(1:idxFinalAliquot)=false;
    endCS_15Nloop(idxFinalAliquot)=false;
end

% dN2/O2 effect on d18O
% Create logical indices to update within the loop as I "check off" the
% different sets of aliquots that make up the different CS experiments
iCS_18Oloop = iCS_18O;
endCS_18Oloop = endCS_18O;

figure;
CS_18O = nan(size(aliquot_means,4),1);
for ii=1:sum(endCS_18O)
    idxFinalAliquot = find(endCS_18Oloop,1);
    iAliquotsToUse = iCS_18Oloop & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    d = squeeze(aliquot_means(1,5,1,iAliquotsToUse)); % response variable = d18O
    G = [ones(size(d)) ((squeeze(aliquot_means(1,9,1,iAliquotsToUse))./1000+1).^-1-1)*1000]; % predictor variable = dN2/O2 
    
    m = (G'*G)\G'*d; % Calculate the 18O CS
    r_sq = corrcoef(G(:,2),d).^2;
    
    CS_18O(idxFinalAliquot) = m(2);
    
    subplot(1,sum(endCS_18O),ii); hold on;
    plot(G(:,2),d,'xk')
    plot(G(:,2),G*m,'-r')
    text(50,-0.08,['CS = ' num2str(m(2)*1000) ' per meg/per mil'])
    text(50,-0.09,['r^2 = ' num2str(r_sq(2,1))])
    xlabel('\deltaN_2/O_2 [per mil]')
    ylabel('\delta^{18}O [per mil]')
    axis([-10 350 -0.1 0.1]);
    
    title(['\delta^{18}O CS: ' datestr(aliquot_metadata.msDatetime(1,1,idxFinalAliquot),'yyyy-mmm-dd')])    
    
    iCS_18Oloop(1:idxFinalAliquot)=false;
    endCS_18Oloop(idxFinalAliquot)=false;
end

% dN2/O2 effect on d17O
% Create logical indices to update within the loop as I "check off" the
% different sets of aliquots that make up the different CS experiments
iCS_18Oloop = iCS_18O;
endCS_18Oloop = endCS_18O;

figure;
CS_17O = nan(size(aliquot_means,4),1);
for ii=1:sum(endCS_18O)
    idxFinalAliquot = find(endCS_18Oloop,1);
    iAliquotsToUse = iCS_18Oloop & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    d = squeeze(aliquot_means(1,6,1,iAliquotsToUse)); % response variable = d17O
    G = [ones(size(d)) ((squeeze(aliquot_means(1,9,1,iAliquotsToUse))./1000+1).^-1-1)*1000]; % predictor variable = dN2/O2 
    
    m = (G'*G)\G'*d; % Calculate the 18O CS
    r_sq = corrcoef(G(:,2),d).^2;
    
    CS_17O(idxFinalAliquot) = m(2);
    
    subplot(1,sum(endCS_18O),ii); hold on;
    plot(G(:,2),d,'xk')
    plot(G(:,2),G*m,'-r')
    text(50,0,['CS = ' num2str(m(2)*1000) ' per meg/per mil'])
    text(50,-0.05,['r^2 = ' num2str(r_sq(2,1))])
    xlabel('\deltaN_2/O_2 [per mil]')
    ylabel('\delta^{17}O [per mil]')
    axis([-10 350 -0.1 1]);
    
    title(['\delta^{17}O CS: ' datestr(aliquot_metadata.msDatetime(1,1,idxFinalAliquot),'yyyy-mmm-dd')])    
    
    iCS_18Oloop(1:idxFinalAliquot)=false;
    endCS_18Oloop(idxFinalAliquot)=false;
end



%%
disp('>> Script Complete')