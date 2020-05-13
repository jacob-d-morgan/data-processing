%% dataImport %%

%% Import Data and Fix Variable Types
% Uncomment this section when running the script for the first time in a
% session. It's easiest to then comment these lines back out and not clear
% the xp2018 etc. variables as they take a long time to load in.

clc;
set(0,'defaultFigureVisible','off');
disp('Run dataImport: Turning Figures Off');
clearvars -EXCEPT xp20*

% clearvars;
% disp({'Loading file 1: Working...'})
% xp2018 = readtable('XP-2018(excelExportIntensityJDM).csv','Delimiter',',');
% disp({'Loading file 1: Complete'}); disp({'Loading file 2: Working...'})
% xp2017 = readtable('XP-2017(excelExportIntensityJDM).csv','Delimiter',',');
% disp({'Loading file 2: Complete'}); disp({'Loading file 3: Working...'})
% xp2016 = readtable('XP-2016(excelExportIntensityJDM).csv','Delimiter',',');
% disp({'Loading file 3: Complete'}); disp({'Loading file 4: Working...'})
% % xp2015 = readtable('XP-2015(excelExportIntensityJDM).csv','Delimiter',',');
% % disp({'Loading file 4: Complete'})

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
[~,idx_blocksAll,~] = unique(cycle_metadata.filename,'stable');
blockLengthsAll = diff(idx_blocksAll); blockLengthsAll(end+1)=length(cycle_deltas)-(idx_blocksAll(end)-1);

% TAKE ONLY THE BLOCKS WITH 16 CYCLES! - This causes me to lose 156 blocks
% (17672 -> 17516, <1%), mostly with <5 cycles in them.
idx_blocks = idx_blocksAll(blockLengthsAll==16); 
blockLengths = blockLengthsAll(blockLengthsAll==16);

numberOfBlocks = length(idx_blocks);
longestBlock = max(blockLengths);

% Fill the Array of Delta Values

% Reshape Delta Values...
% ...as an Array
block_deltas = nan(longestBlock,size(cycle_deltas,2),numberOfBlocks);
tic;
for ii = length(idx_blocks):-1:1
    block_deltas(1:blockLengths(ii),:,ii) = cycle_deltas(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
end
reshape_deltaArr = toc

% ...as a Cell Array of Arrays
% Make a cell array of delta values to test vs the array
cycle_deltasCell = cell(1,10);
for ii = 1:length(cycle_deltasCell)
    cycle_deltasCell{ii} = cycle_deltas(:,ii);
end

block_deltasCellFw = cell(1,size(cycle_deltas,2));
tic;
for ii = 1:length(idx_blocks)
    for jj = 1:length(cycle_deltasCell)
        block_deltasCellFw{jj}(1:blockLengths(ii),ii) = cycle_deltasCell{jj}(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_deltaCellFw = toc

block_deltasCellBw = cell(1,size(cycle_deltas,2));
tic;
for ii = length(idx_blocks):-1:1
    for jj = 1:length(cycle_deltasCell)
        block_deltasCellBw{jj}(1:blockLengths(ii),ii) = cycle_deltasCell{jj}(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_deltaCellBw = toc

% ...as a Table
cycle_deltasTable = array2table(cycle_deltas);
cycle_deltasTable.Properties.VariableNames = delta_cols;

block_deltasTableFw = table();
tic;
for ii = 1:length(idx_blocks)
    for jj = 1:width(cycle_deltasTable)
        block_deltasTableFw.(delta_cols{jj})(1:blockLengths(ii),ii) = cycle_deltasTable.(delta_cols{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_deltaTableFw = toc

block_deltasTableBw = table();
tic;
for ii = length(idx_blocks):-1:1
    for jj = 1:width(cycle_deltasTable)
        block_deltasTableBw.(delta_cols{jj})(1:blockLengths(ii),ii) = cycle_deltasTable.(delta_cols{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_deltaTableBw = toc

% ...as a Structure
% Make a structure of delta values to test vs the array
cycle_deltasStruct = table2struct(cycle_deltasTable,'ToScalar',true);

tic;
for ii = 1:length(idx_blocks)
    for jj = 1:numel(cycle_deltasStruct)
        block_deltasStructFw.(delta_cols{jj})(1:blockLengths(ii),ii) = cycle_deltasStruct.(delta_cols{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_deltaStructFw = toc

tic;
for ii = length(idx_blocks):-1:1
    for jj = 1:numel(cycle_deltasStruct)
        block_deltasStructBw.(delta_cols{jj})(1:blockLengths(ii),ii) = cycle_deltasStruct.(delta_cols{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_deltaStructBw = toc

% Reshape Metadata...

% ...as a Cell Array of Arrays
% Make a cell array to test for speed vs the struct
cycle_metadataCell = cell(1,9);
for ii=1:length(cycle_metadataCell)
    cycle_metadataCell{ii} = cycle_metadata.(metadata_fields{ii});
end

block_metadataCellFw = cell(size(cycle_metadataCell));
tic;
for ii = 1:length(idx_blocks)
    for jj = 1:length(cycle_metadataCell)
        block_metadataCellFw{jj}(1:blockLengths(ii),ii) = cycle_metadataCell{jj}(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_metaCellFw = toc

block_metadataCellBw = cell(size(cycle_metadataCell));
tic;
for ii = length(idx_blocks):-1:1
    for jj = 1:length(cycle_metadataCell)
        block_metadataCellBw{jj}(1:blockLengths(ii),ii) = cycle_metadataCell{jj}(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_metaCellBw = toc

% ...as a Table
% Make a table to test for speed vs struct & cell
cycle_metadataTable = struct2table(cycle_metadata);

block_metadataTableFw = table();
tic;
for ii = 1:length(idx_blocks)
    for jj = 1:width(cycle_metadataTable)
        block_metadataTableFw.(metadata_fields{jj})(1:blockLengths(ii),ii) = cycle_metadataTable.(metadata_fields{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_metaTableFw = toc

block_metadataTableBw = table();
tic;
for ii = length(idx_blocks):-1:1
    for jj = 1:width(cycle_metadataTable)
        block_metadataTableBw.(metadata_fields{jj})(1:blockLengths(ii),ii) = cycle_metadataTable.(metadata_fields{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_metaTableBw = toc
        
% ...as a Structure
tic;
for ii = length(idx_blocks):-1:1
    for jj = 1:numel(metadata_fields)
        block_metadataFw.(metadata_fields{jj})(1:blockLengths(ii),ii) = cycle_metadata.(metadata_fields{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_metaStructFw = toc

tic;
for ii = length(idx_blocks):-1:1
    for jj = 1:numel(metadata_fields)
        block_metadataBw.(metadata_fields{jj})(1:blockLengths(ii),ii) = cycle_metadata.(metadata_fields{jj})(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
    end
end
reshape_metaStructBw = toc


%% Reshape into Cycles-x-Isotope Ratio-x-Block-x-Sample Aliquot
% Define Methods that Are Run at the Start of Each New Sample
sampleStartMethods = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];

idx_working = false(size(block_metadata.msDatetime,2),1);

for ii=1:length(sampleStartMethods)
    idx_working = idx_working | block_metadata.method(1,:)'==sampleStartMethods(ii);
end

% Find new aliquots by finding above methods or by finding the start of a
% new sequence (sequence row decreases from N to 1).
idx_SampleAliquotsAll = find(idx_working | ([0; diff(block_metadata.sequenceRow(1,:)')]<0));
aliquotLengthsAll = diff(idx_SampleAliquotsAll); aliquotLengthsAll(end+1)=length(block_deltas)-(idx_SampleAliquotsAll(end)-1);

% TAKE ONLY THE ALIQUOTS WITH 4 OR 5 BLOCKS! - This causes me to lose 56
% aliquots (4360 -> 4304, ~1%), mostly with only 1 or 2 blocks in them.
idx_SampleAliquots = idx_SampleAliquotsAll(aliquotLengthsAll>3 & aliquotLengthsAll<6);
aliquotLengths = aliquotLengthsAll(aliquotLengthsAll>3 & aliquotLengthsAll<6);

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

figure;
subplot(211); hold on;
plot(cycle_metadata.msDatetime(idx_blocksAll),blockLengthsAll,'-','Color',lineCol(9),'LineWidth',1)
plot(block_metadata.msDatetime(1,:),blockLengths,'-','Color',lineCol(10),'LineWidth',3)
ylabel('Number of Cycles in Block'); 
title('Block Length (by method)');
ylim([0 20])

subplot(212); hold on;
plot(squeeze(block_metadata.msDatetime(1,idx_SampleAliquotsAll)),aliquotLengthsAll,'-','Color',lineCol(9),'LineWidth',1);
plot(squeeze(aliquot_metadata.msDatetime(1,1,:)),aliquotLengths,'-','Color',lineCol(10),'LineWidth',3);
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
% Still need to:
%   1) Weed out the questionable blocks/aliquots that mess up some of the
%      chem slope experiments
%   2) Figure out how I'm going to do the CS Correction for Ar isotopes

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
CS_15N = nan(size(block_means_pisCorr,1),1,size(block_means_pisCorr,3),size(block_means_pisCorr,4));
for ii=1:sum(endCS_15N)
    idxFinalAliquot = find(endCS_15Nloop,1);
    iAliquotsToUse = iCS_15Nloop & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    d = squeeze(aliquot_means(1,4,1,iAliquotsToUse)); % response variable = d15N
    G = [ones(size(d)) squeeze(aliquot_means(1,9,1,iAliquotsToUse))]; % predictor variable = dO2N2
    
    m = (G'*G)\G'*d; % Calculate the 15N CS
    r_sq = corrcoef(G(:,2),d).^2;
    
    CS_15N(:,:,:,idxFinalAliquot) = m(2);
    
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
CS_ArN2 = nan(size(block_means_pisCorr,1),1,size(block_means_pisCorr,3),size(block_means_pisCorr,4));
for ii=1:sum(endCS_15N)
    idxFinalAliquot = find(endCS_15Nloop,1);
    iAliquotsToUse = iCS_15Nloop & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    d = squeeze(aliquot_means(1,10,1,iAliquotsToUse)); % response variable = dArN2
    G = [ones(size(d)) squeeze(aliquot_means(1,9,1,iAliquotsToUse))]; % predictor variable = dO2N2
    
    m = (G'*G)\G'*d; % Calculate the ArN2 CS
    r_sq = corrcoef(G(:,2),d).^2;
    
    CS_ArN2(:,:,:,idxFinalAliquot) = m(2);
    
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
CS_18O = nan(size(block_means_pisCorr,1),1,size(block_means_pisCorr,3),size(block_means_pisCorr,4));
for ii=1:sum(endCS_18O)
    idxFinalAliquot = find(endCS_18Oloop,1);
    iAliquotsToUse = iCS_18Oloop & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    d = squeeze(aliquot_means(1,5,1,iAliquotsToUse)); % response variable = d18O
    G = [ones(size(d)) ((squeeze(aliquot_means(1,9,1,iAliquotsToUse))./1000+1).^-1-1)*1000]; % predictor variable = dN2/O2 
    
    m = (G'*G)\G'*d; % Calculate the 18O CS
    r_sq = corrcoef(G(:,2),d).^2;
    
    CS_18O(:,:,:,idxFinalAliquot) = m(2);
    
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
CS_17O = nan(size(block_means_pisCorr,1),1,size(block_means_pisCorr,3),size(block_means_pisCorr,4));
for ii=1:sum(endCS_18O)
    idxFinalAliquot = find(endCS_18Oloop,1);
    iAliquotsToUse = iCS_18Oloop & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    d = squeeze(aliquot_means(1,6,1,iAliquotsToUse)); % response variable = d17O
    G = [ones(size(d)) ((squeeze(aliquot_means(1,9,1,iAliquotsToUse))./1000+1).^-1-1)*1000]; % predictor variable = dN2/O2 
    
    m = (G'*G)\G'*d; % Calculate the 18O CS
    r_sq = corrcoef(G(:,2),d).^2;
    
    CS_17O(:,:,:,idxFinalAliquot) = m(2);
    
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

% N.B. - These ones are a little different as there are two effects to correct for for d40/36Ar and d40/38Ar
% dN2Ar and dO2Ar effect on d4036
% Create logical indices to update within the loop as I "check off" the
% different sets of aliquots that make up the different CS experiments
iCS_18Oloop = iCS_18O;
endCS_18Oloop = endCS_18O;
iCS_15Nloop = iCS_15N;
endCS_15Nloop = endCS_15N;

figure;
CS_36Ar = nan(size(block_means_pisCorr,1),2,size(block_means_pisCorr,3),size(block_means_pisCorr,4));
for ii=1:max([sum(endCS_15N) sum(endCS_18O)])-1
    idxFinalAliquot = max([find(endCS_18Oloop,1),find(endCS_15Nloop,1)]); % Use the index of whichever CS experiment was done latest
    iAliquotsToUse = (iCS_18Oloop | iCS_15Nloop) & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    d = squeeze(aliquot_means(1,7,1,iAliquotsToUse)); % response variable = d4036Ar
    G = [ones(size(d)) ((squeeze(aliquot_means(1,10,1,iAliquotsToUse))./1000+1).^-1-1)*1000 ((squeeze(aliquot_means(1,9,1,iAliquotsToUse))/1000+1)./((squeeze(aliquot_means(1,10,1,iAliquotsToUse))./1000+1))-1)*1000]; % predictor variables = dN2/Ar AND dO2Ar (= [q_o2n2/qarn2 -1]*1000)
    
    m = (G'*G)\G'*d; % Calculate the 4036Ar CS
    r_sq = corrcoef([G(:,2:3),d]).^2;
    
    CS_36Ar(:,1,:,idxFinalAliquot) = m(2);
    CS_36Ar(:,2,:,idxFinalAliquot) = m(3);
    
    subplot(1,sum(endCS_18O),ii); hold on;
    [X,Y]=meshgrid(G(:,2),G(:,3));
    plot3(G(:,2),G(:,3),d,'xk');
    surf(X,Y,m(1)+X.*m(2)+Y*m(3));
    %text(50,0.05,['CS = ' num2str(m(2)*1000) ' per meg/per mil'])
    %text(50,0,['CS = ' num2str(m(3)*1000) ' per meg/per mil'])
    %text(50,-0.05,['r^2 = ' num2str(r_sq(2,1))])
    xlabel('\deltaN_2/Ar [per mil]')
    ylabel('\deltaO_2/Ar [per mil]')
    zlabel('\delta^{40}/_{36}Ar [per mil]')
    %axis([-10 350 -10 350 -0.1 10]);
    
    title(['\delta^{40}/{36}Ar CS: ' datestr(aliquot_metadata.msDatetime(1,1,idxFinalAliquot),'yyyy-mmm-dd')])    
    
    iCS_18Oloop(1:idxFinalAliquot)=false;
    endCS_18Oloop(1:idxFinalAliquot)=false;
    iCS_15Nloop(1:idxFinalAliquot)=false;
    endCS_15Nloop(1:idxFinalAliquot)=false;
end
%% Make the CS Corrections

CS = nan(size(block_means_pisCorr,1),7,size(block_means_pisCorr,3),size(block_means_pisCorr,4));
CS(:,:,:,:) = [CS_15N CS_18O CS_17O zeros(size(CS_15N)) zeros(size(CS_15N)) zeros(size(CS_15N)) CS_ArN2];
CS = fillmissing(CS,'previous',4);

CS_predictors = [block_means_pisCorr(:,9,:,:) ... % O2N2 CS on d15N
                 ((block_means_pisCorr(:,9,:,:)/1000+1).^-1-1)*1000 ... % N2O2 CS on d18O
                 ((block_means_pisCorr(:,9,:,:)/1000+1).^-1-1)*1000 ... % N2O2 CS on d17O
                 zeros(size(block_means_pisCorr(:,9,:,:))) ... % d4036Ar CS Below
                 zeros(size(block_means_pisCorr(:,9,:,:))) ... % d4038Ar CS Below
                 zeros(size(block_means_pisCorr(:,9,:,:))) ... % No CS Corr for dO2N2 
                 block_means_pisCorr(:,9,:,:)]; % O2N2 CS on dArN2

block_means_pisCorr_csCorr = block_means_pisCorr;
block_means_pisCorr_csCorr(:,4:end,:,:) = block_means_pisCorr_csCorr(:,4:end,:,:) - CS.*CS_predictors;

CS_36Ar = fillmissing(CS_36Ar,'previous',4);
block_means_pisCorr_csCorr(:,7,:,:) = block_means_pisCorr_csCorr(:,7,:,:) - (CS_36Ar(:,1,:,:).*((block_means_pisCorr(:,10,:,:)/1000+1).^-1-1)*1000) - (CS_36Ar(:,2,:,:).*((block_means_pisCorr(:,9,:,:)/1000+1)./(block_means_pisCorr(:,10,:,:)/1000+1)-1)*1000);

aliquot_means_pisCorr_csCorr = nanmean(block_means_pisCorr_csCorr,3);

%% Calculate the LJA Normalization Values
% Still need to:
%   1) Plot the std dev as an error bar on top of the mean values, it won't
%      work right now because errorbar won't accept datetime values...


iLJA = contains(squeeze(aliquot_metadata.ID1(1,1,:)),'LJA');

diffsLJA = duration(nan(3223,3));
diffsLJA(iLJA) = [diff(squeeze(aliquot_metadata.msDatetime(1,1,iLJA))); hours(999)+minutes(59)+seconds(59)];

endLJA = diffsLJA > 5;

iLJAloop = iLJA;
endLJAloop = endLJA;

figure; figNum = get(gcf,'Number');
ljaValues = nan(size(block_means_pisCorr_csCorr,1),size(block_means_pisCorr_csCorr,2)-3,size(block_means_pisCorr_csCorr,3),size(block_means_pisCorr,4));
ljaStd  = nan(size(block_means_pisCorr_csCorr,1),size(block_means_pisCorr_csCorr,2)-3,size(block_means_pisCorr_csCorr,3),size(block_means_pisCorr,4));
for ii = 1:sum(endLJA)
    idxFinalAliquot = find(endLJAloop,1);
    iAliquotsToUse = iLJAloop & squeeze(aliquot_metadata.msDatetime(1,1,:)) <= aliquot_metadata.msDatetime(1,1,idxFinalAliquot);
    
    ljaValues(:,:,:,idxFinalAliquot) = repmat(nanmean(aliquot_means_pisCorr_csCorr(1,4:end,1,iAliquotsToUse),4),[1 1 5]);
    ljaStd(:,:,:,idxFinalAliquot) = repmat(nanstd(aliquot_means_pisCorr_csCorr(1,4:end,1,iAliquotsToUse),0,4),[1 1 5]);
    
    for jj = 1:7
        figure(figNum)
        subplot(7,1,jj); hold on;
        plot(squeeze(aliquot_metadata.msDatetime(1,1,iAliquotsToUse)),squeeze(block_means_pisCorr_csCorr(1,jj+3,:,iAliquotsToUse)),'.','Color',lineCol(jj));
        plot(squeeze(aliquot_metadata.msDatetime(1,1,idxFinalAliquot)),squeeze(ljaValues(1,jj,1,idxFinalAliquot)),'x-k');
        %errorbar(squeeze(aliquot_metadata.msDatetime(1,1,idxFinalAliquot)),squeeze(LJA(1,jj,1,idxFinalAliquot)),squeeze(LJAstd(1,jj,1,idxFinalAliquot)),'x-k');
        ylabel([delta_cols{jj+3} '[per mil]'])
        
        figure(figNum+1)
        subplot(7,1,jj); hold on;
        plot(squeeze(aliquot_metadata.msDatetime(1,1,idxFinalAliquot)),squeeze(ljaStd(1,jj,1,idxFinalAliquot)),'x-k');
        ylabel([delta_cols{jj+3} '[per mil]'])
        
    end
    
    iLJAloop(1:idxFinalAliquot)=false;
    endLJAloop(idxFinalAliquot)=false;
    
end


%% Make the LJA Correction

LJA = nan(size(block_means_pisCorr_csCorr,1),7,size(block_means_pisCorr_csCorr,3),size(block_means_pisCorr_csCorr,4));
LJA(:,:,:,:) = ljaValues;
LJA = fillmissing(LJA,'previous',4);

block_means_pisCorr_csCorr_ljaCorr = block_means_pisCorr_csCorr;
block_means_pisCorr_csCorr_ljaCorr(1,4:end,:,:) = ((block_means_pisCorr_csCorr_ljaCorr(1,4:end,:,:)/1000+1)./(LJA/1000+1)-1)*1000;

aliquot_means_pisCorr_csCorr_ljaCorr = nanmean(block_means_pisCorr_csCorr_ljaCorr,3);

%%
set(0,'defaultFigureVisible','on');
disp('>> Script Complete: Turning Figures Back On');