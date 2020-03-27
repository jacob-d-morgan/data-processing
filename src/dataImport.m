%% dataImport %%

%% Import Data and Fix Variable Types
% clearvars; clc;
% disp({'Loading file 1: Working...'})
% xp2018 = readtable('XP-2018(excelExportIntensityJDM).csv','Delimiter',',');
% disp({'Loading file 1: Complete'}); disp({'Loading file 2: Working...'})
%  xp2017 = readtable('XP-2017(excelExportIntensityJDM).csv','Delimiter',',');
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
cycle_deltas.pressure = (intSA.rIntensity28./intST.rIntensity28 - 1)*1000;

cycle_deltas.d15N = ((intSA.rIntensity29./intSA.rIntensity28)./(intST.rIntensity29./intST.rIntensity28) - 1)*1000;
cycle_deltas.d18O = ((intSA.rIntensity34./intSA.rIntensity32)./(intST.rIntensity34./intST.rIntensity32) - 1)*1000;
cycle_deltas.d17O = ((intSA.rIntensity33./intSA.rIntensity32)./(intST.rIntensity33./intST.rIntensity32) - 1)*1000;
cycle_deltas.d4036Ar = ((intSA.rIntensity40./intSA.rIntensity36)./(intST.rIntensity40./intST.rIntensity36) - 1)*1000;
cycle_deltas.d4038Ar = ((intSA.rIntensity40./intSA.rIntensity38)./(intST.rIntensity40./intST.rIntensity38) - 1)*1000;
cycle_deltas.dO2N2 = ((intSA.rIntensity32./intSA.rIntensity28)./(intST.rIntensity32./intST.rIntensity28) - 1)*1000;
cycle_deltas.dArN2 = ((intSA.rIntensity40./intSA.rIntensity28)./(intST.rIntensity40./intST.rIntensity28) - 1)*1000;

cycle_delta_cols = cycle_deltas.Properties.VariableNames;
cycle_deltas = table2array(cycle_deltas);

%% Compile Useful Metadata
% NOTE: The current approach of subsampling the colums for the ~isRef rows
% excludes the CO2 check rows. I must remember to add these back in later.
% Also, the cycles all have exactly the same metadata so this step is
% somewhat pointless.

cycle_metadata = table();
cycle_metadata.datetime = importedData.TimeCode(~importedData.IsRef__);
cycle_metadata.filename = importedData.FileHeader_Filename(~importedData.IsRef__);
cycle_metadata.sequenceRow = importedData.Row(~importedData.IsRef__);
cycle_metadata.ASInlet = importedData.AS_SIOInlet(~importedData.IsRef__);
cycle_metadata.ID1 = importedData.Identifier1(~importedData.IsRef__);
cycle_metadata.method = importedData.Method(~importedData.IsRef__);
cycle_metadata.scriptName = importedData.ScriptName(~importedData.IsRef__);
cycle_metadata.gasConfig = importedData.GasConfiguration(~importedData.IsRef__);
cycle_metadata.gasName = importedData.GasName(~importedData.IsRef__);

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
for ii = 1:length(idx_blocks)-1
    block_deltas(1:blockLengths(ii),:,ii) = cycle_deltas(idx_blocks(ii):idx_blocks(ii)+blockLengths(ii)-1,:);
end

% Keep One Line of Metadata Per Block, Rather than One Per Cycle
block_metadata = cycle_metadata(idx_blocks,:);

%% Reshape into Cycles-x-Isotope Ratio-x-Block-x-Sample Aliquot
% Define Methods that Are Run at the Start of Each New Sample
sampleStartMethods = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];

idx_working = false(size(block_metadata,1),1);
for ii=1:length(sampleStartMethods)
    idx_working = idx_working | (block_metadata.method==sampleStartMethods(ii));
end

% Find new aliquots by finding above methods or by finding the start of a
% new sequence (sequence row decreases from N to 1).
idx_SampleAliquots = find(idx_working | ([0; diff(block_metadata.sequenceRow)]<0));
aliquotLengths = diff(idx_SampleAliquots); aliquotLengths(end+1)=length(block_deltas)-(idx_SampleAliquots(end)-1);

% TAKE ONLY THE ALIQUOTS WITH 4 OR 5 BLOCKS! - This causes me to lose 56
% aliquots (4360 -> 4304, ~1%), mostly with only 1 or 2 blocks in them.
idx_SampleAliquots = idx_SampleAliquots(aliquotLengths>3 & aliquotLengths<6);
aliquotLengths = aliquotLengths(aliquotLengths>3 & aliquotLengths<6);

numberOfAliquots = length(idx_SampleAliquots);
longestAliquot = max(aliquotLengths);

figure
subplot(211)
plot(block_metadata.datetime,blockLengths)
xlabel('Block Number'); ylabel('Number of Cycles in Block'); 
title('Block Length (by method)');
ylim([0 20])

subplot(212)
plot(block_metadata.datetime(idx_SampleAliquots),aliquotLengths)
xlabel('Aliquot Number'); ylabel('Number of Blocks in Aliquot'); 
title('Aliquot Lengths (by method)')
ylim([0 10])

% Reshape
aliquot_deltas = nan(longestBlock,size(cycle_deltas,2),longestAliquot,numberOfAliquots);
for ii = 1:length(idx_SampleAliquots)-1
    aliquot_deltas(:,:,1:aliquotLengths(ii),ii) = block_deltas(:,:,idx_SampleAliquots(ii):idx_SampleAliquots(ii)+aliquotLengths(ii)-1);
end
% Could pull out metadata into a separate variable but it would have to be
% a cell array. It can't be a table as aliquot_deltas is 4D so the metadata
% would have to be 3D unless I just took one for each aliquot and excluded
% information about the different blocks within the aliquot.

%% Do Some Staistics on the Blocks
% There are some weird looking blocks here that plot off the top of the
% y-axis. I should take a closer look. Set a threshold for inclusion?

figure
subplot(211)
semilogy(block_metadata.datetime,squeeze(std(block_deltas)),'.');
ylabel('Std Dev of Cycles in a Block [per mil]');
ylim([0 150]);
legend(cycle_delta_cols,'Location','N','Orientation','Horizontal');

subplot(212)
semilogy(block_metadata.datetime(idx_SampleAliquots),squeeze(nanstd(mean(aliquot_deltas,1),0,3)),'.');
ylabel('Std Dev of Blocks in an Aliquot [per mil]')
ylim([0 1500])


%% PIS Correction
block_means = nanmean(aliquot_deltas,1);
iPisBlocks = ~isnan(block_means(:,:,5,:));

allDeltasPisCheck = sum(squeeze(iPisBlocks),2);
if length(unique(allDeltasPisCheck)) > 1
    warning('CAUTION: Some delta values are missing a PIS value')
end

pisAliquots = block_means(:,:,:,iPisBlocks(:,1,:,:)); 
% This line takes all of the aliquots that have non-NaN values for the
% fifth block, as indicated by iPisBlocks. It doesn't matter that we only
% use the first column of iPisBlocks because all the columns should be the
% same - if one of the delta values has a value for the fifth block then
% they all should have a value.

for ii=1:size(pisAliquots,4) % loop through all the PIS aliquots
    for jj=4:size(pisAliquots,2) % loop all through delta values, skip the first three columns as these are voltages and pressure imbalance
        d = squeeze(pisAliquots(:,jj,:,ii)); % response variable = the looped delta value from the looped aliquot
        G = [ones(size(d)) squeeze(pisAliquots(:,3,:,ii))]; % predictor variable = the pressure imbalance (col 3) from the looped variable
        
        m = (G'*G)\G'*d;
        pisValues(ii,jj-3)=m(2);
    end
end






%%
disp('>> Script Complete')