%% dataImport %%

%% Import Data and Fix Variable Types
xp2018 = readtable('XP-2018(excelExportIntensityJDM).csv','Delimiter',',');
disp({'Loading file 1...'})
xp2017 = readtable('XP-2017(excelExportIntensityJDM).csv','Delimiter',',');
disp({'Loading file 2...'})
xp2016 = readtable('XP-2016(excelExportIntensityJDM).csv','Delimiter',',');
disp({'Loading file 3...'})
xp2015 = readtable('XP-2015(excelExportIntensityJDM).csv','Delimiter',',');
disp({'Loading file 4...'})
%%
xp2018.MeasurmentErrors = num2cell(xp2018.MeasurmentErrors);
xp2016.MeasurmentErrors = num2cell(xp2016.MeasurmentErrors);
importedData = [xp2018; xp2017; xp2016];

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

% importedData = importedData(contains(importedData.FileHeader_Filename,'SPICE'),:);

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

cycle_deltas = table2array(cycle_deltas);

%% Compile Useful Metadata
% NOTE: The current approach of subsampling the colums for the ~isRef rows
% excludes the CO2 check rows. I must remember to add these back in later.
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
numberOfBlocks = length(idx_blocks);
blockLengths = diff(idx_blocks);
longestBlock = max(blockLengths);

% Fill the Array of Delta Values
block_deltas = nan(longestBlock,size(cycle_deltas,2),numberOfBlocks);
for ii = 1:length(idx_blocks)-1
    block_deltas(1:blockLengths(ii),:,ii) = cycle_deltas(idx_blocks(ii):idx_blocks(ii+1)-1,:);
end

% Keep One Line of Metadata Per Block, Rather than One Per Cycle
block_metadata = cycle_metadata(idx_blocks,:);

% Define Methods that Are Run at the Start of Each New Sample
sampleStartMethods = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];

idx_working = false(size(block_metadata,1),1);
for ii=1:length(sampleStartMethods)
    idx_working = idx_working | (block_metadata.method==sampleStartMethods(ii));
end
idx_SampleAliquots = find(idx_working | ([0; diff(block_metadata.sequenceRow)]<0));

numberOfAliquots = length(idx_SampleAliquots);
aliquotLengths = diff(idx_SampleAliquots);
longestAliquot = max(aliquotLengths);

figure
subplot(211)
plot(blockLengths)
xlabel('Block Number'); ylabel('Number of Cycles in Block'); title('Block Length (by filename)');
xlim([-inf inf]); ylim([0 20])
subplot(212)
plot(aliquotLengths)
xlabel('Aliquot Number'); ylabel('Number of Blocks in Aliquot'); title('Aliquot Lengths (by method)')
xlim([-inf inf]);ylim([0 10])

% Reshape
aliquot_deltas = nan(longestBlock,size(cycle_deltas,2),longestAliquot,numberOfAliquots);
for ii = 1:length(idx_SampleAliquots)-1
    aliquot_deltas(:,:,1:aliquotLengths(ii),ii) = block_deltas(:,:,idx_SampleAliquots(ii):idx_SampleAliquots(ii+1)-1);
end



disp('>> Script Complete')