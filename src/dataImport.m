%% dataImport %%

%% Import Data and Fix Variable Types
xp2018 = readtable('XP-2018(excelExportIntensityJDM).csv');
xp2018.TimeCode = datetime(datenum(xp2018.TimeCode),'ConvertFrom','Datenum'); %This is probably redundant...
xp2018.Date = datestr(datenum(xp2018.Date) + datenum('31 Dec 1999')); %Correct for two-character month '0018'
xp2018.IsRef__ = logical(xp2018.IsRef__);

idxTimeCode = find(cellfun(@(varName) strcmp('TimeCode',varName),xp2018.Properties.VariableNames)==1);
xp2018 = sortrows(xp2018,idxTimeCode);

%% Check the Gas Configuration and Gas Names
% Do something here to highlight the different unique values of the Gas
% Config or the Gas Name (CAREFUL! These can be different) and ask the user
% to select which ones they want to include? For example, discard
% everything that is not 'Air+' or discard everything that is 'Ne'.
%
% Could do something similar here with folder names etc.

% xp2018 = xp2018(contains(xp2018.FileHeader_Filename,'SPICE'),:);

%% Interpolate ST Voltages and Calculate Delta Values
% NOTE: The current method of interpolating onto the ~isRef indices
% excludes the CO2 checks from the blocks I consider as they are recorded 
% in isRef as 'true'. This is probably the right way to do it for now but I
% must remember to add them back in, including their metadata, which is
% also excluded when I subsample the metadata columns for the ~isRef rows.

intSA = xp2018(~xp2018.IsRef__,1:9);
intST = array2table(interp1(find(xp2018.IsRef__), ...
                    xp2018{xp2018.IsRef__,1:9}, ...
                    find(~xp2018.IsRef__)));
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
cycle_metadata.datetime = xp2018.TimeCode(~xp2018.IsRef__);
cycle_metadata.filename = xp2018.FileHeader_Filename(~xp2018.IsRef__);
cycle_metadata.sequenceRow = xp2018.Row(~xp2018.IsRef__);
cycle_metadata.ASInlet = xp2018.AS_SIOInlet(~xp2018.IsRef__);
cycle_metadata.ID1 = xp2018.Identifier1(~xp2018.IsRef__);
cycle_metadata.method = xp2018.Method(~xp2018.IsRef__);
cycle_metadata.scriptName = xp2018.ScriptName(~xp2018.IsRef__);
cycle_metadata.gasConfig = xp2018.GasConfiguration(~xp2018.IsRef__);
cycle_metadata.gasName = xp2018.GasName(~xp2018.IsRef__);

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