%% dataImport %%

%% Import Data and Fix Variable Types
%xp2018 = readtable('XP-2018(excelExportIntensityJDM).csv');
xp2018.TimeCode = datetime(datenum(xp2018.TimeCode),'ConvertFrom','Datenum'); %This is probably redundant...
xp2018.Date = datenum(xp2018.Date) + datenum('31 Dec 1999'); %Correct for two-character month '0018'
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
metadata = table();
metadata.datetime = xp2018.TimeCode(~xp2018.IsRef__);
metadata.filename = xp2018.FileHeader_Filename(~xp2018.IsRef__);
metadata.sequenceRow = xp2018.Row(~xp2018.IsRef__);
metadata.ASInlet = xp2018.AS_SIOInlet(~xp2018.IsRef__);
metadata.ID1 = xp2018.Identifier1(~xp2018.IsRef__);
metadata.method = xp2018.Method(~xp2018.IsRef__);
metadata.scriptName = xp2018.ScriptName(~xp2018.IsRef__);
metadata.gasConfig = xp2018.GasConfiguration(~xp2018.IsRef__);
metadata.gasName = xp2018.GasName(~xp2018.IsRef__);

%% Reshape into Cycles-x-Isotope Ratio-x-Block-x-Sample
% New Sample Every Time:
%    (a) After a PIS, i.e. the row after ID1=PIS
%    (b) ID1 changes, but not if it changes to PIS: ID1(ii)~=ID1(ii-1) & ID1(ii)~=PIS
%    (c) Every time the method is Automation_SA or Automation_SA_Delay


for ii=1:height(metadata)
    findStart = strfind(metadata.filename{ii},'SPICE-2018')+length('SPICE-2018_');
    metadata.sampleNameChar{ii} = metadata.filename{ii}(findStart:end-9);
end
metadata.sampleName = string(metadata.sampleNameChar);
[~,idx_aliquots,~] = unique(metadata.sampleNameChar,'stable');

% metadata.blockDuration = NaT(height(metadata),1)-NaT(1);
% metadata.blockDuration(2:end) = diff(metadata.datetime);
% idx_aliquots = (metadata.blockDuration>minutes(45));

% idx_aliquots = false(height(metadata),1);
% idx_aliquots(1) = true;
% for ii = 2:height(metadata)
%     if strcmp(metadata.ID1(ii-1),'PIS') && ~strcmp(metadata.ID1(ii),'PIS'); idx_aliquots(ii)=true; end
%     if ~strcmp(metadata.ID1(ii-1),metadata.ID1(ii)) && ~strcmp(metadata.ID1(ii),'PIS'); idx_aliquots(ii)=true; end
%     if contains(metadata.method(ii),'Automation_SA') && ~contains(metadata.method(ii-1),'Automation_SA'); idx_aliquots(ii)=true; end
% end
%idx_aliquots = find(idx_aliquots);

numberOfAliquots = length(idx_aliquots);
aliquotLengths = diff(idx_aliquots);
longestAliquot = max(aliquotLengths);

[~,idx_blocks,~] = unique(metadata.fileName,'stable');
numberOfBlocks = length(idx_blocks);
blockLengths = diff(idx_blocks);
longestBlock = max(blockLengths);

block_deltas = nan(longestBlock,size(cycle_deltas,2),numberOfBlocks);

for ii = 1:length(idx_blocks)-1
    block_deltas(1:blockLengths(ii),:,ii) = cycle_deltas(idx_blocks(ii):idx_blocks(ii+1)-1,:);
end

metadata = metadata(idx_blocks,:);




