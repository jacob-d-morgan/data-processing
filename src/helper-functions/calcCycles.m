function [cycle_deltas,cycle_metadata] = calcCycles(importedData,isRefVar)
% CALCCYCLES converts imported ISODAT data to cycle deltas and metadata
%   [CYCLEDELTAS,CYCLEMETADATA] = CALCCYCLES(IMPORTEDDATA,ISREF) calculates
%   delta values from the table of imported data. 
%   ISREF is the variable that indicates whether a given line is a sample
%   (false) or reference (true) measurement. It can be either (1) a logical
%   variable with one column and the same number of rows as importedData,
%   (2) the name of a variable in importedData, (3) the index of a variable
%   in importedData, or (4) a logical index with one row and the same
%   number of columns as imported data.
%   CYCLEMETADATA is the metadata from the sample rows.


% -------------------------------------------------------------------------

% Check inputs
if ~istable(importedData)
    error('Invalid first input: must be a table of importedData')
end

if islogical(isRefVar)
    if size(importedData,1)==size(isRefVar,1)
        iRef = isRefVar;
    elseif size(importedData,2)==size(isRefVar,2)
        iRef = importedData(:,isRefVar);
    else
        error('Invalid logical second input: must contain the same number of columns or rows as the first input.')
    end
elseif ischar(isRefVar) || isStringScalar(isRefVar)
    if ismember(isRefVar,importedData.Properties.VariableNames)
        iRef = importedData.(isRefVar);
    else
        error('Invalid text second input: must be the (case-sepcific) name of a variable in the first input.')
    end
elseif isscalar(isRefVar)
    if ismember(isRefVar,1:size(importedData,2))
        iRef = importedData(:,isRefVar);
    else
        error('Invalid numeric second input: must be the index of a variable in the first input.')
    end
else
    error('Invalid second input: see help for valid forms')
end


% Interpolate the standard measurements onto the sample measurements
colsToUse = contains(importedData.Properties.VariableNames,'rIntensity');

intSA = importedData(~iRef,colsToUse);
intST = array2table(interp1(find(iRef),importedData{iRef,colsToUse},find(~iRef)));
intST.Properties.VariableNames = intSA.Properties.VariableNames;


% Calculate Delta Values
cycle_deltas = table();
cycle_deltas.d15N = ((intSA.rIntensity29./intSA.rIntensity28)./(intST.rIntensity29./intST.rIntensity28) - 1)*1000;
cycle_deltas.d18O = ((intSA.rIntensity34./intSA.rIntensity32)./(intST.rIntensity34./intST.rIntensity32) - 1)*1000;
cycle_deltas.d17O = ((intSA.rIntensity33./intSA.rIntensity32)./(intST.rIntensity33./intST.rIntensity32) - 1)*1000;
cycle_deltas.d4036Ar = ((intSA.rIntensity40./intSA.rIntensity36)./(intST.rIntensity40./intST.rIntensity36) - 1)*1000;
cycle_deltas.d4038Ar = ((intSA.rIntensity40./intSA.rIntensity38)./(intST.rIntensity40./intST.rIntensity38) - 1)*1000;
cycle_deltas.dO2N2 = ((intSA.rIntensity32./intSA.rIntensity28)./(intST.rIntensity32./intST.rIntensity28) - 1)*1000;
cycle_deltas.dArN2 = ((intSA.rIntensity40./intSA.rIntensity28)./(intST.rIntensity40./intST.rIntensity28) - 1)*1000;

descriptions = cycle_deltas.Properties.VariableNames;
descriptions = replace(descriptions,'d','\delta');
descriptions = replace(descriptions,'15','^{15}');
descriptions = replace(descriptions,'17','^{17}');
descriptions = replace(descriptions,'18','^{18}');
descriptions = replace(descriptions,'4036','^{40}/_{36}');
descriptions = replace(descriptions,'4038','^{40}/_{38}');
descriptions = replace(descriptions,'O2N2','O_2/N_2');
descriptions = replace(descriptions,'ArN2','Ar/N_2');
descriptions = append(descriptions,' [',char(8240),']');

cycle_deltas.Properties.VariableDescriptions = descriptions;
cycle_deltas.Properties.VariableUnits = repmat({'per mil'},1,size(cycle_deltas,2));


% Extract Cycle Metadata
cycle_metadata.msDatetime = datetime(importedData.datetime(~importedData.IsRef__));
cycle_metadata.msDatenum = datenum(importedData.datenum(~importedData.IsRef__));
cycle_metadata.filename = importedData.FileHeader_Filename(~importedData.IsRef__);
cycle_metadata.sequenceRow = importedData.SequenceRow(~importedData.IsRef__);
cycle_metadata.ASInlet = importedData.AS_SIOInlet(~importedData.IsRef__);
cycle_metadata.ID1 = importedData.Identifier1(~importedData.IsRef__);
cycle_metadata.method = importedData.Method(~importedData.IsRef__);
cycle_metadata.scriptName = importedData.ScriptName(~importedData.IsRef__);
cycle_metadata.gasConfig = importedData.GasConfiguration(~importedData.IsRef__);
cycle_metadata.gasName = importedData.GasName(~importedData.IsRef__);
cycle_metadata.int28SA = intSA.rIntensity28;
cycle_metadata.int28ST = intST.rIntensity28;
%cycle_metadata.pressureImbal = (intSA.rIntensity28./intST.rIntensity28 - 1)*1000; % N.B. - Ross typically calculates this just as a raw SA - ST rather than a delta value. Does this make a difference?
cycle_metadata.pressureImbal = (intSA.rIntensity28 - intST.rIntensity28);




end