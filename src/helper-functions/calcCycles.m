function cycles = calcCycles(importedData,isRefVar)
% CALCCYCLES converts imported ISODAT data to an array of cycle data
%   CYCLES = CALCCYCLES(IMPORTEDDATA,ISREF) calculates delta values from
%   the table of imported data generated by csvReadIsodat(). ISREFVAR is
%   the variable that indicates whether a given line is a sample (false) or
%   reference (true) measurement. It can be either (1) a logical variable
%   with one column and the same number of rows as importedData, (2) the
%   name of a variable in importedData, (3) the index of a variable in
%   importedData, or (4) a logical index with one row and the same number
%   of columns as imported data.
%
%   CYCLES is a structure containing two tables:
%       1. METADATA - a table of ISODAT metadata from the sample rows.
%       2. DELTAS - a table of delta values calculated for each cycle using
%       the raw (background corrected) voltages.
%
% -------------------------------------------------------------------------

%% Parse Inputs
% Check that the variables provided are of the correct format.

if ~istable(importedData)
    error('Invalid first input: must be a table of importedData generated by csvReadIsodat()')
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


%% Interpolate Standard Voltages onto the Sample Voltage Indices
% Interpolate the standard voltages that bracket each sample measurement
% onto the indices of the sample voltages for delta value calculation.

colsToUse = contains(importedData.Properties.VariableNames,'rIntensity');

intSA = importedData(~iRef,colsToUse);
intST = array2table(interp1(find(iRef),importedData{iRef,colsToUse},find(~iRef)));
intST.Properties.VariableNames = intSA.Properties.VariableNames;


%% Calculate Delta Values
% Calculate various delta values using the interpolated voltages.

deltas = table();
deltas.d15N = ((intSA.rIntensity29./intSA.rIntensity28)./(intST.rIntensity29./intST.rIntensity28) - 1)*1000;
deltas.d18O = ((intSA.rIntensity34./intSA.rIntensity32)./(intST.rIntensity34./intST.rIntensity32) - 1)*1000;
deltas.d17O = ((intSA.rIntensity33./intSA.rIntensity32)./(intST.rIntensity33./intST.rIntensity32) - 1)*1000;
deltas.d4036Ar = ((intSA.rIntensity40./intSA.rIntensity36)./(intST.rIntensity40./intST.rIntensity36) - 1)*1000;
deltas.d4038Ar = ((intSA.rIntensity40./intSA.rIntensity38)./(intST.rIntensity40./intST.rIntensity38) - 1)*1000;
deltas.dO2N2 = ((intSA.rIntensity32./intSA.rIntensity28)./(intST.rIntensity32./intST.rIntensity28) - 1)*1000;
deltas.dArN2 = ((intSA.rIntensity40./intSA.rIntensity28)./(intST.rIntensity40./intST.rIntensity28) - 1)*1000;

descriptions = deltas.Properties.VariableNames;
descriptions = replace(descriptions,'d','\delta');
descriptions = replace(descriptions,'15','^{15}');
descriptions = replace(descriptions,'17','^{17}');
descriptions = replace(descriptions,'18','^{18}');
descriptions = replace(descriptions,'4036','^{40}/_{36}');
descriptions = replace(descriptions,'4038','^{40}/_{38}');
descriptions = replace(descriptions,'O2N2','O_2/N_2');
descriptions = replace(descriptions,'ArN2','Ar/N_2');

deltas.Properties.VariableDescriptions = descriptions;
deltas.Properties.VariableUnits = repmat({char(8240)},1,size(deltas,2));
deltas.Properties.Description = ...
    "Delta values for each cycle of samples analyzed on the XP, calculated using the integrated intensities reported by ISODAT reprocessing.";

%% Calculate Check Raw Ratios
% Calculate the raw ratios for water and carbon dioxide using the
% interpolated voltages.

h2oCheck = table;
h2oCheck.int18SA = intSA.rIntensity18;
h2oCheck.int18ST = intST.rIntensity18;
h2oCheck.rH2O_SA = intSA.rIntensity20./intSA.rIntensity18;
h2oCheck.rH2O_ST = intST.rIntensity20./intST.rIntensity18;

co2Check = table;
co2Check.int44SA = intSA.rIntensity44;
co2Check.int44ST = intST.rIntensity44;
co2Check.rCO2O2_SA = intSA.rIntensity44./intSA.rIntensity32;
co2Check.rCO2O2_ST = intST.rIntensity44./intST.rIntensity32;

h2oCheck.Properties.VariableDescriptions = [
    "Integrated mass/charge 18 sample beam", ...
    "Integrated mass/charge 18 standard beam", ...
    "Mass charge 20/18 sample raw ratio", ...
    "Mass charge 20/18 standard raw ratio", ...
    ];
h2oCheck.Properties.VariableUnits = {'mV','mV','',''};
h2oCheck.Properties.Description = ...
    "Raw water vapour intensities and ratios of the sample and standard gas for each sample analyzed on the XP, calculated using the integrated intensities reported by ISODAT reprocessing.";

co2Check.Properties.VariableDescriptions = [
    "Integrated mass/charge 44 sample beam", ...
    "Integrated mass/charge 44 standard beam", ...
    "Mass charge 44/32 sample raw ratio", ...
    "Mass charge 44/32 standard raw ratio" ...
    ];
co2Check.Properties.VariableUnits = {'mV','mV','',''};
co2Check.Properties.Description = ...
    "Raw carbon dioxide intensities and ratios of the sample and standard gas for each sample analyzed on the XP, calculated using the integrated intensities reported by ISODAT reprocessing.";

%% Extract Cycle Metadata
% Create a table of metadata for each cycle.

metadata = table();
metadata.msDatetime = datetime(importedData.datetime(~importedData.IsRef__));
metadata.msDatenum = datenum(importedData.datenum(~importedData.IsRef__));
metadata.fileCreated = datetime(importedData.fileCreated(~importedData.IsRef__),'InputFormat',"yyyy-MM-dd_HHmm");
metadata.fileCreatedUTC = datetime(importedData.fileCreatedUTC(~importedData.IsRef__),'InputFormat',"yyyy-MM-dd_HHmm");
metadata.fileModifiedUTC = datetime(importedData.fileModifiedUTC(~importedData.IsRef__),'InputFormat',"yyyy-MM-dd_HHmm");
metadata.fileFullName = importedData.FileHeader_Filename(~importedData.IsRef__);
metadata.fileRelPath = importedData.fileRelPath(~importedData.IsRef__);
metadata.fileNameIsodat = importedData.fileNameIsodat(~importedData.IsRef__);
metadata.fileNameUser = importedData.fileNameUser(~importedData.IsRef__);
metadata.fileType = importedData.fileType(~importedData.IsRef__);

metadata.sequenceRow = importedData.SequenceRow(~importedData.IsRef__);
metadata.ASInlet = importedData.AS_SIOInlet(~importedData.IsRef__);
metadata.ID1 = importedData.Identifier1(~importedData.IsRef__);
metadata.method = importedData.Method(~importedData.IsRef__);
metadata.scriptName = importedData.ScriptName(~importedData.IsRef__);
metadata.gasConfig = importedData.GasConfiguration(~importedData.IsRef__);
metadata.gasName = importedData.GasName(~importedData.IsRef__);
metadata.int28SA = intSA.rIntensity28;
metadata.int28ST = intST.rIntensity28;
metadata.pressureImbal = (intSA.rIntensity28 - intST.rIntensity28);


metadata.Properties.VariableDescriptions = [
    "ISODAT timestamp for the block in MATLAB datetime format", ...
    "ISODAT timestamp for the block in MATLAB datenumber format", ...
    "Windows timestamp for the file creation", ...
    "Windows timestamp for the file creation in UTC", ...
    "Windows timestamp for the file modification in UTC", ...
    "ISODAT file name for the block, prepended by the Windows file creation date and relative path from the Dual-Inlet\Results folder", ...
    "Relative path to file from the folder that the PowerShell script operated on", ...
    "Original ISODAT file name", ...
    "Inferred user-entered file name after removing the characters prepended or appended by ISODAT", ...
    "File type", ...
    "Row Number of the block in the ISODAT Sequence", ...
    "SIO Autosampler Parameter Value", ...
    "Sample ID1 provided by the operator in the ISODAT Sequence file", ...
    "ISODAT Method name for the block", ...
    "ISODAT Method script name", ...
    "ISODAT Method gas configuration", ...
    "ISODAT Method gas name", ...
    "Integrated mass/charge 28 sample beam", ...
    "Integrated mass/charge 28 standard beam", ...
    "Difference in the mass/charge 28 beam, sample minus standard"    
    ];
metadata.Properties.VariableUnits = {'date & time','days since 01-Jan-0000','date & time','date & time','date & time','','','','','','','','','','','','','mV','mV','mV'};
metadata.Properties.Description = ...
    "ISODAT metadata for each cycle of samples analyzed on the XP.";

%% Assemble Output
% Combine both tables into one structure.

cycles = struct('metadata',metadata,'h2oCheck',h2oCheck,'co2Check',co2Check,'deltas',deltas);


end