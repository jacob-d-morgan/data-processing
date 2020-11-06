function masterSheet = makeMasterSheet(filesToImport)
% MAKEMASTERSHEET Imports and compiles ISODAT reprocessed .csv files into a
% working dataset organised by analysis date.
%   MASTERSHEET = MAKEMASTERSHEET(FILESTOIMPORT) generates a "Master Sheet"
%   from the data in the ISODAT-reprocessed .csv files in FILESTOIMPORT.
%
%   The master sheet is a compilation of the relevant aliquots from the
%   input files, including data from Standard Can, PIS, Chem Slope, and La
%   Jolla Air analyses, as well as any true ice core or firn air samples.
%   The aliquots are arranged by analysis date.
%
% -------------------------------------------------------------------------

%% Parse Inputs

narginchk(1,inf)


%% Import Raw Data
% Imports the 'standard' raw dataset of aliquot delta values and their
% metadata. Does not include any cycles with a gas configuration different
% to Air+ or where the 28N2 beam is either saturated or absent, or cycles
% with a non-standard (different than the mode) number of cycles or blocks.

[rawDataset,rawPisDataset] = makeRawDataset(filesToImport,'includePIS',true);


%% Make PIS Correction
% Correct the delta values in aliquot_deltas for the effect of imbalance in
% the total pressure of gas in the source.

% Calculate PIS Values
[calcPis,pisStats] = calculatePisValues(rawDataset,rawPisDataset);

% Make PIS Correction
calcPis(pisStats.rejections) = nan;
[pisCorrDataset,PIS] = makePisCorr(rawDataset,calcPis);

%% Make Chemical Slope Correction
% Correct the delta values in aliquot_deltas_pisCorr for the effect of
% different gas ratios in the source.

% Identify the CS Experiment Aliquots
exp2match = '1?[5|8]?\s?(?<CS>[N|O])\s?C?S?\s(?<rep>\d+)';
tokens = regexp(pisCorrDataset.metadata.ID1(:,1,1),exp2match,'names');

iCS = ~cellfun('isempty',tokens);
csType = strings(size(iCS));
csType(iCS) = cellfun(@(x) x.CS,tokens(iCS));
csRep = strings(size(iCS));
csRep(iCS) = cellfun(@(x) x.rep,tokens(iCS));

iCS_AddO2 = false(size(iCS));
iCS_AddO2(iCS & csType=="N") = true;

iCS_AddN2 = false(size(iCS));
iCS_AddN2(iCS & csType=="O") = true;


% == MANUALLY INCLUDE THE ONLY 2016-02-09 REP-0 IN BOTH CS EXPERIMENTS == %
iCS_AddN2(pisCorrDataset.metadata.msDatetime(:,1,1)=={'2016-02-08 13:27:59'}) = true;
% ======================================================================= %

% Calculate the Univariate (N2 & O2 Isotopes, Ar/N2 Ratio) Chem Slopes
[calcCS_15N,csStats_15N] = calculateCsValues(pisCorrDataset.deltas.dO2N2,pisCorrDataset.deltas.d15N,pisCorrDataset.metadata,iCS_AddO2);
[calcCS_ArN2,csStats_ArN2] = calculateCsValues(pisCorrDataset.deltas.dO2N2,pisCorrDataset.deltas.dArN2,pisCorrDataset.metadata,iCS_AddO2);
[calcCS_18O,csStats_18O] = calculateCsValues((1/(pisCorrDataset.deltas.dO2N2/1000+1)-1)*1000,pisCorrDataset.deltas.d18O,pisCorrDataset.metadata,iCS_AddN2);
[calcCS_17O,csStats_17O] = calculateCsValues((1/(pisCorrDataset.deltas.dO2N2/1000+1)-1)*1000,pisCorrDataset.deltas.d17O,pisCorrDataset.metadata,iCS_AddN2);

% Calculate the Bivariate (Ar Isotopes) Chem Slopes
x_temp = [((pisCorrDataset.deltas.dArN2/1000+1).^-1 - 1)*1000 ((pisCorrDataset.deltas.dO2N2/1000+1)./(pisCorrDataset.deltas.dArN2/1000+1)-1)*1000]; % predictor variables = dN2/Ar AND dO2Ar (= [q_o2n2/q_arn2 -1]*1000)
[calcCS_4036Ar,csStats_36Ar] = calculateCsValues(x_temp,pisCorrDataset.deltas.d4036Ar,pisCorrDataset.metadata,iCS_AddN2 | iCS_AddO2,true);
[calcCS_4038Ar,csStats_38Ar] = calculateCsValues(x_temp,pisCorrDataset.deltas.d4038Ar,pisCorrDataset.metadata,iCS_AddN2 | iCS_AddO2,true);

% Make CS Corrections
csValues = [{calcCS_15N} {calcCS_18O} {calcCS_17O} {calcCS_4036Ar} {calcCS_4038Ar} {calcCS_ArN2} ];
csRegressors = pisCorrDataset.deltas(:,{'d15N','d18O','d17O','d4036Ar','d4038Ar','dArN2'});
csPredictors = [
    {pisCorrDataset.deltas.dO2N2}, ...
    {((pisCorrDataset.deltas.dO2N2/1000+1).^-1-1)*1000}, ...
    {((pisCorrDataset.deltas.dO2N2/1000+1).^-1-1)*1000}, ...
    {[((pisCorrDataset.deltas.dArN2/1000+1).^-1-1)*1000 ((pisCorrDataset.deltas.dO2N2/1000+1)./(pisCorrDataset.deltas.dArN2/1000+1)-1)*1000]}, ...
    {[((pisCorrDataset.deltas.dArN2/1000+1).^-1-1)*1000 ((pisCorrDataset.deltas.dO2N2/1000+1)./(pisCorrDataset.deltas.dArN2/1000+1)-1)*1000]}, ...
    {pisCorrDataset.deltas.dO2N2}, ...
    ];

csCorr = cell(size(csValues));
CS = cell(size(csValues));
for ii = 1:length(csValues)
    [csCorr{ii},CS{ii}] = makeCsCorr(csRegressors{:,ii},pisCorrDataset.metadata.msDatetime,csPredictors{ii},csValues{ii});
end

csCorrDataset = pisCorrDataset;
csCorrDataset.deltas{:,{'d15N','d18O','d17O','d4036Ar','d4038Ar','dArN2'}}(:,:,:,:) = [csCorr{:}];
csCorrDataset.deltas.Properties.Description = join([csCorrDataset.deltas.Properties.Description ...
    "Delta values are corrected for the effect of differences in elemental ratios between the sample and standard using an empirically measured chemical slope."]);

%% Make the LJA Correction
% Correct the delta values in aliquot_deltas_pisCorr_csCorr so that they
% are measured relative to La Jolla Air.
%
% Identifies analyses of LJA made during identical MS conditions (i.e. same
% filament, std can etc.) and calculates the mean of these aliquots, after
% rejecting outliers. The mean of all the aliquots is used to normalize all
% aliquots measured under identical MS conditions, unless there is a trend
% in the LJA values. In this case, the extrapolated values are used.

% Calculate LJA Values
iLja = contains(csCorrDataset.metadata.ID1(:,1,1),'LJA');
[ljaValues,ljaStats] = calculateLjaValues(csCorrDataset,iLja);

% Make LJA Correction
[ljaCorrDataset,LJA] = makeLjaCorr(csCorrDataset,ljaStats,ljaValues);


%% Assemble Into Final Output
% Unite all of the useful data and metadata into tables and then output all
% the tables as one master table.

delta_names = string(rawDataset.deltas.Properties.VariableNames);
delta_labels = string(rawDataset.deltas.Properties.VariableDescriptions);

% Make Table of Metadata
metadata = ljaCorrDataset.metadata;

% Make Table of Raw Delta Values
rawDeltas = rawDataset.deltas;
rawDeltas.Properties.DimensionNames = {'Sample Aliquot','Isotope Ratio'};

% Make Table of Fully Corrected Delta Values
corrDeltas = ljaCorrDataset.deltas;
corrDeltas.Properties.DimensionNames = {'Sample Aliquot','Isotope Ratio'};

% Make a Table of Correction Differences
deltaCorrections.PIS = table;
deltaCorrections.PIS{:,:} = pisCorrDataset.deltas{:,:} - rawDataset.deltas{:,:};
deltaCorrections.PIS.Properties = pisCorrDataset.deltas.Properties;
deltaCorrections.PIS.Properties.Description = "Table of PIS correction amplitudes, the difference between the PIS corrected and raw delta values. Add the values in deltaCorrections.PIS to the raw delta values to get the PIS corrected delta values.";

deltaCorrections.CS = table;
deltaCorrections.CS{:,:} = csCorrDataset.deltas{:,:} - pisCorrDataset.deltas{:,:};
deltaCorrections.CS.Properties = csCorrDataset.deltas.Properties;
deltaCorrections.CS.Properties.Description = "Table of CS correction amplitudes, the difference between the CS corrected and PIS corrected delta values. Add the values in deltaCorrections.CS to the PIS corrected delta values to get the PIS and CS corrected delta values. N.B. Adding the values in deltaCorrections to the raw delta values will not give the same result as making the CS correction to the raw delta values as the amplitude of the CS correction will be different if non-PIS corrected delta values are used to calculate the chemical slope and the CS correction.";

deltaCorrections.LJA = table;
deltaCorrections.LJA{:,:} = ljaCorrDataset.deltas{:,:} - csCorrDataset.deltas{:,:};
deltaCorrections.LJA.Properties = ljaCorrDataset.deltas.Properties;
deltaCorrections.LJA.Properties.Description = "Table of LJA correction amplitudes, the difference between the LJA corrected and CS corrected delta values. Add the values in deltaCorrections.LJA to the PIS and CS corrected delta values to get the PIS, CS, and LJA corrected delta values.  N.B. Adding the values in deltaCorrections to the raw delta values will not give the same result as making the LJA correction to the raw delta values as the amplitude of the LJA correction will be different if non-PIS and/or non-CS corrected values are used to calculate the LJA correction.";

% Make Tables of PIS, CS, and LJA Values Used
pisLog = array2table(PIS);
pisLog.Properties.VariableNames = cellstr(delta_names);
pisLog.Properties.VariableUnits = repmat({[char(8240) '/' char(8240)]},1,numel(delta_names));
pisLog.Properties.VariableDescriptions = delta_labels;
pisLog.Properties.DimensionNames = {'Sample Aliquot','Isotope Ratio'};
pisLog.Properties.Description = "Table of Pressure Imbalance Sensitivities, the coefficients used to correct each sample measurement for the effect of total gas pressure imbalance between the sample and standard bellows.";

cs_names = delta_names([1:5 7]);
csLog = table;
for ii = 1:numel(cs_names)
    csLog.(cs_names(ii)) = CS{ii};
end
csLog.Properties.VariableNames = cellstr(cs_names);
csLog.Properties.VariableUnits = repmat({[char(8240) '/' char(8240)]},1,numel(cs_names));
%csLog.Properties.VariableDescriptions = delta_labels;
csLog.Properties.DimensionNames = {'Sample Aliquot','Isotope Ratio'};
csLog.Properties.Description = "Table of Chemical Slopes, the coefficients used to correct each sample measurement for the effect of partial gas pressure imbalances between the sample and standard bellows.";

ljaLog = array2table(LJA,'VariableNames',cellstr(delta_names));
ljaLog.Properties.VariableNames = cellstr(delta_names);
ljaLog.Properties.VariableUnits = repmat({char(8240)},1,numel(delta_names));
ljaLog.Properties.VariableDescriptions = delta_labels;
ljaLog.Properties.DimensionNames = {'Sample Aliquot','Isotope Ratio'};
ljaLog.Properties.Description = "Table of La Jolla Air values, the measured composition of La Jolla Air (relative to a standard can) used to normalize each sample measurement to the composition of the modern atmosphere.";

correctionCoeffs = struct('PIS',pisLog,'CS',csLog,'LJA',ljaLog);

% Make Tables of PIS, CS, and LJA Diagnostics
pisDiagnostics = struct2table(pisStats);
csDiagnostics = table(...
    struct2table(csStats_15N), ...
    struct2table(csStats_18O), ...
    struct2table(csStats_17O), ...
    struct2table(csStats_36Ar), ...
    struct2table(csStats_38Ar), ...
    struct2table(csStats_ArN2), ...
    'VariableNames',cs_names);
ljaDiagnostics = struct2table(ljaStats);

correctionDiagnostics = struct('PIS',pisDiagnostics,'CS',csDiagnostics,'LJA',ljaDiagnostics);
    
% Make the Master Data Structure
masterSheet = struct( ...
    'metadata',metadata, ...
    'rawDeltas',rawDeltas, ...
    'corrDeltas',corrDeltas, ...
    'deltaCorrections',deltaCorrections, ...
    'correctionCoeffs',correctionCoeffs, ...
    'correctionDiagnostics', correctionDiagnostics);

end
