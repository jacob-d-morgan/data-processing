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

[aliquot_deltas,metadata,aliquot_deltas_pis,aliquot_metadata_pis] = makeRawDataset(filesToImport,'includePIS',true);

aliquot_metadata = metadata.metadata;
delta_names = metadata.delta_names;
delta_labels = metadata.delta_labels;
delta_units = metadata.delta_units;

metadata_fields = string(fieldnames(aliquot_metadata))';


%% Make PIS Correction
% Correct the delta values in aliquot_deltas for the effect of imbalance in
% the total pressure of gas in the source.

% Calculate PIS Values
[calcPis,pisStats] = calculatePisValues(aliquot_deltas,aliquot_metadata,aliquot_deltas_pis,aliquot_metadata_pis);

% Make PIS Correction
calcPis(pisStats.rejections) = nan;
[aliquot_deltas_pisCorr,PIS] = makePisCorr(aliquot_deltas,aliquot_metadata.msDatetime,aliquot_metadata.pressureImbal,calcPis);

%% Make Chemical Slope Correction
% Correct the delta values in aliquot_deltas_pisCorr for the effect of
% different gas ratios in the source.

% Identify the CS Experiment Aliquots
iCS = contains(aliquot_metadata.ID1(:,1,1),'CS');
iCS_AddO2 = iCS & contains(aliquot_metadata.ID1(:,1,1),{'15','N'});
iCS_AddN2 = iCS & contains(aliquot_metadata.ID1(:,1,1),{'18','O'});

% == MANUALLY INCLUDE THE ONLY 2016-02-09 REP-0 IN BOTH CS EXPERIMENTS == %
iCS_AddN2(aliquot_metadata.msDatetime(:,1,1)=={'2016-02-08 13:27:59'}) = true;
% ======================================================================= %

% Calculate the Univariate (N2 & O2 Isotopes, Ar/N2 Ratio) Chem Slopes
[calcCS_15N,csStats_15N] = calculateCsValues(aliquot_deltas_pisCorr(:,6,:,:),aliquot_deltas_pisCorr(:,1,:,:),aliquot_metadata,iCS_AddO2);
[calcCS_ArN2,csStats_ArN2] = calculateCsValues(aliquot_deltas_pisCorr(:,6,:,:),aliquot_deltas_pisCorr(:,7,:,:),aliquot_metadata,iCS_AddO2);
[calcCS_18O,csStats_18O] = calculateCsValues((1/(aliquot_deltas_pisCorr(:,6,:,:)/1000+1)-1)*1000,aliquot_deltas_pisCorr(:,2,:,:),aliquot_metadata,iCS_AddN2);
[calcCS_17O,csStats_17O] = calculateCsValues((1/(aliquot_deltas_pisCorr(:,6,:,:)/1000+1)-1)*1000,aliquot_deltas_pisCorr(:,3,:,:),aliquot_metadata,iCS_AddN2);

% Calculate the Bivariate (Ar Isotopes) Chem Slopes
x_temp = [((aliquot_deltas_pisCorr(:,7,:,:)./1000+1).^-1-1)*1000 ((aliquot_deltas_pisCorr(:,6,:,:)/1000+1)./(aliquot_deltas_pisCorr(:,7,:,:)./1000+1)-1)*1000]; % predictor variables = dN2/Ar AND dO2Ar (= [q_o2n2/q_arn2 -1]*1000)
[calcCS_4036Ar,csStats_36Ar] = calculateCsValues(x_temp,aliquot_deltas_pisCorr(:,4,:,:),aliquot_metadata,iCS_AddN2 | iCS_AddO2,true);
[calcCS_4038Ar,csStats_38Ar] = calculateCsValues(x_temp,aliquot_deltas_pisCorr(:,5,:,:),aliquot_metadata,iCS_AddN2 | iCS_AddO2,true);

% Make CS Corrections
csValues = [{calcCS_15N} {calcCS_18O} {calcCS_17O} {calcCS_4036Ar} {calcCS_4038Ar} {calcCS_ArN2} ];
csRegressors = [
    {aliquot_deltas_pisCorr(:,delta_names=='d15N',:,:);}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='d18O',:,:)}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='d17O',:,:)}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='d4036Ar',:,:)}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='d4038Ar',:,:);}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:);}, ...
    ];
csPredictors = [
    {aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:);}, ...
    {((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1).^-1-1)*1000}, ...
    {((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1).^-1-1)*1000}, ...
    {[((aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1).^-1-1)*1000 ((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1)./(aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1)-1)*1000]}, ...
    {[((aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1).^-1-1)*1000 ((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1)./(aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1)-1)*1000]}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:);}, ...
    ];

csCorr = cell(size(csValues));
CS = cell(size(csValues));
for ii = 1:length(csValues)
    [csCorr{ii},CS{ii}] = makeCsCorr(csRegressors{ii},aliquot_metadata.msDatetime,csPredictors{ii},csValues{ii});
end

aliquot_deltas_pisCorr_csCorr = aliquot_deltas_pisCorr;
aliquot_deltas_pisCorr_csCorr(:,[1:5 7],:,:) = [csCorr{:}];


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
iLja = contains(aliquot_metadata.ID1(:,1,1),'LJA');
[ljaValues,ljaStats] = calculateLjaValues(aliquot_deltas,aliquot_metadata,iLja);

% Make LJA Correction
[aliquot_deltas_pisCorr_csCorr_ljaCorr,LJA] = makeLjaCorr(aliquot_deltas_pisCorr_csCorr,aliquot_metadata.msDatetime(:,1,1),ljaStats,ljaValues);


%% Assemble Into Final Output
% Unite all of the useful data and metadata into tables and then output all
% the tables as one master table.

% Make Table of Metadata
metadata = struct2table(aliquot_metadata);

% Make Table of Raw Delta Values
deltas_raw = table();
for ii = 1:numel(delta_names)
    deltas_raw = [deltas_raw table(aliquot_deltas(:,ii,:,:),'VariableNames',cellstr(delta_names(ii)))];
end

% Make Table fo Fully Corrected Delta Values
deltas_corr = table();
for ii = 1:numel(delta_names)
    deltas_corr = [deltas_corr table(aliquot_deltas_pisCorr_csCorr_ljaCorr(:,ii,:,:),'VariableNames',cellstr(delta_names(ii)))];
end

% Make Tables of PIS, CS, and LJA Values Used
cs_names = delta_names([1:5 7]);
csLog = table;
for ii = 1:numel(cs_names)
    csLog = [csLog table(CS{ii},'VariableNames',cellstr(cs_names(ii)))];
end

pisLog = array2table(PIS,'VariableNames',cellstr(delta_names));
ljaLog = array2table(LJA,'VariableNames',cellstr(delta_names));
correctionCoeffs = table(pisLog,csLog,ljaLog,'VariableNames',{'PIS','CS','LJA'});

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

correctionDiagnostics = table(pisDiagnostics,csDiagnostics,ljaDiagnostics,'VariableNames',{'PIS','CS','LJA'});
    
% Make the Master Table
masterSheet = table(metadata,deltas_raw,deltas_corr,correctionCoeffs,correctionDiagnostics,'VariableNames',{'metadata','deltas_raw','deltas_corr','correctionCoeffs','correctionDiagnostics'});

end
