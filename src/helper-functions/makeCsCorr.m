function [deltas_csCorr,CS] = makeCsCorr(deltas_raw,aliquotDates,csPredictor,csValues)
% MAKECSCORR corrects delta values for the chemical slope effect
%   For an array of delta values, DELTAS_RAW, a chemical slope
%   correction is performed using the predictor variable and the
%   calculated chemical slopes given in CSPREDICTOR and
%   CSVALUES respectively.
%   
%   CSVALUES should be an M-by-N matrix where M and N are the size of
%   DELTAS_RAW along the first two dimensions. Missing (NaN) values are
%   filled with the previous non-missing value, except for delta values
%   whose DATETIME indicates that they were measured immediately after a
%   change in filament or refocusing of the machine, as logged in
%   spreadsheet_metadata.xlsx, in which case the next non-missing value is
%   used. The CS values used for each aliquot are given in the second
%   output, CS.
%
%   e.g. [DELTAS_CSCORR,CSS] = MAKEPISCORR(DELTAS_RAW,ALIQUOTDATES,CSPREDICTOR,CSVALUES) 
%
%   CSPREDICTOR should be an array of delta values with the same number of
%   rows as DELTAS_RAW in the first dimension and with the same number of
%   columns as csValues.
%   
% -------------------------------------------------------------------------


%% Make CS Correction

CS = csValues;

massSpecEvents = getMassSpecEvents;
newCorrections = massSpecEvents(massSpecEvents.Event == "New Filament" | massSpecEvents.Event == "Refocus",:);
for ii = 1:height(newCorrections)
    if max(aliquotDates(:,1,1,1) > newCorrections.StartDate(ii))
        idxStart = find(aliquotDates(:,1,1) >= newCorrections.StartDate(ii),1); % Idx of first aliquot after a filament change/refocus
        idxEnd = find(aliquotDates(:,1,1) >= newCorrections.StartDate(ii) & any(~isnan([nan(idxStart-1,size(CS,2)); CS(idxStart:end,:)]),2),1) - 1; % Idx of last aliquot after a filament change/refocus AND before a CS
        CS(idxStart:idxEnd,:) = -99;
    else
        break
    end
end

CS = fillmissing(CS,'previous');
CS = standardizeMissing(CS,-99);
CS = fillmissing(CS,'next');

deltas_csCorr = deltas_raw - sum(CS.*csPredictor,2);

