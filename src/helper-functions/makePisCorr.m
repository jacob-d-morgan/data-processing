function [deltas_pisCorr,PIS] = makePisCorr(deltas_raw,aliquotDates,pImbal,pisValues)
% MAKEPISCORR corrects delta values for the pressure imbalance effect
%   For an array of delta values, DELTAS_RAW, a pressure imbalance
%   correction is performed using the pressure imbalances and the
%   calculated  pressure imbalance sensitivities given in PIMBAL and
%   CALCPIS respectively.
%   
%   CALCPIS should be an M-by-N matrix where M and N are the size of
%   DELTAS_RAW along the first two dimensions. Missing (NaN) values are
%   filled with the previous non-missing value, except for delta values
%   whose DATETIME indicates that they were measured immediately after a
%   change in filament or refocusing of the machine, as logged in
%   spreadsheet_metadata.xlsx, in which case the next non-missing value is
%   used. The PIS values used for each aliquot are given in the second
%   output, PIS.
%
%   e.g. [DELTAS_PISCORR,PIS] = MAKEPISCORR(DELTAS_RAW,ALIQUOTDATES,PIMBAL,PISVALUES) 
%
%   PIMBAL should be an M-by-P-by-Q matrix where M, P, and Q are the size
%   of DELTAS_RAW along the first, third, and fourth dimensions
%   respectively. In the case that DELTAS_RAW is a matrix, P and Q are
%   equal to 1.
%
% -------------------------------------------------------------------------


%% Make PIS Correction
% Use the calculated and filtered Pressure Imbalance Sensitivities to make
% the PIS correction.


PIS = pisValues;

massSpecEvents = readtable('spreadsheet_metadata.xlsx','Sheet',1);
newCorrections = massSpecEvents(massSpecEvents.Event == "New Filament" | massSpecEvents.Event == "Refocus",:);
for ii = 1:height(newCorrections)
idxStart = find(aliquotDates(:,1,1) >= newCorrections.StartDate(ii),1);
idxEnd = find(aliquotDates(:,1,1) >= newCorrections.StartDate(ii) & any(~isnan([nan(idxStart-1,7); PIS(idxStart:end,:)]),2),1) - 1;
PIS(idxStart:idxEnd,:) = -99;
end

PIS = fillmissing(PIS,'previous');
PIS = standardizeMissing(PIS,-99);
PIS = fillmissing(PIS,'next');

PIS = repmat(PIS,[1 1 size(deltas_raw,[3 4])]);
deltas_pisCorr = deltas_raw - permute(pImbal,[1 4 2 3]).*PIS;