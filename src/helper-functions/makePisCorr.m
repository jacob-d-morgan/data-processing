function [pisCorrDataset,PIS] = makePisCorr(rawDataset,pisValues)
% MAKEPISCORR corrects delta values for the pressure imbalance effect
%   For the array of delta values, in the 'deltas' field of RAWDATASET, a
%   pressure imbalance correction is performed using the pressure
%   imbalances in the 'metadata' field of RAWDATASET, and the calculated
%   pressure imbalance sensitivities in PISVALUES.
%   
%   PISVALUES should be an M-by-N matrix where M and N are the size of the
%   'deltas' field along the first two dimensions. Missing (NaN) values are
%   filled with the previous non-missing value, except for delta values
%   whose DATETIME indicates that they were measured immediately after a
%   change in filament or refocusing of the machine, as logged in
%   spreadsheet_metadata.xlsx, in which case the next non-missing value is
%   used. The PIS values used for each aliquot are given in the second
%   output, PIS.
%
%   e.g. [PISCORRDATASET,PIS] = MAKEPISCORR(RAWDATASET,PISVALUES)
%
% -------------------------------------------------------------------------


%% Make PIS Correction
% Use the calculated and filtered Pressure Imbalance Sensitivities to make
% the PIS correction.

PIS = pisValues;

% Flag Samples Run After Filament Change & Before First PIS Experiment
massSpecEvents = getMassSpecEvents;
newCorrections = massSpecEvents(massSpecEvents.Event == "New Filament" | massSpecEvents.Event == "Refocus",:);
for ii = 1:height(newCorrections)
    idxStart = find(rawDataset.metadata.msDatetime(:,1,1) >= newCorrections.StartDate(ii),1);
    idxEnd = find(rawDataset.metadata.msDatetime(:,1,1) >= newCorrections.StartDate(ii) & any(~isnan([nan(idxStart-1,size(PIS,2)); PIS(idxStart:end,:)]),2),1) - 1;
    PIS(idxStart:idxEnd,:) = -99;
end

% Fill Most Samples with Previous PIS Value
PIS = fillmissing(PIS,'previous');

% Fill Flagged Samples with First Available PIS Value
PIS = standardizeMissing(PIS,-99);
PIS = fillmissing(PIS,'next');

% Make PIS Correction
pisCorrDataset = rawDataset;
pisCorrDataset.deltas{:,:} = rawDataset.deltas{:,:} - rawDataset.metadata.pressureImbal.*PIS;
