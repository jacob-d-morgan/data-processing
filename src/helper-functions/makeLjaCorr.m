function [ljaCorrDataset,LJA] = makeLjaCorr(rawDataset,ljaStats,ljaValues)
% MAKELJACORR normalizes delta values to measurements of La Jolla Air
%   For the array of values in the field 'deltas' of the structure DATASET,
%   all measured against a standard can, each value is normalized to La
%   Jolla Air using measurements of LJA against the same standard can.
%   
%   LJAVALUES should be an M-by-N matrix where M and N are the size of the
%   'deltas' field along the first two dimensions. Missing (NaN) values are
%   generally filled with the mean of all LJA measurements against the same
%   standard can, using the same focussing and filament, as logged in
%   spreadsheet_metadata.xlsx. An excpetion to this is if there is a trend
%   in the values of LJA against the standard can through time. If the
%   trend is larger than the scatter in the LJA measurements, LJAVALUES is
%   a linear interpolation/extrapolation of the measured values over the
%   correction period. If no LJA measurements were made during a given
%   correction period, the next available LJA value is used. The LJA values
%   used for each aliquot are given in the second output, LJA.
%
%   e.g. [DELTAS_CSCORR,CSS] = MAKELJACORR(DATASET,LJAVALUES) 
%
% -------------------------------------------------------------------------

LJA = ljaValues;

massSpecEvents = getMassSpecEvents;
newCorrections = massSpecEvents(...
    massSpecEvents.Event == "New Filament" | ...
    massSpecEvents.Event == "Refocus" | ...
    massSpecEvents.Event == "New Std Cans" | ...
    massSpecEvents.Event == "Swap Std Cans" | ...
    massSpecEvents.Event == "MS Change",:);

idxLjas = find(~isnat(ljaStats.ljaDatetime));
for ii = 1:length(idxLjas)
    ljaDate = ljaStats.ljaDatetime(idxLjas(ii));
    startDate = newCorrections.EndDate(find(newCorrections.EndDate < ljaDate,1,'last'));
    endDate = newCorrections.StartDate(find(newCorrections.EndDate > ljaDate,1,'first'));
    firstLjaDate = min(ljaStats.ljaAliquotSetDates{idxLjas(ii)});
    lastLjaDate = max(ljaStats.ljaAliquotSetDates{idxLjas(ii)});
    
    effect = ljaStats.ljaSlope(idxLjas(ii),:).*days(lastLjaDate - firstLjaDate);
    
    alsToUse = ljaStats.ljaAliquotSetDeltas{idxLjas(ii)};
    alsToUse(ljaStats.ljaRejections{idxLjas(ii)})=nan;
    
    ljaAtMidpoint = ljaStats.ljaIntercept(idxLjas(ii),:) + ljaStats.ljaSlope(idxLjas(ii),:).*datenum(mean([startDate endDate]));
    ljaPolyval = ljaStats.ljaIntercept(idxLjas(ii)) + ljaStats.ljaSlope(idxLjas(ii),:).*datenum(ljaStats.ljaAliquotSetDates{idxLjas(ii)});
    
    detrended = alsToUse +  ljaAtMidpoint - ljaPolyval;
    stdDetrended = nanstd(detrended);
    
    effectSize = effect./stdDetrended;
    iTemporalValues = abs(effectSize) > 1.1;
    
    aliquotDates = rawDataset.metadata.msDatetime(:,:,1,1);
    LJA(aliquotDates > startDate & aliquotDates < endDate,~iTemporalValues) = fillmissing(LJA(aliquotDates > startDate & aliquotDates < endDate,~iTemporalValues),'nearest');
    LJA(aliquotDates > startDate & aliquotDates < endDate,iTemporalValues) = ljaAtMidpoint(iTemporalValues) + ljaStats.ljaSlope(idxLjas(ii),iTemporalValues).*(datenum(aliquotDates(aliquotDates > startDate & aliquotDates < endDate) - mean([startDate endDate])));
    
end

LJA = fillmissing(LJA,'next');


ljaCorrDataset = rawDataset;
ljaCorrDataset.deltas{:,:} = ((rawDataset.deltas{:,:}/1000+1)./(LJA/1000+1)-1)*1000;
