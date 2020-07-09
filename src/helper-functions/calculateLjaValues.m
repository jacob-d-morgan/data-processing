function [calcLja, stats] = calculateLjaValues(aliquot_deltas,aliquot_metadata,iLjaToUse)
% CALCULATELJAVALUES calculates the mean of each set of LJA samples run in
% order to calibrate a standard can to the ultimate atmospheric standard.
%
% -------------------------------------------------------------------------
%calculateLjaValues(aliquot_deltas_pisCorr_csCorr, aliquot_metadata, iLja)

% Identify each set of aliquots as separated in time by >5 days
ljaAliquotSeparation = duration(nan(size(aliquot_deltas,1),3));
ljaAliquotSeparation(iLjaToUse) = [diff(aliquot_metadata.msDatetime(iLjaToUse,1,1)); hours(999)+minutes(59)+seconds(59)];

isLastAliquot = ljaAliquotSeparation > 5;
idxLastAliquot = find(isLastAliquot);

% Increment an index by one each time its a new set of aliquots
ljaCount = nan(size(aliquot_deltas,1),1);
ljaCount(isLastAliquot) = cumsum(isLastAliquot(isLastAliquot));

% Assign the same index to all the aliquots from each set of aliquots
ljaGroup = zeros(size(aliquot_deltas,1),1); % Create a vector of zeros
ljaGroup(iLjaToUse)=nan; % Assign NaN to the values I want to replace
ljaGroup(isLastAliquot) = ljaCount(isLastAliquot); % Fill in the indices of the end of each LJA aliquot
ljaGroup(iLjaToUse) = fillmissing(ljaGroup(iLjaToUse),'next'); % Replace the nans for each aliquot with the next index - must fill only the (iLJA) subset otherwise it fills some indices with interspersed zeros that are then reassigned to NaN in the next line
ljaGroup = standardizeMissing(ljaGroup,0); % Change the zeros to nans to properly represent the fact that they are missing values

% Reject Outliers from each set of aliquots
isoutlierAnonFn = @(x){(isoutlier(x,'quartiles'))};
tf=splitapply(isoutlierAnonFn,nanmean(mean(aliquot_deltas(iLjaToUse,:,:,:),4),3),ljaGroup(iLjaToUse));

iLjaAfterRej = repmat(iLjaToUse,[1 size(aliquot_deltas,2)]);
for ii=1:length(tf)
    iLjaAfterRej(ljaGroup==ii,:) = ~tf{ii};
end

% Calculate LJA Values for Normalization and Stats for each Set of Aliquots
calcLja = nan(size(aliquot_deltas,[1 2]));
stats(max(ljaGroup)) = struct();
for jj = 1:max(ljaGroup)
    for ii=1:7
        iAliquotsToUse = iLjaAfterRej(:,ii) & ljaGroup == jj;
        aliquots = mean(mean(aliquot_deltas(iAliquotsToUse,ii,:,:),4),3);
        
        calcLja(idxLastAliquot(jj),ii,:,:) = mean(aliquots);

        stats(jj).means(ii) = mean(aliquots);
        stats(jj).aliquots{ii} = aliquots;
        stats(jj).stdevs(ii) = std(aliquots);
        stats(jj).N(ii) = size(aliquots,1);
    end
end











