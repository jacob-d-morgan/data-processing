function [calcCS, stats] = calculateChemSlope(x,y,aliquot_metadata,iCStoUse,flagAr)
% CALCULATECHEMSLOPE calculates the chemical slope for a given delta value,
% Y using the predictor variable(s) X. The Chem Slope experiment data are 
% identified in X, and Y by iCSTOUSE and are split into the different
% experiments based on their separation in time. A separation of 10 or more
% hours indicates a different chem slope experiment unless FLAGAR is set to
% TRUE, in which case a separation of 24 hours is used instead.
%
% -------------------------------------------------------------------------

%% Parse Inputs
nargoutchk(2,2);
narginchk(4,5);

if nargin==4
    flagAr = false;
end

%% Find the different CS experiments
% Finding the CS aliquots separated by more than 10 hours or with a
% decrease in the replicate identifier (assumes CS 0 - CS 30 are measured
% in that order)
csAliquotSeparation = duration(nan(size(aliquot_metadata.msDatetime,1),3));
csAliquotSeparation(iCStoUse) = [diff(aliquot_metadata.msDatetime(iCStoUse,1,1)); hours(999)+minutes(59)+seconds(59)];

csAliquotReplicateDiff = nan(size(aliquot_metadata.ID1,1),1);
idx_temp = cell2mat(strfind(aliquot_metadata.ID1(iCStoUse,1,1),"CS"));
csAliquotReplicate = extractAfter(aliquot_metadata.ID1(iCStoUse,1,1),idx_temp+2);
csAliquotReplicate = double(erase(csAliquotReplicate,{'a' 'b' 'c' 'd'}));
csAliquotReplicateDiff(iCStoUse) = [diff(csAliquotReplicate); -9999];

if ~flagAr
    isLastAliquot = csAliquotSeparation > 10/24 | csAliquotReplicateDiff < 0; % Separated by more than 10 hours or if the replicate value decreases
else
    isLastAliquot = csAliquotSeparation > 24; % Unless its an Ar chem slope, in which case 24 hours is used instead
end
idxLastAliquot = find(isLastAliquot);

% Increment an index by one for each new CS experiment
csCount = nan(size(iCStoUse,1),1);
csCount(isLastAliquot) = cumsum(isLastAliquot(isLastAliquot));

% Make the grouping variable
% Assign the same index to all the aliquots from each experiment
csGroup = zeros(size(iCStoUse,1),1); % Create a vector of zeros
csGroup(iCStoUse) = nan; % Assign NaN to the values I want to replace
csGroup(isLastAliquot) = csCount(isLastAliquot); % Fill in the indices of the end of each CS experiment
csGroup(iCStoUse) = fillmissing(csGroup(iCStoUse),'next'); % Replace the nans for each aliquot with the next index - must fill only the (iCS_15N) subset otherwise it fills some indices with interspersed zeros that are then reassigned to NaN in the next line
csGroup = standardizeMissing(csGroup,0); % Change the zeros to nans to properly represent the fact that they are missing values

% Calculate the Chem Slope
calcCS = nan(size(x,[1 2]));
stats(max(csGroup)) = struct();
for ii = min(csGroup):max(csGroup)
    x_temp = mean(mean(x(csGroup==ii,:,:,:),4),3);
    y_temp = mean(mean(y(csGroup==ii,:,:,:),4),3);
    
    csFit = [x_temp ones(length(x_temp),1)]\y_temp;
    if size(x,2)==1
        [r,pVal] = corrcoef(x_temp,y_temp);
    else
        r = nan; pVal = nan;
    end
    
    calcCS(idxLastAliquot(ii),:) = csFit(1:end-1,:);
    
    stats(ii).datetime = aliquot_metadata.msDatetime(idxLastAliquot(ii));
    stats(ii).xData = x_temp;
    stats(ii).yData = y_temp;
    stats(ii).slope = csFit(1:end-1,:);
    stats(ii).intercept = csFit(end,:);
    stats(ii).corrcoef = r;
    stats(ii).pValCorrcoef = pVal;
    
end
end % end calculateChemSlope
