function [ljaValues, ljaStats] = calculateLjaValues(aliquot_deltas,aliquot_metadata,iLjaToUse)
% CALCULATELJAVALUES calculates the mean of each set of LJA samples run in
% order to calibrate a standard can to the ultimate atmospheric standard.
%
% -------------------------------------------------------------------------

%% Parse Inputs

narginchk(3,3)
nargoutchk(1,2)

iCanToUse = contains(aliquot_metadata.filename(:,1,1),'Aair','IgnoreCase',true);

%% Find the Different Sets of LJA Aliquots
% Determine which LJA aliquots were measured under identical MS conditions
% (i.e. same filament, focussing, std can etc.)

% Load Information about MS Conditions
massSpecEvents = readtable('spreadsheet_metadata.xlsx','Sheet',1);

% Only Include Changes to MS Conditions
newCorrections = massSpecEvents(...
    massSpecEvents.Event == "New Filament" | ...
    massSpecEvents.Event == "Refocus" | ...
    massSpecEvents.Event == "New Std Cans" | ...
    massSpecEvents.Event == "Swap Std Cans" | ...
    massSpecEvents.Event == "MS Change",:);

% Make Sure The Information on MS Conditions Covers Full Span of Data
edges = newCorrections.EndDate;
if min(aliquot_metadata.msDatetime(:,1,1)) < edges(1)
    warning('Some samples fall earlier than the start of the information from the spreadsheets.')
    edges = [min(aliquot_metadata.msDatetime(:,1,1)); edges];
end

if max(aliquot_metadata.msDatetime(:,1,1)) > edges(end)
   warning('Some samples fall later than the end of the information from the spreadsheets.')
   edges = [edges; max(aliquot_metadata.msDatetime(:,1,1))];
end

% Split LJA into Different Groups
corrPeriod = discretize(aliquot_metadata.msDatetime(:,1,1),edges);

ljaGroup = nan(size(iLjaToUse,1),1); % Create a vector of nans
ljaGroup(iLjaToUse) = corrPeriod(iLjaToUse);

canGroup = nan(size(iCanToUse,1),1); % Create a vector of nans
canGroup(iCanToUse) = corrPeriod(iCanToUse);


%% Calculate the LJA Values
% Calculate the mean of all aliquots measured under identical MS
% conditions. Also calculate std and linear trend of each set of aliquots.

% Pre-allocate variables to be filled in the loop
calcLja = nan(size(aliquot_deltas,[1 2]));

stats = struct();
stats.ljaDatetime = NaT(size(aliquot_deltas,1),1);
stats.ljaValues = nan(size(aliquot_deltas,[1 2]));
stats.ljaAliquotSetDates = cell(size(aliquot_deltas,1),1);
stats.ljaAliquotSetDeltas = cell(size(aliquot_deltas,1),1);
stats.ljaRejections = cell(size(aliquot_deltas,1),1);
stats.ljaSlope = nan(size(aliquot_deltas,[1 2]));
stats.ljaIntercept = nan(size(aliquot_deltas,[1 2]));
stats.ljaRSq = nan(size(aliquot_deltas,[1 2]));
stats.ljaPVal = nan(size(aliquot_deltas,[1 2]));

stats.canAliquotSetDates = cell(size(aliquot_deltas,1),1);
stats.canAliquotSetDeltas = cell(size(aliquot_deltas,1),1);
stats.canRejections = cell(size(aliquot_deltas,1),1);
stats.canSlope = nan(size(aliquot_deltas,[1 2]));
stats.canIntercept = nan(size(aliquot_deltas,[1 2]));
stats.canRSq = nan(size(aliquot_deltas,[1 2]));
stats.canPVal = nan(size(aliquot_deltas,[1 2]));

% Loop through different LJA Sets
ljaGroupIdxs = unique(ljaGroup(~isnan(ljaGroup)));
for ii = 1:length(ljaGroupIdxs)
    idxToMatch = ljaGroupIdxs(ii);
    idxLastAliquot = find(ljaGroup==idxToMatch & iLjaToUse,1,'last');
    
    ljaX_temp = aliquot_metadata.msDatetime(ljaGroup==idxToMatch,1,1);
    ljaY_temp = mean(mean(aliquot_deltas(ljaGroup==idxToMatch,:,:,:),4),3);
    
    canX_temp = aliquot_metadata.msDatetime(canGroup==idxToMatch,1,1);
    canY_temp = mean(mean(aliquot_deltas(canGroup==idxToMatch,:,:,:),4),3);
    
    % Identify aliquots to be rejected
    iLjaRej = isoutlier(ljaY_temp,'quartiles');
    iCanRej = isoutlier(canY_temp,'quartiles');
    
    % Pre allocate variables to be filled in the loop
    meanOfAliquots = nan(1,size(ljaY_temp,2));
    
    ljaP = nan(2,size(ljaY_temp,2));
    ljaMu = nan(2,size(ljaY_temp,2));
    ljaRSq = nan(1,size(ljaY_temp,2));
    ljaPVal = nan(1,size(ljaY_temp,2));
    
    canP = nan(2,size(canY_temp,2));
    canMu = nan(2,size(canY_temp,2));
    canRSq = nan(1,size(canY_temp,2));
    canPVal = nan(1,size(canY_temp,2));
    
    % Loop through the different delta values to calculate...
    for jj = 1:size(ljaY_temp,2)
        ljaX_tempAfterRej = ljaX_temp((~iLjaRej(:,jj)));
        ljaY_tempAfterRej = ljaY_temp((~iLjaRej(:,jj)),jj);
        
        canX_tempAfterRej = canX_temp((~iCanRej(:,jj)));
        canY_tempAfterRej = canY_temp((~iCanRej(:,jj)),jj);
                
        meanOfAliquots(jj) = mean(ljaY_tempAfterRej); % ...the mean of the set of aliquots
        
        [ljaP(:,jj),~,ljaMu(:,jj)] = polyfit(datenum(ljaX_tempAfterRej),ljaY_tempAfterRej,1); % ...the trend in LJA values through time (fit after centering and scaling the nearly repeated x-values)
        [canP(:,jj),~,canMu(:,jj)] = polyfit(datenum(canX_tempAfterRej),canY_tempAfterRej,1); % ...the trend in Can values through time (fit after centering and scaling the nearly repeated x-values)
        [r_corr,p_val] = corrcoef(datenum(ljaX_tempAfterRej),ljaY_tempAfterRej); % ...the correlation and significance of the temporal trend in LJA
        ljaRSq(jj) = r_corr(1,2).^2;
        ljaPVal(jj) = p_val(1,2);
        
        [r_corr,p_val] = corrcoef(datenum(canX_tempAfterRej),canY_tempAfterRej); % ...the correlation and significance of the temporal trend in Cans
        canRSq(jj) = r_corr(1,2).^2;
        canPVal(jj) = p_val(1,2);
    end
    
    % Assign temp variables to outputs
    calcLja(idxLastAliquot,:) = meanOfAliquots;
    
    stats.ljaDatetime(idxLastAliquot) = aliquot_metadata.msDatetime(idxLastAliquot,1,1);
    stats.ljaValues(idxLastAliquot,:) = meanOfAliquots;
    stats.ljaAliquotSetDates{idxLastAliquot} = aliquot_metadata.msDatetime(ljaGroup==idxToMatch,1,1);
    stats.ljaAliquotSetDeltas{idxLastAliquot} = ljaY_temp;
    stats.ljaRejections{idxLastAliquot} = iLjaRej;
    stats.ljaSlope(idxLastAliquot,:) = ljaP(1,:)./ljaMu(2,:); % Convert back from centered and scaled fit parameters
    stats.ljaIntercept(idxLastAliquot,:) = ljaP(2,:) - ljaP(1,:).*ljaMu(1,:)./ljaMu(2,:); % Convert back from centered and scaled fit parameters
    stats.ljaRSq(idxLastAliquot,:) = ljaRSq;
    stats.ljaPVal(idxLastAliquot,:) = ljaPVal;
    
    stats.canAliquotSetDates{idxLastAliquot} = aliquot_metadata.msDatetime(canGroup==idxToMatch,1,1);
    stats.canAliquotSetDeltas{idxLastAliquot} = canY_temp;
    stats.canRejections{idxLastAliquot} = iCanRej;
    stats.canSlope(idxLastAliquot,:) = canP(1,:)./canMu(2,:); % Convert back from centered and scaled fit parameters
    stats.canIntercept(idxLastAliquot,:) = canP(2,:) - canP(1,:).*canMu(1,:)./canMu(2,:); % Convert back from centered and scaled fit parameters
    stats.canRSq(idxLastAliquot,:) = canRSq;
    stats.canPVal(idxLastAliquot,:) = canPVal;
end


%% Assign Outputs
% Assign final outputs: the calculated CS values, the recommended
% rejections, and the regression stats.

ljaValues = calcLja;
ljaStats = stats;

end %end calculateLjaValues
