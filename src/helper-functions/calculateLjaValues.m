function [ljaValues, ljaStats] = calculateLjaValues(dataset,iLjaToUse)
% CALCULATELJAVALUES Calculates the LJA Normalization values
%   Calculates the mean of each delta value for sets of LJA aliquots
%   analyzed in time intervals of identical Mass Spectrometer conditions.
%
%   CALCULATELJAVALUES(DATASET,iLJATOUSE) returns
%   the mean of each delta value, for each set of LJA anaylses. The
%   aliquots identified in iLJATOUSE are split into different sets based on
%   the changes to the Mass Spectrometer, documented in the Excel file
%   'spreadsheet_metadata.xlsx'. Only aliquots measured on the same
%   filament, with the same source focusing, against the same standard can
%   are averaged together.
%
%   [LJAVALUES,LJASTATS] = CALCULATELJAVALUES(...) also outputs a structure
%   of statistics and other information for each set of LJA aliquots,
%   including:
%
%   ljaDatetime - the datetime of the first block of the last aliquot in a set of LJA aliquots
%   ljaValues - the mean of all aliquots
%   ljaAliquotSetDates - the datetime of each aliquot
%   ljaAliquotSetDeltas - the mean delta values of each aliquot
%   ljaRejections - the aliquots that are rejected when calculating the mean LJA values (see below)
%   ljaSlope - the slope of a linear fit to the set aliquot means through time
%   ljaIntercept - the intercept of a linear fit to the set aliquot means through time
%   ljaRSq - the r-squared of the linear fit
%   ljaPVal - the p-value of the corelation coefficient
%   
%   canAliquotSetDates - as above but for the standard can aliquots
%   canAliquotSetDeltas - as above but for the standard can aliquots
%   canRejections - as above but for the standard can aliquots
%   canSlope - as above but for the standard can aliquots
%   canIntercept - as above but for the standard can aliquots
%   canRSq - as above but for the standard can aliquots
%   canPVal - as above but for the standard can aliquots
%
% -------------------------------------------------------------------------

%% Parse Inputs

narginchk(2,2)
nargoutchk(1,2)

iCanToUse = contains(dataset.metadata.fileNameUser(:,1,1),'Aair','IgnoreCase',true);

%% Find the Different Sets of LJA Aliquots
% Determine which LJA aliquots were measured under identical MS conditions
% (i.e. same filament, focussing, std can etc.)

% Load Information about MS Conditions
massSpecEvents = getMassSpecEvents;

% Only Include Changes to MS Conditions
newCorrections = massSpecEvents(...
    massSpecEvents.Event == "New Filament" | ...
    massSpecEvents.Event == "Refocus" | ...
    massSpecEvents.Event == "New Std Cans" | ...
    massSpecEvents.Event == "Swap Std Cans" | ...
    massSpecEvents.Event == "MS Change",:);

% Make Sure The Information on MS Conditions Covers Full Span of Data
edges = newCorrections.EndDate;
if min(dataset.metadata.msDatetime(:,1,1)) < edges(1)
    warning('Some samples fall earlier than the start of the information from the spreadsheets.')
    edges = [min(dataset.metadata.msDatetime(:,1,1)); edges];
end

if max(dataset.metadata.msDatetime(:,1,1)) > edges(end)
   warning('Some samples fall later than the end of the information from the spreadsheets.')
   edges = [edges; max(dataset.metadata.msDatetime(:,1,1))];
end

% Split LJA into Different Groups
corrPeriod = discretize(dataset.metadata.msDatetime(:,1,1),edges);

ljaGroup = nan(size(iLjaToUse,1),1); % Create a vector of nans
ljaGroup(iLjaToUse) = corrPeriod(iLjaToUse);

canGroup = nan(size(iCanToUse,1),1); % Create a vector of nans
canGroup(iCanToUse) = corrPeriod(iCanToUse);


%% Calculate the LJA Values
% Calculate the mean of all aliquots measured under identical MS
% conditions. Also fit a linear trend of each set of aliquots.

% Pre-allocate variables to be filled in the loop
calcLja = nan(size(dataset.deltas{:,:},[1 2]));

stats = struct();
stats.ljaDatetime = NaT(size(dataset.deltas{:,:},1),1);
stats.ljaValues = nan(size(dataset.deltas{:,:},[1 2]));
stats.ljaAliquotSetDates = cell(size(dataset.deltas{:,:},1),1);
stats.ljaAliquotSetDeltas = cell(size(dataset.deltas{:,:},1),1);
stats.ljaRejections = cell(size(dataset.deltas{:,:},1),1);
stats.ljaSlope = nan(size(dataset.deltas{:,:},[1 2]));
stats.ljaIntercept = nan(size(dataset.deltas{:,:},[1 2]));
stats.ljaRSq = nan(size(dataset.deltas{:,:},[1 2]));
stats.ljaPVal = nan(size(dataset.deltas{:,:},[1 2]));

stats.canAliquotSetDates = cell(size(dataset.deltas{:,:},1),1);
stats.canAliquotSetDeltas = cell(size(dataset.deltas{:,:},1),1);
stats.canRejections = cell(size(dataset.deltas{:,:},1),1);
stats.canSlope = nan(size(dataset.deltas{:,:},[1 2]));
stats.canIntercept = nan(size(dataset.deltas{:,:},[1 2]));
stats.canRSq = nan(size(dataset.deltas{:,:},[1 2]));
stats.canPVal = nan(size(dataset.deltas{:,:},[1 2]));

% Loop through different LJA Sets
ljaGroupIdxs = unique(ljaGroup(~isnan(ljaGroup)));
for ii = 1:length(ljaGroupIdxs)
    idxToMatch = ljaGroupIdxs(ii);
    idxLastAliquot = find(ljaGroup==idxToMatch & iLjaToUse,1,'last');
    
    ljaX_temp = dataset.metadata.msDatetime(ljaGroup==idxToMatch,1,1);
    ljaY_temp = mean(dataset.deltas{:,:}(ljaGroup==idxToMatch,:,:,:),[4 3]);
    
    canX_temp = dataset.metadata.msDatetime(canGroup==idxToMatch,1,1);
    canY_temp = mean(dataset.deltas{:,:}(canGroup==idxToMatch,:,:,:),[4 3]);
    
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
    
    stats.ljaDatetime(idxLastAliquot) = dataset.metadata.msDatetime(idxLastAliquot,1,1);
    stats.ljaValues(idxLastAliquot,:) = meanOfAliquots;
    stats.ljaAliquotSetDates{idxLastAliquot} = dataset.metadata.msDatetime(ljaGroup==idxToMatch,1,1);
    stats.ljaAliquotSetDeltas{idxLastAliquot} = ljaY_temp;
    stats.ljaRejections{idxLastAliquot} = iLjaRej;
    stats.ljaSlope(idxLastAliquot,:) = ljaP(1,:)./ljaMu(2,:); % Convert back from centered and scaled fit parameters
    stats.ljaIntercept(idxLastAliquot,:) = ljaP(2,:) - ljaP(1,:).*ljaMu(1,:)./ljaMu(2,:); % Convert back from centered and scaled fit parameters
    stats.ljaRSq(idxLastAliquot,:) = ljaRSq;
    stats.ljaPVal(idxLastAliquot,:) = ljaPVal;
    
    stats.canAliquotSetDates{idxLastAliquot} = dataset.metadata.msDatetime(canGroup==idxToMatch,1,1);
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
