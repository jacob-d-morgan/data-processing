function [calcLja, stats] = calculateLjaValues(aliquot_deltas,aliquot_metadata,iLjaToUse)
% CALCULATELJAVALUES calculates the mean of each set of LJA samples run in
% order to calibrate a standard can to the ultimate atmospheric standard.
%
% -------------------------------------------------------------------------

%% Parse Inputs

narginchk(3,3)
nargoutchk(1,2)

massSpecEvents = readtable('spreadsheet_metadata.xlsx','Sheet',1);
newCorrections = massSpecEvents(...
    massSpecEvents.Event == "New Filament" | ...
    massSpecEvents.Event == "Refocus" | ...
    massSpecEvents.Event == "New Std Cans" | ...
    massSpecEvents.Event == "Swap Std Cans",:);

% Identify each set of aliquots measured against the same std can
edges = newCorrections.EndDate;
if min(aliquot_metadata.msDatetime(:,1,1)) < edges(1)
    warning('Some samples fall earlier than the start of the information from the spreadsheets.')
    edges = [min(aliquot_metadata.msDatetime(:,1,1)); edges];
end

if max(aliquot_metadata.msDatetime(:,1,1)) > edges(end)
   warning('Some samples fall later than the end of the information from the spreadsheets.')
   edges = [edges; max(aliquot_metadata.msDatetime(:,1,1))];
end

corrPeriodCounts = histcounts(aliquot_metadata.msDatetime(iLjaToUse,1,1),edges)';
edgesToUse = [edges(corrPeriodCounts>0); edges(end)];
corrPeriod = discretize(aliquot_metadata.msDatetime(iLjaToUse,1,1),edgesToUse);

for ii = 1:length(edgesToUse)-1
    idxLastAliquot(ii) = find(aliquot_metadata.msDatetime(:,1,1) < edgesToUse(ii+1) & iLjaToUse,1,'last');
end
ljaGroup = nan(size(iLjaToUse,1),1); % Create a vector of zeros
ljaGroup(iLjaToUse) = corrPeriod;

% Pre-allocate variables to be filled in the loop
calcLja = nan(size(aliquot_deltas,[1 2]));

stats = struct();
stats.ljaDatetime = NaT(size(aliquot_deltas,1),1);
stats.lja = nan(size(aliquot_deltas,[1 2]));
stats.aliquotDates = cell(size(aliquot_deltas,1),1);
stats.aliquotMeans = cell(size(aliquot_deltas,1),1);
stats.rejections = cell(size(aliquot_deltas,1),1);
stats.slope = nan(size(aliquot_deltas,[1 2]));
stats.intercept = nan(size(aliquot_deltas,[1 2]));
stats.rSq = nan(size(aliquot_deltas,[1 2]));
stats.pVal = nan(size(aliquot_deltas,[1 2]));

% Loop through different LJA Sets
for ii = min(ljaGroup):max(ljaGroup)
    x_temp = aliquot_metadata.msDatetime(ljaGroup==ii,1,1);
    y_temp = mean(mean(aliquot_deltas(ljaGroup==ii,:,:,:),4),3);
    
    iRej = isoutlier(y_temp,'quartiles');
    
    p_time = nan(2,size(y_temp,2));
    rSq = nan(1,size(y_temp,2));
    pVal = nan(1,size(y_temp,2));
    for jj = 1:size(y_temp,2)
        x_tempAfterRej = x_temp((~iRej(:,jj)));
        y_tempAfterRej = y_temp((~iRej(:,jj)),jj);
                
        meanOfAliquots(jj) = mean(y_tempAfterRej);
        p_time(:,jj) = polyfit(datenum(x_tempAfterRej),y_tempAfterRej,1)';
        [r_corr,p_val] = corrcoef(datenum(x_tempAfterRej),y_tempAfterRej);
        rSq(jj) = r_corr(1,2).^2;
        pVal(jj) = p_val(1,2);
    end
    
    calcLja(idxLastAliquot(ii),:) = meanOfAliquots;
    
    stats.ljaDatetime(idxLastAliquot(ii)) = aliquot_metadata.msDatetime(idxLastAliquot(ii),1,1);
    stats.lja(idxLastAliquot(ii),:) = meanOfAliquots;
    stats.aliquotDates{idxLastAliquot(ii)} = aliquot_metadata.msDatetime(ljaGroup==ii,1,1);
    stats.aliquotMeans{idxLastAliquot(ii)} = y_temp;
    stats.rejections{idxLastAliquot(ii)} = iRej;
    stats.slope(idxLastAliquot(ii),:) = p_time(1,:);
    stats.intercept(idxLastAliquot(ii),:) = p_time(2,:);
    stats.rSq(idxLastAliquot(ii),:) = rSq;
    stats.pVal(idxLastAliquot(ii),:) = pVal;
        
end












