%% Spreadsheet Timeline
% Plots a 'timeline' of the date range covered by each of the Excel
% spreadsheets that record the metadata for each aliquot.
% This allows me to identify which files to look in for useful metadata or
% information about the Mass Spec.

[spreadsheetRanges,missingIntervals] = getSpreadsheetRanges('spreadsheet_metadata.xlsx');

figure; hold on;
for ii = 1:height(spreadsheetRanges)
    plot(datenum([spreadsheetRanges.StartDate(ii) spreadsheetRanges.EndDate(ii)]),[1 1].*ii,'-','LineWidth',5)
end
labels = replace(spreadsheetRanges.ExcelFile,'_','\_');
ax = gca;
ax.YTick = 1:height(spreadsheetRanges);
ax.YTickLabel = labels;
grid on

xvals = nan(4,length(missingIntervals));
xvals([1 2],1:end) = repmat(datenum(missingIntervals(:,1))',2,1);
xvals([3 4],1:end) = repmat(datenum(missingIntervals(:,2))',2,1);

yvals = nan(4,length(missingIntervals));
yvals([1 4],1:end) = ax.YLim(1);
yvals([2 3],1:end) = ax.YLim(2);

patch(xvals,yvals,[1 0.7 0.7])
datetick('x')

title('Timespan of Excel Files')