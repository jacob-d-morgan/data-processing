%% Spreadsheet Timeline
% Plots a timeline of the date range covered by each Excel spreadsheet
%
% Plots the time interval covered by each Excel spreadsheet, using the
% information logged in spreadsheet_metadata.xlsx. This allows me to
% identify which files to look in for useful metadata or information about
% the Mass Spec. The timeline is only as good as the information in
% spreadsheet_metadata.xlsx. As I compiled this information by hand, some
% end or, especially, start dates may be inaccurate as there are often
% extra samples included at the start of a file to record the values from
% the most recent CS or LJA experiment that may be many weeks or months
% before the start of the samples recorded in the file.

clearvars;
cluk; clc;

%% Import Start and End Dates
% Get the start and end dates of each spreadsheet, as recorded in the file
% spreadsheet_metadata.xlsx. Also generates the start and end dates of time
% intervals not covered by any of the logged spreadsheets.

[spreadsheetRanges,missingIntervals] = getSpreadsheetRanges('spreadsheet_metadata.xlsx');

%% Plot Timeline
% Plots the time interval covered by each spreadsheet and shades the
% intervals not covered by any of the logged spreadsheets.

% Plot Time Interval of each Spreadsheet
figure; hold on;
for ii = 1:height(spreadsheetRanges)
    plot(datenum([spreadsheetRanges.StartDate(ii) spreadsheetRanges.EndDate(ii)]),[1 1].*ii,'-','LineWidth',5)
end

ax = gca;
ax.YGrid = 'on';
ax.YTick = 1:height(spreadsheetRanges);
labels = replace(spreadsheetRanges.ExcelFile,'_','\_');
ax.YTickLabel = labels;

% Plot Missing Intervals
xvals = nan(4,length(missingIntervals));
xvals([1 2],1:end) = repmat(datenum(missingIntervals(:,1))',2,1);
xvals([3 4],1:end) = repmat(datenum(missingIntervals(:,2))',2,1);

yvals = nan(4,length(missingIntervals));
yvals([1 4],1:end) = ax.YLim(1);
yvals([2 3],1:end) = ax.YLim(2);

patch(xvals,yvals,[1 0.7 0.7])
datetick('x')

title('Timespan of Excel Files')