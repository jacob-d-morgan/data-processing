function [spreadsheetRanges, missingIntervals] = getSpreadsheetRanges(fileName)
% GETSPREADSHEETRANGES - Import ranges of spreadsheets and missing time intervals
%   Reads in the Excel file that contains the ranges of the different
%   spreadsheets and spits out the table in a MATLAB table. 
%   Also outputs the start and end dates of time intervals (days) that do
%   not fall between the start and end date of any of the spreadsheets
%   logged in the Excel file.
%
% -------------------------------------------------------------------------

% Read in and Arrange the Metadata from the Excel File
spreadsheetRanges = readtable(fileName,'Sheet',2);
spreadsheetRanges(ismissing(spreadsheetRanges.StartDate) & ismissing(spreadsheetRanges.EndDate),:) = [];
spreadsheetRanges = sortrows(spreadsheetRanges,[4 5],"ascend");

% Determine which days are missing from the spreadsheets
totalRange = (min(spreadsheetRanges.StartDate):max(spreadsheetRanges.EndDate))';
isTimeCovered = false(size(totalRange));
for ii = 1:height(spreadsheetRanges)
    isTimeCovered(totalRange>=spreadsheetRanges.StartDate(ii) & totalRange<=spreadsheetRanges.EndDate(ii)) = true;
end

% Find the first and last day in each sequence of consecutive missing days
threshold = 1; % Minimum number of consecutive days missing from spreadsheets needed to define a new missing interval
d = diff([1; isTimeCovered; 1]);
startIndex = find(d < 0);
endIndex = find(d > 0)-1;
durations = endIndex - startIndex + 1;

startIndex = startIndex(durations >= threshold);
endIndex = endIndex(durations >= threshold);

missingIntervals = [totalRange(startIndex) totalRange(endIndex)];