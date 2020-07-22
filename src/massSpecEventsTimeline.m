%% Mass Spec Events Timeline
% Plots a timeline of the different 'events' that happened to the Mass Spec
%
% Indicated the timing of 'events' logged in spreadsheet_metadata.xlsx,
% such as a venting, source cleaning, change in filament, refocussing etc.
% Also indicated are LJA and CS measurements. This data was gathered from
% the Excel sheets available to me from the XP machine and from
% 'cryomagic', Ross' desktop.

clearvars;
cluk; clc;

%% Import Event Dates
% Get the start and end dates of each event, as recorded in the file
% spreadsheet_metadata.xlsx

massSpecEvents = readtable('spreadsheet_metadata.xlsx','Sheet',1);
massSpecEvents.Event = categorical(massSpecEvents.Event);

%% Plot Timline
% Plots the time interval covered by each spreadsheet and shades the
% intervals not covered by any of the logged spreadsheets.

figure; hold on
iPlot = massSpecEvents.Event == "New Filament" | massSpecEvents.Event == "Refocus";
plot([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)'],repmat(get(gca,'YLim')',1,sum(iPlot)),'-r');
text(massSpecEvents.EndDate(iPlot),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color','r')

iPlot = massSpecEvents.Event == "New Std Cans" | massSpecEvents.Event == "Swap Std Cans";
plot([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)'],repmat(get(gca,'YLim')',1,sum(iPlot)),':g');
text(massSpecEvents.EndDate(iPlot),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color','g')

xlim(datetime([2014 2019],[1 1],[1 1]))