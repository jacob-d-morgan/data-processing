%% Sample Type Timeline
% Plots a timeline of each sample type in the imported data. 
% 
% This allows me to identify the time intervals in which each ice core was
% analyzed and to identify the different batches of CS, LJA, and see which
% ice cores they are best applied to. I can compare the timing of ice core
% analyses and CS/LJA experiments to the dates of significant changes to
% the MS, such as a refocussing, or significant adjustments to the method,
% such as the implementation of the O2-consumption correction.

% clearvars;
cluk; clc;

%% Import Data
% Loads in all the data from the .csv files and calculates the delta values
% and metadata. Also generates the 'standard' raw dataset in order to
% determine which cycles are omitted.
% Rather than reading in the files each time, it's easier to read them in
% once using the csvRead(), calcCycles(), and makeRawDataset() lines and
% then use the clearvars -EXCEPT line to preserve the variables read in by
% those lines for future runs. Toggle the comments on the indicated lines
% to adjust whether or not the files are re-read.

clearvars -EXCEPT table_deltas cycles.metadata aliquot_deltas metadata ;

filesToImport = [
    "XP-2020(exportJDM).csv"; ...
    "XP-2019(exportJDM).csv"; ...
    "XP-2018(exportJDM).csv"; ...
    "XP-2017(exportJDM).csv"; ...
    "XP-2016(exportJDM).csv"; ...
    "XP-2015(exportJDM).csv"; ...
    "XP-2014(exportJDM).csv"; ...
    "XP-2013(exportJDM).csv"; ...
    "XP-2012(exportJDM).csv"; ...
    "XP-2011(exportJDM).csv"; ...
    "XP-2010(exportJDM).csv"; ...
    "XP-2009(exportJDM).csv"; ...
    "XP-2008(exportJDM).csv"; ...
    "XP-2007(exportJDM).csv"; ...
    "XP-2006(exportJDM).csv"; ...
    ];

% Load in Data from .csv Files
importedData = csvReadIsodat(filesToImport);
cycles = calcCycles(importedData,'IsRef__');

% Generate 'Standard' Raw Dataset
rawDataset = makeRawDataset(filesToImport);

%% Define and Assign Categories
% Assign a sample type category to each sample and also determine whether
% or not it is included in the 'standard' raw dataset. Samples from
% cycle_deltas that are not included are those with a gas configuration
% different to Air+ or where the 28N2 beam is saturated or absent and those
% with a non-standard number of cycles or blocks.

% Define Categories
sampleInclusionCats = {'Included','Not Included','PIS'};
sampleTypeCats = {'Automation Test','Misc. Test','Jeff Misc.','Keeling Tanks Misc.','Sarah Misc.','Alan Misc.','Kenji Misc.','Vas Misc.','Wendy Misc.','Jinho Misc.','Jeff Neon','Std Can','Chem Slope','LJA','SPC','Dome Fuji','WDC','GISP','RICE','NEEM','BYRD','SDMA','Vostok','Bruce Plateau','Taylor Glacier','Law Dome','Minna Bluff','Summit'};

% Make Empty Categorical Arrays to Fill
sampleInclusion = setcats(repmat(categorical(missing),height(cycles.metadata),1),sampleInclusionCats);
sampleType = setcats(repmat(categorical(missing),height(cycles.metadata),1),sampleTypeCats);

% Assign Sample Inclusion Category to Each Cycle
sampleInclusion(ismember(cycles.metadata.msDatetime,rawDataset.metadata.msDatetime))='Included';
sampleInclusion(~ismember(cycles.metadata.msDatetime,rawDataset.metadata.msDatetime))='Not Included';
sampleInclusion(contains(cycles.metadata.ID1,'PIS'))='PIS';

isIncluded = sampleInclusion=='Included' | sampleInclusion=='PIS';

% Assign Sample Type Category to Each Cycle
% Order of these definitions is very important as successive definitions
% overwrite any previous definition for a sample. For example, the 'SPC'
% assignment includes all of the cans and LJA samples in the SPICE folder.
% This is then overwritten for those samples by the 'Std Can' and 'LJA'
% assignments. Also, the 'Jeff Neon' assignment at the end ensures that his
% Ne-Cans are not included with the rest of the Std Can data.
sampleType(contains(cycles.metadata.fileRelPath,'SPICE','IgnoreCase',true))='SPC'; % 'SPICE' is a folder name so this line will also capture all the LJA and Cans etc. in the folder
sampleType(contains(cycles.metadata.fileRelPath,'WDC','IgnoreCase',true) | contains(cycles.metadata.fileRelPath,'WAIS','IgnoreCase',true))='WDC';
sampleType(contains(cycles.metadata.fileRelPath,'Kenji') & contains(cycles.metadata.fileNameUser,'WDC')) = 'WDC';
sampleType(contains(cycles.metadata.fileRelPath,'NEEM','IgnoreCase',true)) = 'NEEM';
sampleType(contains(cycles.metadata.fileRelPath,'BYRD','IgnoreCase',true)) = 'BYRD';
sampleType(contains(cycles.metadata.fileRelPath,'SDMA','IgnoreCase',true)) = 'SDMA';
sampleType(contains(cycles.metadata.fileRelPath,'GISP','IgnoreCase',true)) = 'GISP';
sampleType(contains(cycles.metadata.fileRelPath,'Kenji','IgnoreCase',true) & contains(cycles.metadata.fileNameUser,'GISP','IgnoreCase',true)) ='GISP';
sampleType(contains(cycles.metadata.fileRelPath,'Daniel','IgnoreCase',true) & contains(cycles.metadata.fileNameUser,'GISP','IgnoreCase',true)) ='GISP';
sampleType(contains(cycles.metadata.fileRelPath,'Melissa','IgnoreCase',true) & contains(cycles.metadata.fileNameUser,'GISP','IgnoreCase',true)) ='GISP';
sampleType(contains(cycles.metadata.fileRelPath,'-ACQ-Results','IgnoreCase',true) & contains(cycles.metadata.fileNameUser,'GISP','IgnoreCase',true)) ='GISP';
sampleType(contains(cycles.metadata.fileRelPath,'Melissa','IgnoreCase',true) & contains(cycles.metadata.fileNameUser,'Vostok','IgnoreCase',true)) ='Vostok';
sampleType(contains(cycles.metadata.fileRelPath,'Kenji','IgnoreCase',true) & contains(cycles.metadata.fileNameUser,'DF','IgnoreCase',true)) ='Dome Fuji';
sampleType(contains(cycles.metadata.fileRelPath,'RICE','IgnoreCase',true))='RICE';

sampleType(contains(cycles.metadata.fileRelPath,'Bruce','IgnoreCase',true))='Bruce Plateau';
sampleType(contains(cycles.metadata.fileRelPath,'Taylor','IgnoreCase',true) | contains(cycles.metadata.fileRelPath,'TG','IgnoreCase',true))='Taylor Glacier';
sampleType(contains(cycles.metadata.fileRelPath,'Law-Dome','IgnoreCase',true))='Law Dome';
sampleType(contains(cycles.metadata.fileRelPath,'Minna-Bluff','IgnoreCase',true))='Minna Bluff';
sampleType(contains(cycles.metadata.fileRelPath,'Summit','IgnoreCase',true))='Summit';
sampleType(contains(cycles.metadata.fileRelPath,'Melissa','IgnoreCase',true) & contains(cycles.metadata.fileNameUser,'Summit','IgnoreCase',true)) ='Summit';


sampleType(contains(cycles.metadata.fileNameUser,'LJA','IgnoreCase',true) | contains(cycles.metadata.fileNameUser,'Pier','IgnoreCase',true))='LJA';
sampleType(contains(cycles.metadata.fileNameUser,'LJ Air','IgnoreCase',true) | contains(cycles.metadata.fileNameUser,'La Jolla Air','IgnoreCase',true))='LJA';
sampleType(contains(cycles.metadata.fileNameUser,'CS','IgnoreCase',true) | contains(cycles.metadata.fileNameUser,'Chem','IgnoreCase',true) | contains(cycles.metadata.fileNameUser,'Slope','IgnoreCase',true))='Chem Slope';
sampleType(contains(cycles.metadata.fileNameUser,'Aair','IgnoreCase',true) | contains(cycles.metadata.fileNameUser,'Can','IgnoreCase',true))='Std Can';
sampleType(count(cycles.metadata.fileNameUser,'air','IgnoreCase',true)==2 & contains(cycles.metadata.fileNameUser,'vs','IgnoreCase',true)) = 'Std Can';

sampleType(contains(cycles.metadata.fileNameUser,'Test','IgnoreCase',true) | contains(cycles.metadata.fileRelPath,'Test','IgnoreCase',true))='Misc. Test';
sampleType(contains(cycles.metadata.fileNameUser,'Test','IgnoreCase',true) & contains(cycles.metadata.fileNameUser,'Automation','IgnoreCase',true))='Automation Test';
sampleType(isundefined(sampleType) & contains(cycles.metadata.fileRelPath,'Jeff_Results','IgnoreCase',true))='Jeff Misc.';
sampleType(isundefined(sampleType) & contains(cycles.metadata.fileRelPath,'Keeling-Tanks','IgnoreCase',true))='Keeling Tanks Misc.';
sampleType(isundefined(sampleType) & contains(cycles.metadata.fileRelPath,'Sarah','IgnoreCase',true))='Sarah Misc.';
sampleType(isundefined(sampleType) & contains(cycles.metadata.fileRelPath,'Alan','IgnoreCase',true))='Alan Misc.';
sampleType(isundefined(sampleType) & contains(cycles.metadata.fileRelPath,'Kenji-Intercomparison','IgnoreCase',true))='Kenji Misc.';
sampleType(isundefined(sampleType) & contains(cycles.metadata.fileRelPath,'Vas','IgnoreCase',true))='Vas Misc.';
sampleType(isundefined(sampleType) & contains(cycles.metadata.fileRelPath,'Wendy','IgnoreCase',true))='Wendy Misc.';
sampleType(isundefined(sampleType) & contains(cycles.metadata.fileRelPath,'Jinho','IgnoreCase',true))='Jinho Misc.';
sampleType(contains(cycles.metadata.fileNameUser,'Neon','IgnoreCase',true) | contains(cycles.metadata.fileRelPath,'Neon','IgnoreCase',true) | (contains(cycles.metadata.fileRelPath,'ACQ','IgnoreCase',true) & cycles.metadata.msDatetime > datetime(2018,1,1)))='Jeff Neon';

idxUndef = find(isundefined(sampleType));

%% Plot Timeline
% Create the timeline showing when each sample type was analyzed and
% whether or not it is included in the 'standard' raw dataset.

% Prepare Plotting Variables for the for-loop
sampleTypeToPlot = mergecats(sampleType,["Jeff Misc." "Keeling Tanks Misc." "Sarah Misc." "Alan Misc." "Kenji Misc." "Vas Misc." "Wendy Misc." "Jinho Misc." "Jeff Neon"],"Misc. Samples");
sampleTypeCatsToPlot = string(categories(sampleTypeToPlot));
yVal = [
    0.2:0.05:0.3 ... % Y Values for tests and misc samples
    0.4:0.05:0.5 ... % Y Values for std can, CS, and LJA
    0.6:0.05:1.25 ... % Y Values for Ice Cores
    ];
cols = [
    lineCol(9)*0.7; lineCol(9)*0.7; lineCol(9)*0.7; ... % Colors for tests and misc samples
    lineCol(3); lineCol(8); lineCol(6); ... % Colors for std can, CS, and LJA
    lineCol(2); lineCol(3); lineCol(1); lineCol(7); lineCol(1)*0.7; lineCol(2)*0.7; lineCol(7)*0.7; lineCol(8)*0.85; lineCol(2)*1.1; lineCol(5); lineCol(4); lineCol(1)*1.3; lineCol(4)*1.3; lineCol(5)*0.7 ... % Colors for Ice Cores
    ];

% Plot a Few Extra Categories not in the Plotting Array
figure; hold on;
plot(cycles.metadata.msDatenum(~isIncluded),zeros(sum(~isIncluded),1),'^k','MarkerFaceColor','k','MarkerSize',3)

plot(cycles.metadata.msDatenum(isIncluded & isundefined(sampleType)),0.1*ones(sum(isIncluded & isundefined(sampleType)),1),'.','Color',lineCol(9),'MarkerSize',20)
plot(cycles.metadata.msDatenum(~isIncluded & isundefined(sampleType)),0.1*ones(sum(~isIncluded & isundefined(sampleType)),1),'x','Color',lineCol(9)*0.7,'MarkerSize',5)

% Plot the Data in the Categorical Array
for ii = 1:length(sampleTypeCatsToPlot)
    plot(cycles.metadata.msDatenum(isIncluded & sampleTypeToPlot==sampleTypeCatsToPlot(ii)),yVal(ii)*ones(sum(isIncluded & sampleTypeToPlot==sampleTypeCatsToPlot(ii)),1),'.','Color',cols(ii,:),'MarkerSize',20)
    plot(cycles.metadata.msDatenum(~isIncluded & sampleTypeToPlot==sampleTypeCatsToPlot(ii)),yVal(ii)*ones(sum(~isIncluded & sampleTypeToPlot==sampleTypeCatsToPlot(ii)),1),'x','Color',cols(ii,:)*0.7,'MarkerSize',5)
end

ax=gca;
ax.YGrid = 'on';
ax.YTick = [0 0.1 yVal];
ax.YTickLabel = ["Excluded Samples" "Unassigned" sampleTypeCatsToPlot'];
ax.XLim = datenum([2006 2021],[01 01],[01 01]);
ax.XTick = datenum(2006:2020,ones(1,1,15),ones(1,1,15));
ax.XTickLabel = string(2006:2020);
ax.YLim = ax.YLim.*[1 1.05];
ax.YLimMode = 'manual';

% Plot the Date Ranges that Do Not Fall within any of the Known Spreadsheets
[spreadsheetRanges,missingIntervals] = getSpreadsheetRanges('spreadsheet_metadata.xlsx');
missingIntervals = [missingIntervals; max(spreadsheetRanges.EndDate) datetime('now')];

xvals = nan(4,length(missingIntervals));
xvals([1 2],1:end) = repmat(datenum(missingIntervals(:,1))',2,1);
xvals([3 4],1:end) = repmat(datenum(missingIntervals(:,2))',2,1);

yvals = nan(4,length(missingIntervals));
yvals([1 4],1:end) = ax.YLim(1);
yvals([2 3],1:end) = ax.YLim(2);

patch(xvals,yvals,[1 0.7 0.7])

ax.Children = flipud(ax.Children);
% xlim(datenum([2006 2020],[01 01],[01 01]))
% datetick('x','KeepLimits')


% Plot Mass Spec Events
massSpecEvents = readtable('spreadsheet_metadata.xlsx','Sheet',1);
massSpecEvents.Event = categorical(massSpecEvents.Event);

iPlot = massSpecEvents.Event == "New Filament" | massSpecEvents.Event == "Refocus";
plot(datenum([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)']),repmat(get(gca,'YLim')',1,sum(iPlot)),'-r','LineWidth',2);
text(datenum(massSpecEvents.EndDate(iPlot)),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color','r')

iPlot = massSpecEvents.Event == "New Std Cans" | massSpecEvents.Event == "Swap Std Cans";
plot(datenum([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)']),repmat(get(gca,'YLim')',1,sum(iPlot)),':','Color',lineCol(6),'LineWidth',3);
text(datenum(massSpecEvents.EndDate(iPlot)),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color',lineCol(6))

% xlim(datetime([2014 2019],[1 1],[1 1]))