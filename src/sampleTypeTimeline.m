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

clearvars -EXCEPT table_deltas cycle_metadata aliquot_deltas metadata ;

filesToImport = [
    "XP-2018(excelExportIntensityJDM).csv"; ...
    "XP-2017(excelExportIntensityJDM).csv"; ...
    "XP-2016(excelExportIntensityJDM).csv"; ...
    "XP-2015(excelExportIntensityJDM-REMAKE).csv"; ...
    "XP-2014(excelExportIntensityJDM-REMAKE).csv"; ...
    "XP-2013(excelExportIntensityJDM-REMAKE).csv"; ...
    ];

% Load in Data from .csv Files
% importedData = csvReadIsodat(filesToImport);                             % <-- TOGGLE MY COMMENT
% [table_deltas,cycle_metadata] = calcCycles(importedData,'IsRef__');      % <-- TOGGLE MY COMMENT

cycle_deltas = table2array(table_deltas);

% Generate 'Standard' Raw Dataset
% [aliquot_deltas,metadata] = makeRawDataset(filesToImport);                % <-- TOGGLE MY COMMENT

aliquot_metadata = metadata.metadata;
delta_names = metadata.delta_names;
delta_labels = metadata.delta_labels;
delta_units = metadata.delta_units;

metadata_fields = string(fieldnames(aliquot_metadata))';


%% Define and Assign Categories
% Assign a sample type category to each sample and also determine whether
% or not it is included in the 'standard' raw dataset. Samples from
% cycle_deltas that are not included are those with a gas configuration
% different to Air+ or where the 28N2 beam is saturated or absent and those
% with a non-standard number of cycles or blocks.

% Define Categories
sampleInclusionCats = {'Included','Not Included','PIS'};
sampleTypeCats = {'Automation Test','Misc. Test','Sarah Misc.','Alan Misc.','Jinho Misc.','Jeff Neon','Std Can','Chem Slope','LJA','SPC','WDC','GISP','RICE','Bruce Plateau','Taylor Glacier','Summit'};

% Make Empty Categorical Arrays to Fill
sampleInclusion = setcats(repmat(categorical(missing),length(cycle_metadata.filename),1),sampleInclusionCats);
sampleType = setcats(repmat(categorical(missing),length(cycle_metadata.filename),1),sampleTypeCats);

% Assign Sample Inclusion Category to Each Cycle
sampleInclusion(ismember(cycle_metadata.msDatetime,aliquot_metadata.msDatetime))='Included';
sampleInclusion(~ismember(cycle_metadata.msDatetime,aliquot_metadata.msDatetime))='Not Included';
sampleInclusion(contains(cycle_metadata.ID1,'PIS'))='PIS';

isIncluded = sampleInclusion=='Included' | sampleInclusion=='PIS';

% Assign Sample Type Category to Each Cycle
% Order of these definitions is very important as successive definitions
% overwrite any previous definition for a sample. For example, the 'SPC'
% assignment includes all of the cans and LJA samples in the SPICE folder.
% This is then overwritten for those samples by the 'Std Can' and 'LJA'
% assignments. Also, the 'Jeff Neon' assignment at the end ensures that his
% Ne-Cans are not included with the rest of the Std Can data.
sampleType(contains(cycle_metadata.filename,'SPICE','IgnoreCase',true))='SPC'; % 'SPICE' is a folder name so this line will also capture all the LJA and Cans etc. in the folder
sampleType(contains(cycle_metadata.filename,'WDC','IgnoreCase',true))='WDC'; % 'WDC' is a folder name so this line will also capture all the LJA and Cans etc. in the folder
sampleType(contains(cycle_metadata.filename,'GISP','IgnoreCase',true))='GISP';
sampleType(contains(cycle_metadata.filename,'RICE','IgnoreCase',true))='RICE';
sampleType(contains(cycle_metadata.filename,'Bruce','IgnoreCase',true))='Bruce Plateau';
sampleType(contains(cycle_metadata.filename,'Taylor','IgnoreCase',true) | contains(cycle_metadata.filename,'TG','IgnoreCase',true))='Taylor Glacier';
sampleType(contains(cycle_metadata.filename,'Summit','IgnoreCase',true))='Summit';

sampleType(contains(cycle_metadata.filename,'LJA','IgnoreCase',true) | contains(cycle_metadata.filename,'Pier','IgnoreCase',true))='LJA';
sampleType(contains(cycle_metadata.filename,'CS','IgnoreCase',true) | contains(cycle_metadata.filename,'Chem','IgnoreCase',true) | contains(cycle_metadata.filename,'Slope','IgnoreCase',true))='Chem Slope';
sampleType(contains(cycle_metadata.filename,'Aair','IgnoreCase',true) | contains(cycle_metadata.filename,'Can','IgnoreCase',true))='Std Can';

sampleType(contains(cycle_metadata.filename,'Test','IgnoreCase',true))='Misc. Test';
sampleType(contains(cycle_metadata.filename,'Test','IgnoreCase',true) & contains(cycle_metadata.filename,'Automation','IgnoreCase',true))='Automation Test';
sampleType(isundefined(sampleType) & contains(cycle_metadata.filename,'Sarah','IgnoreCase',true))='Sarah Misc.';
sampleType(isundefined(sampleType) & contains(cycle_metadata.filename,'Alan','IgnoreCase',true))='Alan Misc.';
sampleType(isundefined(sampleType) & contains(cycle_metadata.filename,'Jinho','IgnoreCase',true))='Jinho Misc.';
sampleType(contains(cycle_metadata.filename,'Neon','IgnoreCase',true) | (contains(cycle_metadata.filename,'ACQ','IgnoreCase',true) & cycle_metadata.msDatetime > datetime(2018,1,1)))='Jeff Neon';

idxUndef = find(isundefined(sampleType));

%% Plot Timeline
% Create the timeline showing when each sample type was analyzed and
% whether or not it is included in the 'standard' raw dataset.

% Prepare Plotting Variables for the for-loop
sampleTypeToPlot = mergecats(sampleType,["Sarah Misc." "Alan Misc." "Jinho Misc." "Jeff Neon"],"Misc. Samples");
sampleTypeCatsToPlot = string(categories(sampleTypeToPlot));
yVal = [
    0.2:0.05:0.3 ... % Y Values for tests and misc samples
    0.4:0.05:0.5 ... % Y Values for std can, CS, and LJA
    0.6:0.05:0.9 ... % Y Values for Ice Cores
    ];
cols = [
    lineCol(9)*0.7; lineCol(9)*0.7; lineCol(9)*0.7; ... % Colors for tests and misc samples
    lineCol(3); lineCol(8); lineCol(6); ... % Colors for std can, CS, and LJA
    lineCol(2); lineCol(1); lineCol(7); lineCol(1)*0.7; lineCol(5); lineCol(4); lineCol(5)*0.7 ... % Colors for Ice Cores
    ];

% Plot a Few Extra Categories not in the Plotting Array
figure; hold on;
plot(cycle_metadata.msDatenum(~isIncluded),zeros(sum(~isIncluded),1),'^k','MarkerFaceColor','k','MarkerSize',3)

plot(cycle_metadata.msDatenum(isIncluded & isundefined(sampleType)),0.1*ones(sum(isIncluded & isundefined(sampleType)),1),'.','Color',lineCol(9),'MarkerSize',20)
plot(cycle_metadata.msDatenum(~isIncluded & isundefined(sampleType)),0.1*ones(sum(~isIncluded & isundefined(sampleType)),1),'x','Color',lineCol(9)*0.7,'MarkerSize',5)

% Plot the Data in the Categorical Array
for ii = 1:length(sampleTypeCatsToPlot)
    plot(cycle_metadata.msDatenum(isIncluded & sampleTypeToPlot==sampleTypeCatsToPlot(ii)),yVal(ii)*ones(sum(isIncluded & sampleTypeToPlot==sampleTypeCatsToPlot(ii)),1),'.','Color',cols(ii,:),'MarkerSize',20)
    plot(cycle_metadata.msDatenum(~isIncluded & sampleTypeToPlot==sampleTypeCatsToPlot(ii)),yVal(ii)*ones(sum(~isIncluded & sampleTypeToPlot==sampleTypeCatsToPlot(ii)),1),'x','Color',cols(ii,:)*0.7,'MarkerSize',5)
end

ax=gca;
ax.YGrid = 'on';
ax.YTick = [0 0.1 yVal];
ax.YTickLabel = ["Excluded Samples" "Unassigned" sampleTypeCatsToPlot'];
ax.XLim = datenum([2013 2019],[01 01],[01 01]);

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
xlim(datenum([2013 2019],[01 01],[01 01]))
datetick('x','KeepLimits')
