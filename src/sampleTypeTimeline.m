%% Sample Type Timeline
% Plots a 'timeline' of the different categories of the imported data. This
% allows me to identify the different batches of CS, LJA, and which ice
% cores they are applied to.
% It also allows me to manually plot the date of significant changes to the
% MS, such as a refocussing, or significant adjustments to the method, such
% as the implementation of the O2-consumption correction.

% clearvars;
cluk; clc;

%% Import Data
% Loads in the data from the .csv files and calculates  the delta values
% andmetadata. Also generates the standard raw dataset in order to
% determine which cycles are omitted. Can comment out the csvReadIsodat,
% calcCycles, and makeRawDataset lines, and use the clearvars -EXCEPT to
% save time if re-running a lot.

clearvars -EXCEPT table_deltas cycle_metadata data metadata ;

% Load in data from .csv files
% importedData = csvReadIsodat([
%     "XP-2018(excelExportIntensityJDM).csv"; ...
%     "XP-2017(excelExportIntensityJDM).csv"; ...
%     "XP-2016(excelExportIntensityJDM).csv"; ...
%     "XP-2015(excelExportIntensityJDM-REMAKE).csv"; ...
%     "XP-2014(excelExportIntensityJDM-REMAKE).csv"; ...
%     "XP-2013(excelExportIntensityJDM-REMAKE).csv"; ...
%     ]);
% [table_deltas,cycle_metadata] = calcCycles(importedData,'IsRef__');

delta_names = string(table_deltas.Properties.VariableNames);
delta_labels = string(table_deltas.Properties.VariableDescriptions);
delta_units = string(table_deltas.Properties.VariableUnits);
cycle_deltas = table2array(table_deltas);

% Generate raw dataset
% [data,metadata] = makeRawDataset([
%     "XP-2018(excelExportIntensityJDM).csv"; ...
%     "XP-2017(excelExportIntensityJDM).csv"; ...
%     "XP-2016(excelExportIntensityJDM).csv"; ...
%     "XP-2015(excelExportIntensityJDM-REMAKE).csv"; ...
%     "XP-2014(excelExportIntensityJDM-REMAKE).csv"; ...
%     "XP-2013(excelExportIntensityJDM-REMAKE).csv"; ...
%     ]);

aliquot_deltas = data;
aliquot_metadata = metadata.metadata;
delta_names = metadata.delta_names;
delta_labels = metadata.delta_labels;
delta_units = metadata.delta_units;

metadata_fields = string(fieldnames(aliquot_metadata))';



%% Define and Assign Categories
% Assign a sample type category to each sample and also determine whether
% or not it is included in the calculated cycle deltas. Samples that are
% not included are those with a gas configuration different to Air+ or
% where the 28N2 beam is saturated or absent.

% Define Categories
sampleInclusionCats = {'Included','Not Included','PIS'};
sampleTypeCats = {'Automation Test','Misc. Test','Sarah Misc.','Alan Misc.','Jinho Misc.','Jeff Neon','Std Can','Chem Slope','LJA','SPC','WDC','GISP','RICE','Bruce Plateau','Taylor Glacier','Summit'};

% Make Empty Categorical Arrays to Fill
sampleInclusion = setcats(repmat(categorical(missing),length(cycle_metadata.filename),1),sampleInclusionCats);
sampleType = setcats(repmat(categorical(missing),length(cycle_metadata.filename),1),sampleTypeCats);

% Assign Sample Inclusion Category
sampleInclusion(ismember(cycle_metadata.msDatetime,aliquot_metadata.msDatetime))='Included';
sampleInclusion(~ismember(cycle_metadata.msDatetime,aliquot_metadata.msDatetime))='Not Included';
sampleInclusion(contains(cycle_metadata.ID1,'PIS'))='PIS';

isIncluded = sampleInclusion=='Included' | sampleInclusion=='PIS';

% Assign Sample Type Category
% Order of these definitions is very important as each successive
% definition overwrites any previous definition for a sample. For example,
% the 'SPC' assignment includes all of the cans and LJA samples in the
% SPICE folder. This is then overwritten by the 'Std Can' and 'LJA'
% assignments for the relevant samples. Also, the 'Jeff Neon' assignment at
% the end ensures that his Ne-Cans are not included with the rest of the
% Std Can data.
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
% Create the timeline showing when each sample type was analyzed.

% Prepare plotting variables for the for-loop
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

% Plot a few extra categories that aren't in the plotting array
figure; hold on;
plot(cycle_metadata.msDatetime(~isIncluded),zeros(sum(~isIncluded),1),'^k','MarkerFaceColor','k')

plot(cycle_metadata.msDatetime(isIncluded & isundefined(sampleType)),0.1*ones(sum(isIncluded & isundefined(sampleType)),1),'.','Color',lineCol(9),'MarkerSize',20)
plot(cycle_metadata.msDatetime(~isIncluded & isundefined(sampleType)),0.1*ones(sum(~isIncluded & isundefined(sampleType)),1),'x','Color',lineCol(9)*0.7,'MarkerSize',5)

% Plot the categorical array
for ii = 1:length(sampleTypeCatsToPlot)
    plot(cycle_metadata.msDatetime(isIncluded & sampleType==sampleTypeCatsToPlot(ii)),yVal(ii)*ones(sum(isIncluded & sampleType==sampleTypeCatsToPlot(ii)),1),'.','Color',cols(ii,:),'MarkerSize',20)
    plot(cycle_metadata.msDatetime(~isIncluded & sampleType==sampleTypeCatsToPlot(ii)),yVal(ii)*ones(sum(~isIncluded & sampleType==sampleTypeCatsToPlot(ii)),1),'x','Color',cols(ii,:)*0.7,'MarkerSize',5)
end

ax=gca;
ax.YTick = [0 0.1 yVal];
ax.YTickLabel = ["Excluded Samples" "Unassigned" sampleTypeCatsToPlot'];
ax.XLim = datetime([2013 2019],[01 01],[01 01]);

