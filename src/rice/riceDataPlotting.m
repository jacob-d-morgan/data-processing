%% RICE Data Plotting
% Plots figures from the RICE (Roosevelt Island Climate Experiment) XP
% dataset.

clearvars; cluk;
clc;

%% Import Data
% Load the relevant csv files and compile the data into the RICE master
% sheet and ice core dataset.

% Generate Master Sheet
filesToImport = [
    "XP-2014(exportJDM).csv";...
    "XP-2015(exportJDM).csv";...
    ];

riceMasterSheet = makeMasterSheet(filesToImport);

%% Generate Ice Core Dataset
% Identify the ice samples within the relevant folder and extract the
% bottom depth from the metadata to use as an input to makeIceCoreDataset.

% Extract Bottom Depths from Metadata
exp2match = '^\s?(?<depth>\d{2,3}([\.|_]\d{1,2})?)(?<rep>\s\w)?-\d{4,4}$';
tokens = regexp(riceMasterSheet.metadata.fileNameIsodat(:,1),exp2match,'tokens');
% This expression matches all of the 'Rep A' samples and also some of the
% 'Rep B' samples that are labelled with bottom depths. Many of the 'Rep B'
% samples are not labelled with bottom depth. These are captured separately
% using the expression below.

exp2match = '^b-\d{4,4}$';
tokensRepB = regexp(riceMasterSheet.metadata.fileNameIsodat(:,1),exp2match);
% This expression matches the 'Rep B' samples that are not labelled with
% any bottom depth information. The bottom depth is the same as the
% previous sample.

% Identify RICE Ice Samples
iRICE = contains(riceMasterSheet.metadata.fileRelPath(:,1,1),'RICE') & ...
        (~cellfun(@isempty,tokens) | ~cellfun(@isempty,tokensRepB));

% Calculate Bottom Depth
bD_temp = repmat("",size(tokens));
bD_temp(~cellfun(@isempty,tokens)) = cellfun(@(x) x{1}(1),tokens(~cellfun(@isempty,tokens)));
bD_temp(~cellfun(@isempty,tokensRepB)) = bD_temp(find(~cellfun(@isempty,tokensRepB))-1);
bottomDepth = str2double(replace(bD_temp,"_","."));

% Generate Ice Core Dataset
rice = makeIceCoreDataset(riceMasterSheet,iRICE,bottomDepth);
