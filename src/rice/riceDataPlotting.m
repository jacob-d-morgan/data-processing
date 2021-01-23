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

%%
% Generate Ice Core Dataset
iceCoreID = 'RICE';
exp2match = '^\s?(?<depth>\d{2,3}([\.|_]\d{1,2})?)(?<rep>\s\w)?-\d{4,4}$';
tokens = regexp(riceMasterSheet.metadata.fileNameIsodat(:,1),exp2match,'tokens');

exp2match = '^b-\d{4,4}$';
tokensRepB = regexp(riceMasterSheet.metadata.fileNameIsodat(:,1),exp2match);

bD = repmat("",size(tokens));
bD(~cellfun(@isempty,tokens)) = cellfun(@(x) x{1}(1),tokens(~cellfun(@isempty,tokens)));
bD(~cellfun(@isempty,tokensRepB)) = bD(find(~cellfun(@isempty,tokensRepB))-1);

iRICE = contains(riceMasterSheet.metadata.fileRelPath(:,1,1),iceCoreID) & bD~="";

bottomDepth = str2double(replace(bD,"_","."));
bottomDepth = bottomDepth(~isnan(bottomDepth));

rice = makeIceCoreDataset(riceMasterSheet,iRICE,bottomDepth);
