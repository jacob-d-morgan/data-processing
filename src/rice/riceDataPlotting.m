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
% For RICE, I choose to use the ISODAT File Name because the User File Name
% is trimmed incorrectly by the code that makes a guess at the User File
% Name by removing line numbers (i.e. the d{1,3}_* pattern). 
% I could use ID1 instead but there are a few cases where the ID1 seems to
% not have been updated by the user (see 2015-03-26). Presumably the
% samples were run manually, rather than with the automation, so the user
% manually entered a new file name but forgot to clear the ID1 from the
% previous sample. There is one other sample (2015-02-03 21:38) that is
% missing an ID1, presumable also because it was run manually.

% Extract Bottom Depths from Metadata
exp2match = '^\s?(?<depth>\d{2,3}([\.|_]\d{1,2})?)(?<rep>\s\w)?_?-\d{4,4}$';
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
rice = makeIceCoreDataset(riceMasterSheet,iRICE & ~isnan(bottomDepth),bottomDepth(iRICE & ~isnan(bottomDepth)));
% Including the isnan clause filters out one problematic NaN value that
% arises from assigning a NaN depth to a Rep B that has no corresponding
% Rep A (so the bottom depth is NaN). The Rep A is NaN because it has two
% decimal places/underscores so is not captured by the regexp.
