%% SDMA Data Plotting
% Plots figures from the SDMA (Siple Dome A) XP dataset.

clearvars; cluk;
clc;

%% Import Data
% Load the relevant csv files and compile the data into the NEEM master
% sheet and ice core dataset.

% Generate Master Sheet
filesToImport = [
    "XP-2006(exportJDM).csv";...
    "XP-2007(exportJDM).csv";...
    ];

sdmaMasterSheet = makeMasterSheet(filesToImport);

%% Generate Ice Core Dataset
% Identify the ice samples within the relevant folder and extract the
% bottom depth from the metadata to use as an input to makeIceCoreDataset.

% Extract Bottom Depths from Metadata
exp2match = 'SDMA\s*(?<depth>\d{2,3}\.\d{0,3})\s*?(?<rep>[a-zA-Z])?';
tokens = regexp(sdmaMasterSheet.metadata.fileNameUser(:,1,1),exp2match,'tokens');

iSDMA = contains(sdmaMasterSheet.metadata.fileRelPath(:,1,1),'SDMA') & ~cellfun(@isempty,tokens);

bottomDepth = str2double(cellfun(@(x) x{1}(1),tokens(iSDMA)));

sdma = makeIceCoreDataset(sdmaMasterSheet,iSDMA,bottomDepth);
sdma.metadata.Properties.VariableNames{string(sdma.metadata.Properties.VariableNames)=='depth'} = 'bottomDepth';