%% NEEM Data Plotting
% Plots figures from the NEEM (North Greenland Eemian Ice Drilling) XP
% dataset.

clearvars; cluk;
clc;

%% Import Data
% Load the relevant csv files and compile the data into the NEEM master
% sheet and ice core dataset.

% Generate Master Sheet
filesToImport = [
    "XP-2009(exportJDM).csv";...
    "XP-2010(exportJDM).csv";...
    ];

neemMasterSheet = makeMasterSheet(filesToImport);
