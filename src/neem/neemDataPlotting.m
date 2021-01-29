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

%% Generate Ice Core Dataset
% Identify the ice samples within the relevant folder and extract the
% bottom depth from the metadata to use as an input to makeIceCoreDataset.

% Extract Bottom Depths from Metadata
exp2match = '^(NEEM\s)?(?<bagNum>(X\s*)?\d{3,4})\s?(?<rep>\w)?';
tokens = regexpi(neemMasterSheet.metadata.fileNameUser(:,1,1),exp2match,'tokens');
% This expression captures all regular samples. The first token ('NEEM\s') 
% is empty for all but 4 samples and the Bag Number is the second token.
% The expression also captures a few BYRD samples, which are filtered out
% below using the fileRelPath.

exp2match = '(?<bagNum>T\s\d{1,4})(\s(?<bagDepths>\d{1,2}-\d{1,2}))?';
tokensT = regexp(neemMasterSheet.metadata.fileNameUser(:,1,1),exp2match,'tokens');
% This expression captures the samples from the Bolling Transition, which
% have a different bag number convention and also include the depths of the
% sample within the bag separated by a hyphen.


% Identify NEEM Ice Samples
iRegular = contains(neemMasterSheet.metadata.fileRelPath(:,1,1),'NEEM') & ...
        ~cellfun(@isempty,tokens);
iTransition = contains(neemMasterSheet.metadata.fileRelPath(:,1,1),'NEEM') & ...
               ~cellfun(@isempty,tokensT); 
iNEEM = iRegular | iTransition;


% Combine Bag Numbers and Depths
bagNum = repmat("",size(tokens));
bagNum(iRegular) = cellfun(@(x) x{1}(2),tokens(iRegular));
bagNum(iTransition) = cellfun(@(x) x{1}(1),tokensT(iTransition));
bagNum = str2double(erase(bagNum,{'x','X','T',' '}));

bagDepthRange = repmat("",size(tokensT));
bagDepthRange(iTransition) = cellfun(@(x) x{1}(2),tokensT(iTransition));
bagTopDepths = nan(size(tokensT));
bagTopDepths(iTransition) = str2double(extractBefore(bagDepthRange(iTransition),'-'));
bagTopDepths(isnan(bagTopDepths) & iTransition) = 0;

% Assign Depths
% N.B. The Danes report Top Depths for NEEM that are just the sum of the
% number of 55 cm bags that have come before. The SIO piece seems to always
% comprise the top 5 or 10 cm of each bag (unless noted, e.g. for the 
% transition pieces). The bags for the Transition samples were renumbered
% starting from bag 2580, which has a top depth of 1418.85 m.
topDepth = nan(size(bagNum));
topDepth(iRegular) = (bagNum(iRegular)-1)*0.55;
topDepth(iTransition) = 1418.45 + (bagNum(iTransition)-1)*0.55 + bagTopDepths(iTransition)/100;
topDepth = topDepth(iNEEM);


% Generate Ice Core Dataset
neem = makeIceCoreDataset(neemMasterSheet,iNEEM,topDepth);