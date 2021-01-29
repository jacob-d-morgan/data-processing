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
exp2match = '^(?<bagNum>\d{3,4})\s?(?<rep>\w)?';
tokens = regexpi(neemMasterSheet.metadata.fileNameUser(:,1,1),exp2match,'tokens');
% This expression captures all regular samples. The first token ('NEEM\s') 
% is empty for all but 4 samples and the Bag Number is the second token.
% The expression also captures a few BYRD samples, which are filtered out
% below using the fileRelPath.

exp2match = '^(?:NEEM\s)?(?<bagNum>X\s*\d{3,4})\s?(?<rep>\w)?';
tokensX = regexpi(neemMasterSheet.metadata.fileNameUser(:,1,1),exp2match,'tokens');
% This expression captures the samples from bags labelled with an X. It is,
% as yet, unclear what depths these samples correspond to.

exp2match = '(?<bagNum>T\s\d{1,4})(\s(?<bagDepths>\d{1,2}-\d{1,2}))?';
tokensT = regexp(neemMasterSheet.metadata.fileNameUser(:,1,1),exp2match,'tokens');
% This expression captures the samples from the Bolling Transition, which
% have a different bag number convention and also include the depths of the
% sample within the bag separated by a hyphen.


% Identify NEEM Ice Samples
iRegular = contains(neemMasterSheet.metadata.fileRelPath(:,1,1),'NEEM') & ...
           ~cellfun(@isempty,tokens);
iX = contains(neemMasterSheet.metadata.fileRelPath(:,1,1),'NEEM') & ...
     ~cellfun(@isempty,tokensX);
iTransition = contains(neemMasterSheet.metadata.fileRelPath(:,1,1),'NEEM') & ...
              ~cellfun(@isempty,tokensT); 
iNEEM = iRegular | iX | iTransition;


% Combine Bag Numbers and Depth Ranges
bagNum = repmat("",size(tokens));
bagNum(iRegular) = cellfun(@(x) x{1}(1),tokens(iRegular));
bagNum(iX) = cellfun(@(x) x{1}(1),tokensX(iX));
bagNum(iTransition) = cellfun(@(x) x{1}(1),tokensT(iTransition));
bagNum = str2double(erase(bagNum,{'x','X','T',' '}));

bagDepthRange = repmat("",size(tokensT));
bagDepthRange(iTransition) = cellfun(@(x) x{1}(2),tokensT(iTransition));
bagTopDepths = nan(size(tokensT));
bagTopDepths(iTransition) = str2double(extractBefore(bagDepthRange(iTransition),'-'));
bagTopDepths(isnan(bagTopDepths) & iTransition) = 0; % Some transition samples have no depth range info included. Assume these samples were cut from the top (z = 0 cm).

% Assign Depths
% N.B. The Danes report Top Depths for NEEM that are just the sum of the
% number of 55 cm bags that have come before. The SIO piece seems to often
% comprise the top 5 or 10 cm of each bag (unless noted, e.g. for the 
% transition pieces). 
% The Transition samples were renumbered such that T1 = Bag 2580, which has
% a top depth of 1418.85 m (see BigSampleBook_16112011_bestguessTB.xls).
% The X-labelled samples were renumbered such that X228 = Bag 3471, which
% has a top depth of 1908.5 m (see BigSampleBook_16112011_bestguessTB.xls).
topDepth = nan(size(bagNum));
topDepth(iRegular) = (bagNum(iRegular)-1)*0.55;
topDepth(iX) = (bagNum(iX)-228)*0.55 + 1908.5;
topDepth(iTransition) = 1418.45 + (bagNum(iTransition)-1)*0.55 + bagTopDepths(iTransition)/100;


% Generate Ice Core Dataset
neem = makeIceCoreDataset(neemMasterSheet,iNEEM,topDepth(iNEEM));
neem.metadata.Properties.VariableNames{string(neem.metadata.Properties.VariableNames)=='depth'} = 'topDepth';

