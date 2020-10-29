function [rawDataset, varargout] = makeRawDataset(filesToImport, varargin)
% MAKERAWDATASET generates the delta values and metadata from FILES
%   RAWDATASET = makeRawDataset(filesToImport) imports the raw
%   cycle integrations from the .csv files in FILESTOIMPORT and performs
%   several manipulations to output the raw dataset as a structure.
%
%   FILESTOIMPORT must be a string array of filenames or paths to files
%   that MATLAB can access.
%   
%   RAWDATASET is a structure with two fields. The 'metadata' field is a
%   table of metadata arranged as aliquot-by-variable-by-block-by-cycle.
%   The 'deltas' field is a table of delta values arranged as
%   aliquot-by-delta value-by-block-by-cycle.
%
%   The manipulations performed to produce RAWDATASET are as follows:
%       1. Import the data from the .csv files in FILESTOIMPORT using
%       helper function csvReadIsodat.m
%       2. Calculate cycle delta values using helper function calcCycles.m
%       3. Reject cycles with a gas configuration that we are not
%       interested in or where the m/z = 28 beam is saturated or absent.
%       4. Reshape the array of cycle deltas and metadata to the dimensions
%       of the output (described above). This step excludes blocks with
%       more or fewer cycles than the mode, and aliquots with more or fewer
%       blocks than the mode. These can be optionally outputted in separate
%       variables if desired (see below)
%
%   Optionally, additional outputs can be provided by including one or more
%   of the following flags, and by including the appropriate number of
%   output variabes:
%    - [...,RAWDATASETPIS] = makeRawDataset(...,'includePIS')
%      includes delta values and their metadata from the PIS experiment
%      blocks in the same shape and structure as the primary outputs.
%    - [...,RAWDATASETALLBLOCKS] = makeRawDataset(...,'includeAllBlocks')
%      outputs a cell array containing the aliquots composed of the mode
%      number of blocks of any number of cycles.
%    - [...,RAWDATASETALLALIQUOTS] = makeRawDataset(...,'includeAllAliquots')
%      outputs a cell array containing aliquots of any number of blocks of
%      the mode number of cycles.
%    - [...,RAWDATASETALLDATA] = makeRawDataset(...,'includeAllData')
%      outputs a cell array containing aliquots of any number of blocks of
%      any number of cycles.
%
% The optional flags can be combined to output both the PIS data and the
% variables containing all the blocks/aliquots:
%    - [...,RAWDATASAETPIS,RAWDATASETALLDATE] = makeRawDataset(...,'includePIS','includeAllData')
%      outputs both the PIS data and a cell array containing aliquots of 
%      any number of blocks of any number of cycles.
%
% -------------------------------------------------------------------------

%% Parse Inputs
% Check the user-provided inputs are of the correct form.

narginchk(1,inf)

% Define Input Parser Scheme
p = inputParser;
p.FunctionName = mfilename;

addRequired(p,'filesToImport');
addParameter(p,'includePIS',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'includeAllBlocks',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'includeAllAliquots',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'includeAllData',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

% Parse Inputs and Assign Results
     parse(p,filesToImport,varargin{:});

filesToImport = p.Results.filesToImport;
includePIS = p.Results.includePIS;
includeAllBlocks = p.Results.includeAllBlocks;
includeAllAliquots = p.Results.includeAllAliquots;
if p.Results.includeAllData
    includeAllBlocks = true;
    includeAllAliquots = true;
end

% Check for Correct Number of Outputs
if includePIS && (includeAllBlocks || includeAllAliquots)
    nargoutchk(3,3);
elseif includePIS && ~(includeAllBlocks || includeAllAliquots)
    nargoutchk(2,2)
elseif ~includePIS && (includeAllBlocks || includeAllAliquots)
    nargoutchk(2,2)
else
    nargoutchk(1,1);
end

clear p


%% Import Data, Calculate Delta Values, and Compile Useful Metadata
% NOTE: Interpolating onto the ~isRef indices excludes the CO2 checks from
% the blocks I consider as they are recorded in isRef as 'true'. I have to
% manually add them back in below, including their metadata.

importedData = csvReadIsodat(filesToImport);

cycles = calcCycles(importedData,'IsRef__');

%% Reject Erroneous Cycles
% Get rid of some data that is of little use for evaluating the health of
% the mass spec and which causes problems with PIS etc.
%
% First, discard all data with a gas config different to 'Air+'
% This is mostly Jeff's neon experiments.
    % N.B. Note that there is a difference between Gas Config and Gas Name!
    % N.B. This does not reject the data from the CO2 check blocks as they
    % are already misisng from the cycle deltas and metadata.

iAirPlus = cycles.metadata.gasConfig == 'Air+';
cycles.deltas(~iAirPlus,:) = [];
cycles.metadata(~iAirPlus,:) = [];

% Could do something similar here with folder names e.g. only use data from
% the "SPICE" Folders:
% importedData = importedData(contains(importedData.FileHeader_Filename,'SPICE'),:);


% Second, reject cycles where the mass 28 beam is saturated on the cup
% (>9000 mV) or where it is less than 10 mV, presumably indicating some 
% problem with the source.
    % N.B. This also happens sometimes when measuring in an unusual gas
    % configuration that doesn't focus the 28 beam into a cup but collects 
    % the signal anyway. If the 28 beam isn't included in the gas config 
    % then its value is recorded as NaN.

iBeamReject = cycles.metadata.int28SA > 9000 | cycles.metadata.int28SA < 10 | cycles.metadata.int28ST > 9000 | cycles.metadata.int28ST < 10;
cycles.deltas(iBeamReject,:) = [];
cycles.metadata(iBeamReject,:) = [];

cycles.deltas.Properties.Description = join([cycles.deltas.Properties.Description ...
    "Cycles where the mass/charge 28 beam is either absent or saturated are omitted, as are samples measured in a gas configuration other than ""Air+""."]);
cycles.metadata.Properties.Description = join([cycles.metadata.Properties.Description ...
    "Cycles where the mass/charge 28 beam is either absent or saturated are omitted, as are samples measured in a gas configuration other than ""Air+""."]);

%% Reshape to Isotope Ratio x Delta Value x Block x Cycle
% Reshape the cycle deltas and their metadata into a more useable format.

varargout = {};
if includePIS && (includeAllBlocks || includeAllAliquots)
    [aliquots,aliquotsPis,aliquotsAll] = reshapeCycles(cycles,'includePIS',includePIS,'includeAllBlocks',includeAllBlocks,'includeAllAliquots',includeAllAliquots);
    varargout = [varargout {aliquotsPis aliquotsAll}];
elseif includePIS && ~(includeAllBlocks || includeAllAliquots)
    [aliquots,aliquotsPis] = reshapeCycles(cycles,'includePIS',includePIS,'includeAllBlocks',includeAllBlocks,'includeAllAliquots',includeAllAliquots);
    varargout = [varargout {aliquotsPis}];
elseif ~includePIS && (includeAllBlocks || includeAllAliquots)
    [aliquots,aliquotsAll] = reshapeCycles(cycles,'includePIS',includePIS,'includeAllBlocks',includeAllBlocks,'includeAllAliquots',includeAllAliquots);
    varargout = [varargout {aliquotsAll}];
else
    aliquots = reshapeCycles(cycles,'includePIS',includePIS,'includeAllBlocks',includeAllBlocks,'includeAllAliquots',includeAllAliquots);
end


%% Add the CO2 Check Data Back In
% Find the Timestamp for the CO2 Check Blocks
% I do this by interpolating to find the index/cycle number (y) of the
% block nearest to the CO2 blocks in datetime (x). I then use this index to
% find the datetime of these blocks.
% I have to use the index as the interpolant as the 'nearest' neighbour is
% calculated in the x variable and I want to interpolate to the nearest
% neightbour in datetime, not in cycle number.

iCO2 = importedData.Method=="CO2nonB" | importedData.Method=="CO2_non_B"; % Identify the CO2 blocks using the method name
rCO2N2_SA = importedData.x1_CycleInt_Samp_40(iCO2)./importedData.x1_CycleInt_Samp_29(iCO2);
rCO2N2_ST = importedData.x1_CycleInt_Ref_40(iCO2)./importedData.x1_CycleInt_Ref_29(iCO2);

x_temp = importedData.datenum + cumsum(ones(size(importedData.datenum))).*importedData.datenum*eps; % Generate an x-variable with unique points by adding a quasi-infinitesimal increment (eps) to each row
y_temp = 1:length(x_temp);

idxCO2NearestBlocks = interp1(x_temp(~iCO2),y_temp(~iCO2),x_temp(iCO2),'nearest'); % Interpolate to find the index of the nearest blocks
datetimeCO2NearestBlocks = importedData.datetime(idxCO2NearestBlocks); % Use the index to find the datetime of those blocks

[iA,idxB]=ismember(datetimeCO2NearestBlocks,aliquots.metadata.msDatetime); % Find the linear index of these blocks in the aliquot_metadata structure, not all will be present as some have <16 cycles or 6 < blocks < 4
[idxCO2Aliquots,~,~] = ind2sub(size(aliquots.metadata.msDatetime),idxB(iA)); % Convert to a subscript so that I know which aliquot they correspond to

aliquots.metadata.rCO2N2_SA = nan(size(aliquots.metadata.msDatetime,1),1);
aliquots.metadata.rCO2N2_ST = nan(size(aliquots.metadata.msDatetime,1),1);
aliquots.metadata.dCO2N2 = nan(size(aliquots.metadata.msDatetime,1),1);

aliquots.metadata.rCO2N2_SA(idxCO2Aliquots) = rCO2N2_SA(iA);
aliquots.metadata.rCO2N2_ST(idxCO2Aliquots) = rCO2N2_ST(iA);
aliquots.metadata.dCO2N2(idxCO2Aliquots) = (rCO2N2_SA(iA)./rCO2N2_ST(iA)-1)*1000;

aliquots.metadata.Properties.VariableDescriptions(end-2:end) = {
    'Sample CO2/N2 ratio for the aliquot', ...
    'Standard CO2/N2 ratio for the aliquot', ...
    'dCO2/N2 measured via manually peak jumping using the high voltage'};
aliquots.metadata.Properties.VariableUnits(end-2:end) = {'','',char(8240)};
    

%% Parse Outputs
% Organise the manipulated data into the requested outputs.

rawDataset = aliquots;


end