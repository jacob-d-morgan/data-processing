function [deltaValues, metadata, varargout] = makeRawDataset(filesToImport, varargin)
% MAKERAWDATASET generates the delta values and metadata from FILES
%   [deltaValues, metadata] = makeRawDataset(filesToImport) imports the raw
%   cycle integrations from the .csv files in FILESTOIMPORT and performs
%   several manipulations to output an array of delta values and a
%   structure of arrays of metadata.
%
%   FILESTOIMPORT must be a string array of filenames or paths to files
%   that MATLAB can access.
%   DELTAVALUES is a 4 dimensional array of delta values arranged as
%   aliquot-by-delta value-by-block-by-cycle.
%   METADATA is a structure of 3 dimensional arrays arranged as
%   aliquot-by-block-by-cycle.
%
%   The manipulations performed to produce the outputs are as follows:
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
%    - [...,deltaValuesPIS, metadataPIS] = makeRawDataset(...,'includePIS')
%      includes delta values and their metadata from the PIS experiment
%      blocks in the same shape and structure as the primary outputs.
%    - [...,deltaValuesAllBlocks, metadataAllBlocks] = makeRawDataset(...,'includeAllBlocks')
%      outputs a cell array containing the aliquots composed of the mode
%      number of blocks of any number of cycles.
%    - [...,deltaValuesAllAliquots, metadataAllAliquots] = makeRawDataset(...,'includeAllAliquots')
%      outputs a cell array containing aliquots of any number of blocks of
%      the mode number of cycles.
%    - [...,deltaValuesAllAliquots, metadataAllAliquots] = makeRawDataset(...,'includeAllData')
%      outputs a cell array containing aliquots of any number of blocks of
%      any number of cycles.
%
% The optional flags can be combined to output both the PIS data and the
% variables containing all the blocks/aliquots:
%    - [...,PISDELTAS,PISMETADATA,ALLDATADELTAS,ALLDATAMETADATA] = 
%   makeRawDataset(...,'includePIS',includeAllData) outputs both the PIS
%   data and a cell array containing aliquots of any number of blocks of 
%   any number of cycles.
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
    nargoutchk(6,6);
elseif includePIS && ~(includeAllBlocks || includeAllAliquots)
    nargoutchk(4,4)
elseif ~includePIS && (includeAllBlocks || includeAllAliquots)
    nargoutchk(4,4)
else
    nargoutchk(2,2);
end

clear p


%% Import Data, Calculate Delta Values, and Compile Useful Metadata
% NOTE: Interpolating onto the ~isRef indices excludes the CO2 checks from
% the blocks I consider as they are recorded in isRef as 'true'. I have to
% manually add them back in below, including their metadata.

importedData = csvReadIsodat(filesToImport);

[cycle_deltas,cycle_metadata] = calcCycles(importedData,'IsRef__');

delta_names = string(cycle_deltas.Properties.VariableNames);
delta_labels = string(cycle_deltas.Properties.VariableDescriptions);
delta_units = string(cycle_deltas.Properties.VariableUnits);
cycle_deltas = table2array(cycle_deltas);

metadata_fields = string(fieldnames(cycle_metadata))'; % Transpose it to a row vector so that it works as a loop index


%% Reject Erroneous Cycles
% Get rid of some data that is of little use for evaluating the health of
% the mass spec and which causes problems with PIS etc.
%
% First, discard all data with a gas config different to 'Air+'
% This is mostly Jeff's neon experiments.
    % N.B. Note that there is a difference between Gas Config and Gas Name!
    % N.B. This does not reject the data from the CO2 check blocks as they
    % are already misisng from the cycle deltas and metadata.

iAirPlus = cycle_metadata.gasConfig == 'Air+';
cycle_deltas(~iAirPlus,:) = [];
for ii_field = metadata_fields
    cycle_metadata.(ii_field)(~iAirPlus,:) = [];
end

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

iBeamReject = cycle_metadata.int28SA > 9000 | cycle_metadata.int28SA < 10 | cycle_metadata.int28ST > 9000 | cycle_metadata.int28ST < 10;

cycle_deltas = cycle_deltas(~iBeamReject,:);
for ii_field = metadata_fields
    cycle_metadata.(ii_field) = cycle_metadata.(ii_field)(~iBeamReject);
end


%% Reshape to Isotope Ratio x Delta Value x Block x Cycle
% Reshape the cycle deltas and their metadata into a more useable format.

varargout = {};
if includePIS && (includeAllBlocks || includeAllAliquots)
    [aliquot_deltas,aliquot_metadata,aliquot_deltas_pis,aliquot_metadata_pis,aliquot_deltas_all,aliquot_metadata_all] = reshapeCycles(cycle_deltas,cycle_metadata,'includePIS',includePIS,'includeAllBlocks',includeAllBlocks,'includeAllAliquots',includeAllAliquots);
    varargout = [varargout {aliquot_deltas_pis aliquot_metadata_pis aliquot_deltas_all aliquot_metadata_all}];
elseif includePIS && ~(includeAllBlocks || includeAllAliquots)
    [aliquot_deltas,aliquot_metadata,aliquot_deltas_pis,aliquot_metadata_pis] = reshapeCycles(cycle_deltas,cycle_metadata,'includePIS',includePIS,'includeAllBlocks',includeAllBlocks,'includeAllAliquots',includeAllAliquots);
    varargout = [varargout {aliquot_deltas_pis aliquot_metadata_pis}];
elseif ~includePIS && (includeAllBlocks || includeAllAliquots)
    [aliquot_deltas,aliquot_metadata,aliquot_deltas_all,aliquot_metadata_all] = reshapeCycles(cycle_deltas,cycle_metadata,'includePIS',includePIS,'includeAllBlocks',includeAllBlocks,'includeAllAliquots',includeAllAliquots);
    varargout = [varargout {aliquot_deltas_all aliquot_metadata_all}];
else
    [aliquot_deltas,aliquot_metadata] = reshapeCycles(cycle_deltas,cycle_metadata,'includePIS',includePIS,'includeAllBlocks',includeAllBlocks,'includeAllAliquots',includeAllAliquots);
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

[iA,idxB]=ismember(datetimeCO2NearestBlocks,aliquot_metadata.msDatetime); % Find the linear index of these blocks in the aliquot_metadata structure, not all will be present as some have <16 cycles or 6 < blocks < 4
[idxCO2Aliquots,~,~] = ind2sub(size(aliquot_metadata.msDatetime),idxB(iA)); % Convert to a subscript so that I know which aliquot they correspond to

aliquot_metadata.rCO2N2_SA = nan(size(aliquot_metadata.msDatetime,1),1);
aliquot_metadata.rCO2N2_ST = nan(size(aliquot_metadata.msDatetime,1),1);
aliquot_metadata.dCO2N2 = nan(size(aliquot_metadata.msDatetime,1),1);

aliquot_metadata.rCO2N2_SA(idxCO2Aliquots) = rCO2N2_SA(iA);
aliquot_metadata.rCO2N2_ST(idxCO2Aliquots) = rCO2N2_ST(iA);
aliquot_metadata.dCO2N2(idxCO2Aliquots) = (rCO2N2_SA(iA)./rCO2N2_ST(iA)-1)*1000;


%% Parse Outputs
% Organise the manipulated data into the requested outputs.

deltaValues = aliquot_deltas;

metadata.delta_names = delta_names;
metadata.delta_labels = delta_labels;
metadata.delta_units = delta_units;
metadata.metadata = aliquot_metadata;


end