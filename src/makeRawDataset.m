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






end