function [aliquotDeltas,aliquotMetadata,varargout] = reshapeCycles(cycle_deltas,cycle_metadata,varargin)
% RESHAPECYCLES converts cycle deltas to an aliquot-delta-block-cycle array
%   ALIQUOTDELTAS = RESHAPECYCLES(CYCLEDELTAS,CYCLEMETADATA) reshapes
%   CYCLEDELTAS and CYCLEMETADATA to the array ALIQUOTDELTAS and structure
%   of arrays ALIQUOTMETADATA, which contains most of the elements of
%   CYCLEDELTAS and CYCLEMETADATA in a aliquot-by-delta-by-block-by-cycle
%   structure.
%   The elements included in ALIQUOTDELTAS are the aliquots with the most
%   common number of regular (non-PIS) blocks, and the blocks with the most
%   common number of cycles. This greatly reduces the time the function
%   takes to run as the array can be constructed by indexing and reshaping,
%   rather than by filling a cell array in a for loop with many iterations.
%
%   The data ommitted from ALIQUOTDELTAS can be optionally outputted by
%   setting one or more of four name-value parameters to true:
%       - 'includePIS': If true, outputs the PIS blocks and their metadata 
%                       as separate variables.
%       - 'includeAllBlocks': If true, outputs two cell arrays containing
%                             all blocks and their metadata, regardless of 
%                             the number of cycles they are made up of.
%       - 'includeAllAliquots': If true, outputs two cell arrays containing
%                               all aliquots and their metadata, regardless
%                               of the number of blocks they are made up of.
%       - 'includeAllData': If true, outputs two cell arrays containing all
%                           aliquots and blocks, and their metadata, 
%                           regardless of the number of blocks or cycles 
%                           they are made up of.
%
% [...,PISDELTAS,PISMETADATA] = RESHAPECYCLES(...,'includePIS',true) 
% outputs the PIS experiment delta values and their metadata.
%
% [...,ALLBLOCKSDELTAS,ALLBLOCKSMETADATA] = RESHAPECYCLES(...,'includeAllBlocks',true)
% outputs a cell array containing the aliquots composed of the mode number 
% of blocks of any number of cycles.
%
% [...,ALLALIQUOTDELTAS,ALLALIQUOTMETADATA] = RESHAPECYCLES(...,'includeAllAliquots',true)
% outputs a cell array containing aliquots of any number of blocks of the
% mode number of cycles.
%
% [...,ALLDATADELTAS,ALLDATAMETADATA] = RESHAPECYCLES(...,'includeAllData',true)
% outputs a cell array containing aliquots of any number of blocks of any
% number of cycles.
%
% The optional flags can be combined to output both the PIS data and the
% cell arrays of all the blocks/aliquots, e.g.:
%
% [...,PISDELTAS,PISMETADATA,ALLDATADELTAS,ALLDATAMETADATA] = 
%   RESHAPECYCLES(...,'includePIS',true,'includeAllData',true) outputs both
%   the pis data and a cell array containing aliquots of any number of 
%   blocks of any number of cycles.


% -------------------------------------------------------------------------

%% Check Input/Outputs
% Check the user-provided inputs are of the correct form.

% Define Input Parser Scheme
p = inputParser;
p.StructExpand = false;
p.FunctionName = mfilename;

addRequired(p,'cycle_deltas',@(x) validateattributes(x,{'numeric'},{'2d','real'}));
addRequired(p,'cycle_metadata',@(x) validateattributes(x,{'struct'},{'scalar'}));
addParameter(p,'includePIS',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'includeAllBlocks',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'includeAllAliquots',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'includeAllData',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

% Parse Inputs and Assign Results
parse(p,cycle_deltas,cycle_metadata,varargin{:});

cycle_deltas = p.Results.cycle_deltas;
cycle_metadata = p.Results.cycle_metadata;
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


%% Calculate and Assign Outputs
varargout = {};

% Identify the Blocks and Aliquots from the Cycle Metadata
%[idxStartNewBlocks,idxStartNewAliquots] = detectBlocksAliquots(cycle_metadata);

% Fast Reshape Most of the Data
% Reshapes only the blocks and aliquots with the most common number of
% cycles and blocks respectively. This is the default behaviour of
% reshapeCycles and is faster than including all of data.
% Also, include the PIS blocks if requested by the function call.
if includePIS
    [aliquotDeltas,aliquotMetadata,pis_aliquot_deltas,pis_aliquot_metadata] = reshapeModesToAliquots(cycle_deltas,cycle_metadata,includePIS);
    varargout = [varargout {pis_aliquot_deltas pis_aliquot_metadata}];
else
    [aliquotDeltas,aliquotMetadata] = reshapeModesToAliquots(cycle_deltas,cycle_metadata,false);
end

% Additionally, Reshape all the Blocks and/or Aliquots if Requested
% Reshapes all the blocks and/or aliquots, not just those with the most
% common number of cycles and blocks respectively. This is an optional
% output that can be requested by the function call and can be
% significantly slower than the default behaviour, depending on the size of
% the inputs.
if includeAllAliquots || includeAllBlocks
    [aliquotDeltasAll,aliquotMetadataAll] = reshapeAllToCell(cycle_deltas,cycle_metadata,includeAllBlocks,includeAllAliquots);
    if ~exist('aliquotDeltasAll','var')
        warning('All blocks and aliquots were the same size. No cell array necessary.')
    end
    varargout = [varargout {aliquotDeltasAll aliquotMetadataAll}];
end
end % end main

%% Local Functions

%  ------------------------------------------------------------------------

function [idxStartNewBlocks,idxStartNewAliquots,block_metadata] = detectBlocksAliquots(cycle_metadata)

% Identify the New Blocks
% Each block has a unique filename: find their indices.
[~,idxStartNewBlocks,~] = unique(cycle_metadata.filename,'stable');

% Generate Block Metadata
for fields = string(fieldnames(cycle_metadata))'
    block_metadata.(fields) = cycle_metadata.(fields)(idxStartNewBlocks,:);
end

% Identify the New Aliquots
% Define Methods that Are Run at the Start of Each New Sample
SAMPLE_START_METHODS = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];

% Identify aliquots that were run with one of the three methods
iWorking = false(length(block_metadata.method(:,1)),1); iWorking(1) = true;
for ii=1:length(SAMPLE_START_METHODS)
    iWorking = iWorking | block_metadata.method(:,1)==SAMPLE_START_METHODS(ii);
end

% Identify aliquots also by finding the start of a new sequence (value of sequence row decreases).
idxStartNewAliquots = find(iWorking | ([0; diff(block_metadata.sequenceRow(:,1))]<0));

end

%  ------------------------------------------------------------------------

function [aliquot_deltas,aliquot_metadata,varargout] = reshapeModesToAliquots(cycle_deltas,cycle_metadata,flagPIS)

% Identify the Blocks and Aliquots
[idxStartNewBlocks,idxStartNewAliquots,block_metadata] = detectBlocksAliquots(cycle_metadata);

% Calculate Block and Aliquot Lengths
blockLengthsAll = diff(idxStartNewBlocks);
blockLengthsAll(end+1) = length(cycle_metadata.filename) - (idxStartNewBlocks(end)-1);

aliquotLengthsAll = diff(idxStartNewAliquots);
aliquotLengthsAll(end+1)=length(block_metadata.filename)-(idxStartNewAliquots(end)-1);

% Identify Which Data To Use
% Find the aliquots that are composed of [four] blocks of [sixteen] cycles (or [five] blocks but have a PIS blocks)
MODE_BLOCK_LENGTH = mode(blockLengthsAll);
MODE_ALIQUOT_LENGTH = mode(aliquotLengthsAll);
iModeCycles = blockLengthsAll==MODE_BLOCK_LENGTH;
iModeBlocks = aliquotLengthsAll==MODE_ALIQUOT_LENGTH;
iPisAliquots = aliquotLengthsAll==MODE_ALIQUOT_LENGTH + 1 & block_metadata.ID1(idxStartNewAliquots+aliquotLengthsAll-1)=='PIS'; % Is the last block in each aliquot have an ID1 equal to 'PIS'
% Is this a valid way to ID PIS blocks? It will miss blocks that were maybe not explicitly tagged as PIS. Could see if looking at the actual P Imbalance improves this?

% Find the relevant aliquots
iAliquotsToUse = false(size(idxStartNewAliquots));
for ii = 1:length(idxStartNewAliquots)
    iAliquotsToUse(ii) = all(iModeCycles(idxStartNewAliquots(ii):idxStartNewAliquots(ii)+aliquotLengthsAll(ii)-1)) & (iModeBlocks(ii) | iPisAliquots(ii));
end

% Find the relevant blocks
iBlocksToUse = false(length(blockLengthsAll),1);
iBlocksToUse(idxStartNewAliquots) = iAliquotsToUse;
idxBlocksToUse = find(iBlocksToUse) + (0:MODE_ALIQUOT_LENGTH-1); % Include only the four measurement blocks, not the PIS block at the end

% Find the relevant cycles
iCyclesToUse = ismember(cycle_metadata.filename,block_metadata.filename(idxBlocksToUse)); 

% Build the Array of Deltas and Metadata
% Simply reshape the relevant rows of cycle deltas into an array of the
% right dimensions. Make sure to fill the cycle, block, aliquot dimensions
% in that order and then permute after filling.
aliquot_deltas = reshape(cycle_deltas(iCyclesToUse,:),MODE_BLOCK_LENGTH,MODE_ALIQUOT_LENGTH,[],size(cycle_deltas,2));
aliquot_deltas = permute(aliquot_deltas,[3 4 2 1]);
for fields = string(fieldnames(cycle_metadata))'
    aliquot_metadata.(fields) = reshape(cycle_metadata.(fields)(iCyclesToUse,:),MODE_BLOCK_LENGTH,MODE_ALIQUOT_LENGTH,[]);
    aliquot_metadata.(fields) = permute(aliquot_metadata.(fields),[3 2 1]);
end

varargout = {};
if flagPIS
    % Find the relevant PIS blocks
    iPisAliquotPisBlocks = false(length(blockLengthsAll),1);
    iPisAliquotPisBlocks(idxStartNewAliquots + MODE_ALIQUOT_LENGTH) = iAliquotsToUse & iPisAliquots; % Find where the PIS blocks are in the block_metadata array
    
    % Find the relevant PIS cycles
    iCyclesToUsePIS = ismember(cycle_metadata.filename,block_metadata.filename(iPisAliquotPisBlocks)); % Find the cycles corresponding to the PIS blocks
    
    % Find where the PIS aliquots should go
    % Alternative approach is to make a NaN array and put the PIS blocks in
    % the position of the first block of the PIS aliquots and then reshape
    % as we do for the non-PIS blocks.
    iPisAliquotFirstBlocks = false(length(blockLengthsAll),1);
    iPisAliquotFirstBlocks(idxStartNewAliquots) = iAliquotsToUse & iPisAliquots;
    [~,idxPisAliquotsAfterReshape] = ismember(block_metadata.filename(iPisAliquotFirstBlocks),aliquot_metadata.filename(:,1,1)); % Find where the PIS aliquots ended up after the reshape
    
    % Build the array of PIS deltas and metadata
    pis_aliquot_deltas = nan(size(aliquot_deltas));
    pis_aliquot_deltas(idxPisAliquotsAfterReshape,:,1,:) = permute(reshape(cycle_deltas(iCyclesToUsePIS,:),MODE_BLOCK_LENGTH,1,[],size(cycle_deltas,2)),[3 4 2 1]);
    pis_aliquot_deltas = pis_aliquot_deltas(:,:,1,:); % Only one PIS block per aliquot, discard the rest

    pis_aliquot_metadata = aliquot_metadata;
    for fields = string(fieldnames(cycle_metadata))'
        pis_aliquot_metadata.(fields)(:,:) = missing;
        pis_aliquot_metadata.(fields)(idxPisAliquotsAfterReshape,1,:) = permute(reshape(cycle_metadata.(fields)(iCyclesToUsePIS,:),MODE_BLOCK_LENGTH,1,[]),[3 2 1]);
        pis_aliquot_metadata.(fields) = pis_aliquot_metadata.(fields)(:,1,:);
    end
    
    varargout = {pis_aliquot_deltas pis_aliquot_metadata};
end
end % end reshapeModesToArr()

%  ------------------------------------------------------------------------

function [aliquot_deltas_all,aliquot_metadata_all] = reshapeAllToCell(cycle_deltas,cycle_metadata,iUseAllBlocks,iUseAllAliquots)

% Identify the Blocks and Aliquots
[idxStartNewBlocks,idxStartNewAliquots,block_metadata] = detectBlocksAliquots(cycle_metadata);

% Calculate Block and Aliquot Lengths
blockLengths = diff(idxStartNewBlocks);
blockLengths(end+1) = length(cycle_metadata.filename) - (idxStartNewBlocks(end)-1);

aliquotLengths = diff(idxStartNewAliquots);
aliquotLengths(end+1)=length(block_metadata.method)-(idxStartNewAliquots(end)-1);

MODE_BLOCK_LENGTH = mode(blockLengths);
MODE_ALIQUOT_LENGTH = mode(aliquotLengths);
iModeCycles = blockLengths==MODE_BLOCK_LENGTH; % Blocks of [sixteen] cycles
iModeBlocks = aliquotLengths==MODE_ALIQUOT_LENGTH; % Aliquots of [four] blocks
iPisAliquots = (aliquotLengths==MODE_ALIQUOT_LENGTH+1) & (block_metadata.ID1(idxStartNewAliquots+aliquotLengths-1)=='PIS'); % Aliquots where there is one extra block at the end that is labelled as a 'PIS' block
% Is this a valid way to ID PIS blocks? It will miss blocks that were maybe not explicitly tagged as PIS. Could see if looking at the actual P Imbalance improves this?

% Check if it's necessary to make a cell array or not
% If all the aliquots are four blocks of sixteen cycles then there is no
% need to re-do the reshaping we've already done
if all(blockLengths==MODE_BLOCK_LENGTH) || all(aliquotLengths==MODE_ALIQUOT_LENGTH)
    return;
end

% Identify Which Data to Use
% Find the aliquots that are composed of the right number of blocks of the
% right number of cycles (or one too many blocks but have a PIS block)
if iUseAllBlocks && iUseAllAliquots % Use BLOCKS of ALL LENGTHS and ALIQUOTS of ALL LENGTHS
    iAliquotsToUse = true(length(aliquotLengths),1);
    iWorking = true(length(blockLengths),1);
    idxBlocksToUse = find(iWorking);
    
elseif iUseAllBlocks && ~iUseAllAliquots % Use BLOCKS of ALL LENGTHS but ONLY ALIQUOTS of [four] blocks
    % Identify Which Aliquots To Use
    iAliquotsToUse = iModeBlocks | iPisAliquots;
    
    % Find the relevant blocks
    iWorking = false(length(blockLengths),1);
    iWorking(idxStartNewAliquots) = iAliquotsToUse;
    idxBlocksToUse = find(iWorking) + (0:MODE_ALIQUOT_LENGTH-1); % Include only the four measurement blocks, not the PIS block at the end
    idxBlocksToUse = sort(idxBlocksToUse(:));
    
elseif ~iUseAllBlocks && iUseAllAliquots % Use ONLY BLOCKS of [sixteen cycles] but ALIQUOTS of ALL LENGTHS
    % Identify Which Aliquots To Use
    iAliquotsToUse = false(size(idxStartNewAliquots));
    for ii = 1:length(idxStartNewAliquots)
        iAliquotsToUse(ii) = all(iModeCycles(idxStartNewAliquots(ii):idxStartNewAliquots(ii)+aliquotLengths(ii)-1));
    end
    
    % Find the relevant blocks
    iWorking = false(length(blockLengths),1);
    for ii = 1:length(idxStartNewAliquots)
        iWorking(idxStartNewAliquots(ii):idxStartNewAliquots(ii) + aliquotLengths(ii)-1) = iAliquotsToUse(ii);
    end
    idxBlocksToUse = find(iWorking);
    
end

% Build the Arrays of Deltas and Metadata
NUM_ALIQUOTS = sum(iAliquotsToUse);
LONGEST_ALIQUOT = max(aliquotLengths(iAliquotsToUse));
if iUseAllBlocks && ~iUseAllAliquots
    LONGEST_ALIQUOT = LONGEST_ALIQUOT - 1; % Because we won't inlude any PIS blocks
end

idxStartNewBlocksToUse = idxStartNewBlocks(idxBlocksToUse);
blockLengthsToUse = blockLengths(idxBlocksToUse);

aliquotNum = 0;
blockNum = 0;
aliquot_deltas_all = cell(NUM_ALIQUOTS,LONGEST_ALIQUOT);
for ii = 1:length(idxStartNewBlocksToUse)
    if ismember(idxStartNewBlocksToUse(ii),idxStartNewBlocks(idxStartNewAliquots(iAliquotsToUse)))
        blockNum = 1;
        aliquotNum = aliquotNum + 1;
    else
        blockNum = blockNum + 1;
    end
    
    aliquot_deltas_all{aliquotNum,blockNum} = cycle_deltas(idxStartNewBlocksToUse(ii):idxStartNewBlocksToUse(ii) + blockLengthsToUse(ii)-1,:);
    for fields = string(fieldnames(cycle_metadata))'
        aliquot_metadata_all.(fields){aliquotNum,blockNum} = cycle_metadata.(fields)(idxStartNewBlocksToUse(ii):idxStartNewBlocksToUse(ii) + blockLengthsToUse(ii)-1,:);
    end
end
end % end reshapeAllToCell()
