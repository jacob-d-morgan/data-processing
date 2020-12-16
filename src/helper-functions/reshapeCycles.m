function [aliquots,varargout] = reshapeCycles(cycles,varargin)
% RESHAPECYCLES converts cycle deltas to an aliquot-delta-block-cycle array
%   ALIQUOTS = RESHAPECYCLES(CYCLES) reshapes the data in CYCLES to the
%   structure ALIQUOTS, which contains most of the elements of CYCLES in
%   the structure fields 'metadata' and 'deltas', tables with variables
%   arranged as aliquot-by-block-by-cycle.
%   
%   The elements included in ALIQUOTS are the aliquots with the most
%   common number of regular (non-PIS) blocks, and the blocks with the most
%   common number of cycles. This greatly reduces the time the function
%   takes to run as the array can be constructed by indexing and reshaping,
%   rather than by filling a cell array in a for loop with many iterations.
%
%   The data ommitted from ALIQUOTS can be optionally outputted by
%   setting one or more of four name-value parameters to true:
%       - 'includePIS': If true, outputs the PIS blocks and their metadata 
%                       in a separate variable.
%       - 'includeAllBlocks': If true, outputs aliquots with the mode 
%                             number of blocks in a separate variable, 
%                             regardless of the number of cycles the blocks
%                             are made up of.
%       - 'includeAllAliquots': If true, outputs aliquots in a separate
%                               variable, regardless of the number of 
%                               blocks they are made up of.
%       - 'includeAllData': If true, outputs all aliquots in a separate
%                           varibale, regardless of the number of blocks or
%                           cycles they are made up of.
%
% The outputs that include all blocks or aliquots store data in cell arrays
% due to the different sizes of different blocks or aliquots.
%
% [...,ALIQUOTSPIS] = RESHAPECYCLES(...,'includePIS',true) 
% outputs the PIS experiment delta values and their metadata.
%
% [...,ALIQUOTSALLBLOCKS] = RESHAPECYCLES(...,'includeAllBlocks',true)
% outputs a cell array containing the aliquots composed of the mode number 
% of blocks of any number of cycles.
%
% [...,ALIQUOTSALLALIQUOTS] = RESHAPECYCLES(...,'includeAllAliquots',true)
% outputs a cell array containing aliquots of any number of blocks of the
% mode number of cycles.
%
% [...,ALLIQUOTSALLDATA] = RESHAPECYCLES(...,'includeAllData',true)
% outputs a cell array containing aliquots of any number of blocks of any
% number of cycles.
%
% The optional flags can be combined to output both the PIS data and the
% cell arrays of all the blocks/aliquots, e.g.:
%
% [...,ALIQUOTSpis,ALIQUOTSALLDATA] = 
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

addRequired(p,'cycles',@(x) validateattributes(x,{'struct'},{'scalar'}));
addParameter(p,'includePIS',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'includeAllBlocks',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'includeAllAliquots',false,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'includeAllData',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

% Parse Inputs and Assign Results
parse(p,cycles,varargin{:});

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
    [aliquots,aliquotsPis] = reshapeModesToAliquots(cycles,true);
    varargout = [varargout {aliquotsPis}];
else
    aliquots = reshapeModesToAliquots(cycles,false);
end

% Additionally, Reshape all the Blocks and/or Aliquots if Requested
% Reshapes all the blocks and/or aliquots, not just those with the most
% common number of cycles and blocks respectively. This is an optional
% output that can be requested by the function call and can be
% significantly slower than the default behaviour, depending on the size of
% the inputs.
if includeAllAliquots || includeAllBlocks
    aliquotsAll = reshapeAllToCell(cycles,includeAllBlocks,includeAllAliquots);
    if ~exist('aliquotsAll','var')
        warning('All blocks and aliquots were the same size. No cell array necessary.')
    end
    varargout = [varargout {aliquotsAll}];
end
end % end main

%% Local Functions

%  ------------------------------------------------------------------------

function [idxStartNewBlocks,idxStartNewAliquots] = detectBlocksAliquots(cycles)

% Identify the New Blocks
% Each block has a unique time stamp: find their indices.
[~,idxStartNewBlocks,~] = unique(cycles.metadata.msDatetime,'stable');

% Generate Block Metadata
blockMetadata = cycles.metadata(idxStartNewBlocks,:);

% Identify the New Aliquots
% Identify aliquots that were run with one of the three methods
SAMPLE_START_METHODS = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];
iStartMethods = false(length(blockMetadata.method(:,1)),1); iStartMethods(1) = true;
for ii=1:length(SAMPLE_START_METHODS)
    iStartMethods = iStartMethods | blockMetadata.method(:,1)==SAMPLE_START_METHODS(ii);
end

% Identify Changes in the User-Entered File Name
iNameChange = [1; ~strcmp(blockMetadata.fileNameUser(1:end-1),blockMetadata.fileNameUser(2:end))] & ...
              ~(contains(blockMetadata.fileNameUser,'PIS','IgnoreCase',true) | ...
                contains(blockMetadata.ID1,'PIS','IgnoreCase',true));

% Identify aliquots also by finding the start of a new sequence (value of sequence row decreases).
idxStartNewAliquots = find(iStartMethods | iNameChange | ([0; diff(blockMetadata.sequenceRow(:,1))]<0));

end

%  ------------------------------------------------------------------------

function [aliquots,varargout] = reshapeModesToAliquots(cycles,flagPIS)

% Identify the Blocks and Aliquots
[idxStartNewBlocks,idxStartNewAliquots] = detectBlocksAliquots(cycles);

% Calculate Block and Aliquot Lengths
blockLengthsAll = diff(idxStartNewBlocks);
blockLengthsAll(end+1) = height(cycles.metadata) - (idxStartNewBlocks(end)-1);

blockMetadata = cycles.metadata(idxStartNewBlocks,:);
aliquotLengthsAll = diff(idxStartNewAliquots);
aliquotLengthsAll(end+1) = height(blockMetadata)-(idxStartNewAliquots(end)-1);

% Identify Which Data To Use
% Find the aliquots that are composed of [four] blocks of [sixteen] cycles (or [five] blocks but have a PIS blocks)
MODE_BLOCK_LENGTH = mode(blockLengthsAll);
MODE_ALIQUOT_LENGTH = mode(aliquotLengthsAll);
iModeCycles = blockLengthsAll==MODE_BLOCK_LENGTH;
iModeBlocks = aliquotLengthsAll==MODE_ALIQUOT_LENGTH;
iPisAliquots = aliquotLengthsAll==MODE_ALIQUOT_LENGTH + 1 & (blockMetadata.ID1(idxStartNewAliquots+aliquotLengthsAll-1)=='PIS' | contains(blockMetadata.fileNameUser(idxStartNewAliquots+aliquotLengthsAll-1),'PIS','IgnoreCase',true)); % Is the last block in each aliquot have an ID1 equal to 'PIS'
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
iCyclesToUse = ismember(cycles.metadata.fileFullName,blockMetadata.fileFullName(idxBlocksToUse));

% Build the Array of Deltas and Metadata
% Simply reshape the relevant rows of cycle deltas into an array of the
% right dimensions. Make sure to fill the cycle, block, aliquot dimensions
% in that order and then permute after filling.
aliquot_deltas = table();
aliquot_deltas{:,:} = permute(reshape(cycles.deltas{iCyclesToUse,:},MODE_BLOCK_LENGTH,MODE_ALIQUOT_LENGTH,[],size(cycles.deltas,2)),[3 4 2 1]);
aliquot_deltas.Properties = cycles.deltas.Properties;
aliquot_deltas.Properties.Description = join([aliquot_deltas.Properties.Description, ...
    "Delta values are structured as aliquot-by-delta value-by-block-by-cycle. Only aliquots with the mode number of blocks of the mode number of cycles are included here."]);
    

aliquot_metadata = table();
for ii_var = string(cycles.metadata.Properties.VariableNames)
    aliquot_metadata.(ii_var) = permute(reshape(cycles.metadata.(ii_var)(iCyclesToUse,:),MODE_BLOCK_LENGTH,MODE_ALIQUOT_LENGTH,[],size(cycles.metadata.(ii_var),2)),[3 4 2 1]);
end
aliquot_metadata.Properties = cycles.metadata.Properties;
aliquot_metadata.Properties.Description = join([aliquot_metadata.Properties.Description, ...
    "Metadata are structured as aliquot-by-delta value-by-block-by-cycle. Only aliquots with the mode number of blocks of the mode number of cycles are included here."]);

aliquots = struct('metadata',aliquot_metadata,'deltas',aliquot_deltas);

varargout = {};
if flagPIS
    % Find the relevant PIS blocks
    iPisAliquotPisBlocks = false(length(blockLengthsAll),1);
    iPisAliquotPisBlocks(idxStartNewAliquots + MODE_ALIQUOT_LENGTH) = iAliquotsToUse & iPisAliquots; % Find where the PIS blocks are in the block_metadata array
    
    % Find the relevant PIS cycles
    iCyclesToUsePIS = ismember(cycles.metadata.fileFullName,blockMetadata.fileFullName(iPisAliquotPisBlocks)); % Find the cycles corresponding to the PIS blocks
    
    % Find where the PIS aliquots should go
    % Alternative approach is to make a NaN array and put the PIS blocks in
    % the position of the first block of the PIS aliquots and then reshape
    % as we do for the non-PIS blocks.
    iPisAliquotFirstBlocks = false(length(blockLengthsAll),1);
    iPisAliquotFirstBlocks(idxStartNewAliquots) = iAliquotsToUse & iPisAliquots;
    [~,idxPisAliquotsAfterReshape] = ismember(blockMetadata.fileFullName(iPisAliquotFirstBlocks),aliquot_metadata.fileFullName(:,1,1)); % Find where the PIS aliquots ended up after the reshape
    
    % Build the array of PIS deltas and metadata
    pis_aliquot_deltas = aliquots.deltas;
    pis_aliquot_deltas{:,:}(:) = 0;
    for ii_var = string(pis_aliquot_deltas.Properties.VariableNames)
        pis_aliquot_deltas.(ii_var)(:,:,2:end,:) = [];
    end
    pis_aliquot_deltas{idxPisAliquotsAfterReshape,:} = permute(reshape(cycles.deltas{iCyclesToUsePIS,:},MODE_BLOCK_LENGTH,1,[],size(cycles.deltas,2)),[3 4 2 1]);
    pis_aliquot_deltas = standardizeMissing(pis_aliquot_deltas,0);
    pis_aliquot_deltas.Properties = cycles.deltas.Properties;
    
    pis_aliquot_metadata = aliquot_metadata;
    for ii_var = string(cycles.metadata.Properties.VariableNames)
        pis_aliquot_metadata.(ii_var)(:,:) = missing;
        pis_aliquot_metadata.(ii_var)(idxPisAliquotsAfterReshape,:,1,:) = permute(reshape(cycles.metadata.(ii_var)(iCyclesToUsePIS,:),MODE_BLOCK_LENGTH,1,[],size(cycles.metadata.(ii_var),2)),[3 4 2 1]);
        pis_aliquot_metadata.(ii_var) = pis_aliquot_metadata.(ii_var)(:,:,1,:);
    end
    
    aliquotsPis = struct('metadata',pis_aliquot_metadata,'deltas',pis_aliquot_deltas);
    
    varargout = {aliquotsPis};
    
end
end % end reshapeModesToAliquot()

%  ------------------------------------------------------------------------

function aliquotsAll = reshapeAllToCell(cycles,iUseAllBlocks,iUseAllAliquots)

% Identify the Blocks and Aliquots
[idxStartNewBlocks,idxStartNewAliquots] = detectBlocksAliquots(cycles);

% Calculate Block and Aliquot Lengths
blockLengths = diff(idxStartNewBlocks);
blockLengths(end+1) = height(cycles.metadata) - (idxStartNewBlocks(end)-1);

blockMetadata = cycles.metadata(idxStartNewBlocks,:);
aliquotLengths = diff(idxStartNewAliquots);
aliquotLengths(end+1) = height(blockMetadata) - (idxStartNewAliquots(end)-1);

MODE_BLOCK_LENGTH = mode(blockLengths);
MODE_ALIQUOT_LENGTH = mode(aliquotLengths);
iModeCycles = blockLengths==MODE_BLOCK_LENGTH; % Blocks of [sixteen] cycles
iModeBlocks = aliquotLengths==MODE_ALIQUOT_LENGTH; % Aliquots of [four] blocks
iPisAliquots = (aliquotLengths==MODE_ALIQUOT_LENGTH+1) & (blockMetadata.ID1(idxStartNewAliquots+aliquotLengths-1)=='PIS'); % Aliquots where there is one extra block at the end that is labelled as a 'PIS' block
% Is this a valid way to ID PIS blocks? It will miss blocks that were maybe not explicitly tagged as PIS. Could see if looking at the actual P Imbalance improves this?

% Check if it's necessary to make a cell array or not
% If all the aliquots are four blocks of sixteen cycles then there is no
% need to re-do the reshaping we've already done
if all(blockLengths==MODE_BLOCK_LENGTH) && all(aliquotLengths==MODE_ALIQUOT_LENGTH)
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
aliquotsAll.metadata = table;
aliquotsAll.deltas = table;
for ii_var = string(cycles.metadata.Properties.VariableNames)
        aliquotsAll.metadata.(ii_var){NUM_ALIQUOTS,LONGEST_ALIQUOT} = [];
end
aliquotsAll.metadata.Properties = cycles.metadata.Properties;
aliquotsAll.metadata.Properties.Description = join([aliquotsAll.metadata.Properties.Description, ...
    "Metadata are structured as aliquot-by-delta value-by-block-by-cycle. Aliquots are included regardless of the number of blocks or the number of cycles."]);

for ii_var = string(cycles.deltas.Properties.VariableNames)
    aliquotsAll.deltas.(ii_var){NUM_ALIQUOTS,LONGEST_ALIQUOT} = [];
end
aliquotsAll.deltas.Properties = cycles.deltas.Properties;
aliquotsAll.deltas.Properties.Description = join([aliquotsAll.deltas.Properties.Description, ...
    "Delta values are structured as aliquot-by-delta value-by-block-by-cycle. Aliquots are included regardless of the number of blocks or the number of cycles."]);
    
%aliquot_deltas_all = cell(NUM_ALIQUOTS,LONGEST_ALIQUOT);
for ii = 1:length(idxStartNewBlocksToUse)
    if ismember(idxStartNewBlocksToUse(ii),idxStartNewBlocks(idxStartNewAliquots(iAliquotsToUse)))
        blockNum = 1;
        aliquotNum = aliquotNum + 1;
    else
        blockNum = blockNum + 1;
    end
    
    for ii_var = string(cycles.deltas.Properties.VariableNames)
        aliquotsAll.deltas.(ii_var){aliquotNum,blockNum} = cycles.deltas.(ii_var)(idxStartNewBlocksToUse(ii):idxStartNewBlocksToUse(ii) + blockLengthsToUse(ii)-1,:);
        %aliquot_deltas_all{aliquotNum,blockNum} = cycle_deltas(idxStartNewBlocksToUse(ii):idxStartNewBlocksToUse(ii) + blockLengthsToUse(ii)-1,:);
    end
    
    for ii_var = string(cycles.metadata.Properties.VariableNames)
        aliquotsAll.metadata.(ii_var){aliquotNum,blockNum} = cycles.metadata.(ii_var)(idxStartNewBlocksToUse(ii):idxStartNewBlocksToUse(ii) + blockLengthsToUse(ii)-1,:);
        %aliquot_metadata_all.(ii_var){aliquotNum,blockNum} = cycle_metadata.(fields)(idxStartNewBlocksToUse(ii):idxStartNewBlocksToUse(ii) + blockLengthsToUse(ii)-1,:);
    end
end
end % end reshapeAllToCell()
