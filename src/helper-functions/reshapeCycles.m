function [aliquotDeltas,aliquotMetadata,varargout] = reshapeCycles(cycle_deltas,cycle_metadata,varargin)
% RESHAPECYCLES converts cycle deltas to an aliquot-delta-block-cycle array
%   ALIQUOTDELTAS = RESHAPECYCLES(CYCLEDELTAS) reshapes CYCLEDELTAS to the
%   array ALIQUOTDELTAS, which contains most of the elements of CYCLEDELTAS
%   in a aliquot-by-delta-by-block-by-cycle structure.
%   The elements included in ALIQUOTDELTAS are the aliquots with the most
%   common number of regular (non-PIS) blocks, and the blocks with the most
%   common number of cycles. This greatly reduces the time the function
%   takes to run as the array can be constructed by indexing and reshaping,
%   rather than by filling a cell array in a for loop with many iterations.
%
%   The data ommitted from ALIQUOTDELTAS can be optionally outputted by 
%   including one or more of four optional flags (case-insensitive):
%       - 'includePIS': Outputs the PIS blocks and their metadata as 
%                       separate variables.
%       - 'includeAllBlocks': Outputs two cell arrays containing all blocks
%                             and their metadata, regardless of the number
%                             of cycles they are made up of.
%       - 'includeAllAliquots': Outputs two cell arrays containing all 
%                               aliquots and their metadata, regardless of
%                               the number of blocks they are made up of.


% -------------------------------------------------------------------------

%% Parse Inputs
nargoutchk(2,6);

flagPIS = false;
flagUseAllBlocks = false;
flagUseAllAliquots = false;

% Ouptut PIS Blocks?
if any(contains(varargin,'includePIS','IgnoreCase',true))
    flagPIS = true;
    nargoutchk(4,6);
end

% Include All Blocks?
if any(contains(varargin,'includeAllBlocks','IgnoreCase',true) | contains(varargin,'includeAllData','IgnoreCase',true))
    flagUseAllBlocks = true;
    nargoutchk(4,6)
end

% Include All Aliquots?
if any(contains(varargin,'includeAllAliquots','IgnoreCase',true) | contains(varargin,'includeAllData','IgnoreCase',true))
    flagUseAllAliquots = true;
    nargoutchk(4,6)
end

if ~flagPIS && ~flagUseAllBlocks && ~flagUseAllAliquots
    nargoutchk(2,2);
end

%% Calculate and Assign Outputs
varargout = {};

% Fast Reshape Most of the Data
% Reshapes only the blocks and aliquots with the most common number of
% cycles and blocks respectively. This is the default behaviour of
% reshapeCycles and is faster than including all of data.
% Also, include the PIS blocks if requested by the function call.
[blockDeltas,blockMetadata] = reshapeToBlocks(cycle_deltas,cycle_metadata,false);

if flagPIS
    [aliquotDeltas,aliquotMetadata,pisAliquotDeltas,pisAliquotMetadata] = reshapeArrToAliquots(blockDeltas,blockMetadata,false,true);
    varargout = [varargout {pisAliquotDeltas pisAliquotMetadata}];
else
    [aliquotDeltas,aliquotMetadata] = reshapeArrToAliquots(blockDeltas,blockMetadata,false,false);
end

% Additionally, Reshape all the Blocks and/or Aliquots if Requested
% Reshapes all the blocks and/or aliquots, not just those with the most
% common number of cycles and blocks respectively. This is an optional
% output that can be requested by the function call and can be
% significantly slower than the default behaviour, depending on the size of
% the inputs.
if ~flagUseAllBlocks && flagUseAllAliquots
    [aliquotDeltasAll,aliquotMetadataAll] = reshapeArrToAliquots(blockDeltas,blockMetadata,true,false);
    varargout = [varargout {aliquotDeltasAll aliquotMetadataAll}];
    
elseif flagUseAllBlocks && ~flagUseAllAliquots
    [blockDeltasAll,blockMetadataAll] = reshapeToBlocks(cycle_deltas,cycle_metadata,true);
    [aliquotDeltasAll,aliquotMetadataAll] = reshapeCellToAliquots(blockDeltasAll,blockMetadataAll,false);
    varargout = [varargout {aliquotDeltasAll aliquotMetadataAll}];
    
elseif flagUseAllBlocks && flagUseAllAliquots
    [blockDeltasAll,blockMetadataAll] = reshapeToBlocks(cycle_deltas,cycle_metadata,true);
    [aliquotDeltasAll,aliquotMetadataAll] = reshapeCellToAliquots(blockDeltasAll,blockMetadataAll,true);
    varargout = [varargout {aliquotDeltasAll aliquotMetadataAll}];
end


%% ------------------------------------------------------------------------
% Nested Functions
    
    function [block_deltas,block_metadata] = reshapeToBlocks(cycle_deltas,cycle_metadata,flagAllBlocks)
        % Reshapes the cycle_deltas and cycle_metadata arrays to the
        % block_deltas and block_metadata outputs. The output is a cell
        % array of all blocks if the flag is set to true and is an array of
        % only the blocks with the most common number of cycles if the flag
        % is set to false.
        %
        % -----------------------------------------------------------------
        
        % Identify the New Blocks
        % Each block has a unique filename: find their indices.
        [~,idxStartNewBlockAll,~] = unique(cycle_metadata.filename,'stable');
        blockLengthsAll = diff(idxStartNewBlockAll);
        blockLengthsAll(end+1) = length(cycle_deltas) - (idxStartNewBlockAll(end)-1);
        
        if flagAllBlocks % Use all the blocks
            NUM_ALL_BLOCKS = length(idxStartNewBlockAll);
            
            block_deltas = cell(NUM_ALL_BLOCKS,1);
            for ii_BLOCK_TO_FILL = 1:NUM_ALL_BLOCKS
                block_deltas{ii_BLOCK_TO_FILL} = cycle_deltas(idxStartNewBlockAll(ii_BLOCK_TO_FILL):idxStartNewBlockAll(ii_BLOCK_TO_FILL)+blockLengthsAll(ii_BLOCK_TO_FILL)-1,:);
                
                for fields = string(fieldnames(cycle_metadata))'
                    block_metadata.(fields){ii_BLOCK_TO_FILL,1} = cycle_metadata.(fields)(idxStartNewBlockAll(ii_BLOCK_TO_FILL):idxStartNewBlockAll(ii_BLOCK_TO_FILL)+blockLengthsAll(ii_BLOCK_TO_FILL)-1);
                end
            end
            
        else % Take only the blocks with the most common number of cycles!
            MODE_BLOCK_LENGTH = mode(blockLengthsAll);
            idxStartNewBlockUse = idxStartNewBlockAll(blockLengthsAll==MODE_BLOCK_LENGTH);
            NUM_BLOCKS = length(idxStartNewBlockUse);
            
            iCyclesToUse = ismember(cycle_metadata.filename,cycle_metadata.filename(idxStartNewBlockUse));
            
            % Fill the Array of Delta Values
            block_deltas = reshape(cycle_deltas(iCyclesToUse,:),MODE_BLOCK_LENGTH,NUM_BLOCKS,[]);
            
            for fields = string(fieldnames(cycle_metadata))'
                block_metadata.(fields) = reshape(cycle_metadata.(fields)(iCyclesToUse,:),MODE_BLOCK_LENGTH,NUM_BLOCKS);
            end
        end
    end

    function [aliquot_deltas,aliquot_metadata,pis_aliquot_deltas,pis_aliquot_metadata] = reshapeArrToAliquots(block_deltas,block_metadata,flagAllAliquots,flagPIS)
        % Reshapes the block_deltas and block_metadata inputs to the
        % aliquot_deltas and aliquot_metadata outputs. 
        % If the input is a cell array (or structure of cell arrays in the
        % case of the metadata), the output is also a cell array of all
        % aliquots if the flag is set to true or is a cell array of only
        % the aliquots with the most common number of blocks if the flag is
        % set to false.
        % If the input is a numeric array (or structure of numeric arrays
        % in the case of the metadata), the output is a cell array of
        % all the aliquots if the flag is set to true, or is a numberic
        % array of only the aliquots with the most common number of blocks
        % if the flag is set to true.
        %
        % -----------------------------------------------------------------
        
        % Check Number of Output Arguments
        if flagPIS
            nargoutchk(4,4)
        else
            nargoutchk(2,2)
        end
        
        % Identify the New Aliquots 
        % Define Methods that Are Run at the Start of Each New Sample
        SAMPLE_START_METHODS = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];
               
        % Identify aliquots that were run with one of the three methods
        iWorking = false(1,size(block_metadata.method(1,:),2));
        for ii=1:length(SAMPLE_START_METHODS)
            iWorking = iWorking | block_metadata.method(1,:)==SAMPLE_START_METHODS(ii);
        end
        
        % Identify aliquots also by finding the start of a new sequence (value of sequence row decreases).
        idxStartNewAliquotAll = find(iWorking | ([0 diff(block_metadata.sequenceRow(1,:))]<0));
        aliquotLengthsAll = diff(idxStartNewAliquotAll);
        aliquotLengthsAll(end+1)=length(block_deltas)-(idxStartNewAliquotAll(end)-1);
        
        if flagAllAliquots % Use all aliquots
            NUM_ALL_ALIQUOTS = length(idxStartNewAliquotAll);
            
            % Fill the Cell Array of Delta Values
            aliquot_deltas = cell(NUM_ALL_ALIQUOTS,1);
            for ii_ALIQUOT_TO_FILL = 1:NUM_ALL_ALIQUOTS
                aliquot_deltas{ii_ALIQUOT_TO_FILL} = block_deltas(:,idxStartNewAliquotAll(ii_ALIQUOT_TO_FILL):idxStartNewAliquotAll(ii_ALIQUOT_TO_FILL)+aliquotLengthsAll(ii_ALIQUOT_TO_FILL)-1,:);
                aliquot_deltas{ii_ALIQUOT_TO_FILL} = permute(aliquot_deltas{ii_ALIQUOT_TO_FILL},[3 2 1]);
                
                for fields = string(fieldnames(block_metadata))'
                    aliquot_metadata.(fields){ii_ALIQUOT_TO_FILL,1} = squeeze(block_metadata.(fields)(:,idxStartNewAliquotAll(ii_ALIQUOT_TO_FILL):idxStartNewAliquotAll(ii_ALIQUOT_TO_FILL) + aliquotLengthsAll(ii_ALIQUOT_TO_FILL)-1));
                    aliquot_metadata.(fields){ii_ALIQUOT_TO_FILL,1} = permute(aliquot_metadata.(fields){ii_ALIQUOT_TO_FILL,1},[2 1]);
                end
            end
            
        else % Take only the aliquots with the most common numner of blocks! 
        % NB: Also include aliquots with 1 block more than the mode, where
        % this extra block is a PIS block. We still want to include the
        % main measurement blocks in the aliquot_delta output, just not the
        % PIS block. The PIS blocks are not included in the aliquot_deltas
        % output (they are filtered into a different variable instead).
            
            iPIS = block_metadata.ID1(1,idxStartNewAliquotAll+aliquotLengthsAll-1)=='PIS'; % Is the last block in each aliquot have an ID1 equal to 'PIS'
            % Is this a valid way to ID PIS blocks? It will miss blocks that were maybe not explicitly tagged as PIS. Could see if looking at the actual P Imbalance improves this?
            
            MODE_ALIQUOT_LENGTH = mode(aliquotLengthsAll);
            idxStartNewAliquotUse = idxStartNewAliquotAll(aliquotLengthsAll==MODE_ALIQUOT_LENGTH | ((aliquotLengthsAll==MODE_ALIQUOT_LENGTH+1) & iPIS));
            NUM_ALIQUOTS = length(idxStartNewAliquotUse);
            
            % Identify the Blocks to Inlcude
            idxBlocksToUse = idxStartNewAliquotUse + (0:MODE_ALIQUOT_LENGTH-1)'; % PIS blocks are excluded here: + 0:MODE_ALIQUOT_LENGTH-1 excludes the final block
            %block_deltas = block_deltas(:,idxBlocksToUse,:);
            
            % Fill the Array of Delta Values
            aliquot_deltas = reshape(block_deltas(:,idxBlocksToUse,:),size(block_deltas,1),MODE_ALIQUOT_LENGTH,NUM_ALIQUOTS,size(block_deltas,3));
            aliquot_deltas = permute(aliquot_deltas,[3 4 2 1]);
            
            
            for fields = string(fieldnames(block_metadata))'
                aliquot_metadata.(fields) = reshape(block_metadata.(fields)(:,idxBlocksToUse,:),size(block_deltas,1),MODE_ALIQUOT_LENGTH,NUM_ALIQUOTS);
                aliquot_metadata.(fields) = permute(aliquot_metadata.(fields),[3 2 1]);
            end
            
            if flagPIS
                idxPisAliquots = idxStartNewAliquotAll(aliquotLengthsAll==MODE_ALIQUOT_LENGTH+1 & iPIS);
                idxPisBlocks = idxPisAliquots + aliquotLengthsAll(aliquotLengthsAll==MODE_ALIQUOT_LENGTH+1 & iPIS)-1;
                pis_block_deltas = nan(size(block_deltas));
                pis_block_deltas(:,idxPisAliquots,:) = block_deltas(:,idxPisBlocks,:);
                                
                pis_aliquot_deltas = reshape(pis_block_deltas(:,idxBlocksToUse,:),size(block_deltas,1),MODE_ALIQUOT_LENGTH,NUM_ALIQUOTS,size(block_deltas,3));
                pis_aliquot_deltas = permute(pis_aliquot_deltas,[3 4 2 1]);
                pis_aliquot_deltas = pis_aliquot_deltas(:,:,1,:);
                
                pis_block_metadata = block_metadata;
                for fields = string(fieldnames(block_metadata))'
                    pis_block_metadata.(fields)(:,:) = missing;
                    pis_block_metadata.(fields)(:,idxPisAliquots) = block_metadata.(fields)(:,idxPisBlocks);
                    pis_aliquot_metadata.(fields) = permute(reshape(pis_block_metadata.(fields)(:,idxBlocksToUse,:),size(block_deltas,1),MODE_ALIQUOT_LENGTH,NUM_ALIQUOTS),[3 2 1]);
                    pis_aliquot_metadata.(fields) = pis_aliquot_metadata.(fields)(:,1,:);
                end
            end
        end
    end  
        
    function [aliquot_deltas,aliquot_metadata] = reshapeCellToAliquots(block_deltas,block_metadata,flagAllAliquots)
        % Reshapes the block_deltas and block_metadata inputs to the
        % aliquot_deltas and aliquot_metadata outputs.
        % If the input is a cell array (or structure of cell arrays in the
        % case of the metadata), the output is also a cell array of all
        % aliquots if the flag is set to true or is a cell array of only
        % the aliquots with the most common number of blocks if the flag is
        % set to false.
        % If the input is a numeric array (or structure of numeric arrays
        % in the case of the metadata), the output is a cell array of
        % all the aliquots if the flag is set to true, or is a numberic
        % array of only the aliquots with the most common number of blocks
        % if the flag is set to true.
        %
        % -----------------------------------------------------------------
        
        % Identify the New Aliquots
        % Define Methods that Are Run at the Start of Each New Sample
        SAMPLE_START_METHODS = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];
        
        blockMethods = cellfun(@(x) x(1),block_metadata.method)';
        blockSequenceRow = cellfun(@(x) x(1),block_metadata.sequenceRow)';
        
        % Identify aliquots that were run with one of the three methods
        iWorking = false(1,size(blockMethods,2));
        for ii=1:length(SAMPLE_START_METHODS)
            iWorking = iWorking | blockMethods==SAMPLE_START_METHODS(ii);
        end
        
        % Identify aliquots also by finding the start of a new sequence (value of sequence row decreases).
        idxStartNewAliquotAll = find(iWorking | ([0 diff(blockSequenceRow)]<0));
        aliquotLengthsAll = diff(idxStartNewAliquotAll);
        aliquotLengthsAll(end+1)=length(block_deltas)-(idxStartNewAliquotAll(end)-1);
        
        
        if flagAllAliquots % Use All Aliquots
            NUM_ALL_ALIQUOTS = length(idxStartNewAliquotAll);
            
            aliquot_deltas = cell(NUM_ALL_ALIQUOTS,1);
            for ii_ALIQUOT_TO_FILL = 1:NUM_ALL_ALIQUOTS
                aliquot_deltas{ii_ALIQUOT_TO_FILL} = block_deltas(idxStartNewAliquotAll(ii_ALIQUOT_TO_FILL):idxStartNewAliquotAll(ii_ALIQUOT_TO_FILL) + aliquotLengthsAll(ii_ALIQUOT_TO_FILL)-1);
                for fields = string(fieldnames(block_metadata))'
                    aliquot_metadata.(fields){ii_ALIQUOT_TO_FILL,1} = block_metadata.(fields)(idxStartNewAliquotAll(ii_ALIQUOT_TO_FILL):idxStartNewAliquotAll(ii_ALIQUOT_TO_FILL) + aliquotLengthsAll(ii_ALIQUOT_TO_FILL)-1);
                end
            end
            
        else % Take only the aliquots with the most common numner of blocks!
        % NB: Also include aliquots with 1 block more than the mode, where
        % this extra block is a PIS block. We still want to include the
        % main measurement blocks in the aliquot_delta output, just not the
        % PIS block. The PIS blocks are not included in the aliquot_deltas
        % output (they are filtered out below).
            
            blockID1 = cellfun(@(x) x(1), block_metadata.ID1)';
            iPIS = blockID1(idxStartNewAliquotAll+aliquotLengthsAll-1)=='PIS'; % Is the last block in each aliquot have an ID1 equal to 'PIS'
            % Is this a valid way to ID PIS blocks? It will miss blocks that were maybe not explicitly tagged as PIS. Could see if looking at the actual P Imbalance improves this?
            
            MODE_ALIQUOT_LENGTH = mode(aliquotLengthsAll);
            idxStartNewAliquotUse = idxStartNewAliquotAll(aliquotLengthsAll==MODE_ALIQUOT_LENGTH | ((aliquotLengthsAll==MODE_ALIQUOT_LENGTH+1) & iPIS));
            NUM_ALIQUOTS = length(idxStartNewAliquotUse);
            
            % Identify the Blocks to Inlcude
            idxBlocksToUse = idxStartNewAliquotUse + (0:MODE_ALIQUOT_LENGTH-1)'; % PIS blocks are excluded here: + 0:MODE_ALIQUOT_LENGTH-1 excludes the final block
            block_deltas = block_deltas(idxBlocksToUse,:);
            
            % Fill the Cell Array of Delta Values
            aliquot_deltas = reshape(block_deltas,MODE_ALIQUOT_LENGTH,NUM_ALIQUOTS)';
                        
            for fields = string(fieldnames(block_metadata))'
                block_metadata.(fields) = block_metadata.(fields)(idxBlocksToUse);
                aliquot_metadata.(fields) = reshape(block_metadata.(fields),MODE_ALIQUOT_LENGTH,NUM_ALIQUOTS)';
            end
            
        end
        
    end
end