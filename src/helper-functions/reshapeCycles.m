function [aliquot_deltas,aliquot_metadata] = reshapeCycles(cycle_deltas,cycle_metadata)
% RESHAPECYCLES converts cycle deltas to an aliquot-delta-block-cycle array
%   ALIQUOTDELTAS = RESHAPECYCLES(CYCLEDELTAS) reshapes CYCLEDELTAS to the
%   array ALIQUOTDELTAS, which contains the elements of CYCLEDELTAS in a
%   aliquot-by-delta-by-block-by-cycle structure.


% -------------------------------------------------------------------------

% Identify the New Blocks
[~,idxStartNewBlockAll,~] = unique(cycle_metadata.filename,'stable');
blockLengthsAll = diff(idxStartNewBlockAll);
blockLengthsAll(end+1) = length(cycle_deltas) - (idxStartNewBlockAll(end)-1);

% TAKE ONLY THE BLOCKS WITH THE MOST COMMON NUMBER OF CYCLES!
MODE_BLOCK_LENGTH = mode(blockLengthsAll);
idxStartNewBlockUse = idxStartNewBlockAll(blockLengthsAll==MODE_BLOCK_LENGTH);
NUM_BLOCKS = length(idxStartNewBlockUse);

% Identify the Cycles to Include
iCyclesToUse = ismember(cycle_metadata.filename,cycle_metadata.filename(idxStartNewBlockUse));
cycle_deltas = cycle_deltas(iCyclesToUse,:);

% Fill the Array of Delta Values
block_deltas = reshape(cycle_deltas,MODE_BLOCK_LENGTH,NUM_BLOCKS,size(cycle_deltas,2));

metadata_fields = fieldnames(cycle_metadata);
for ii = 1:numel(metadata_fields)
    cycle_metadata.(metadata_fields{ii}) = cycle_metadata.(metadata_fields{ii})(iCyclesToUse,:);
    block_metadata.(metadata_fields{ii}) = squeeze(reshape(cycle_metadata.(metadata_fields{ii}),MODE_BLOCK_LENGTH,NUM_BLOCKS,size(cycle_metadata.(metadata_fields{ii}),2)));
end



% Identify the New Aliquots
% Define Methods that Are Run at the Start of Each New Sample
SAMPLE_START_METHODS = ["can_v_can"; "Automation_SA_Delay"; "Automation_SA"];

% Identify aliquots that were run with one of these three methods
iWorking = false(1,size(block_metadata.method,2));
for ii=1:length(SAMPLE_START_METHODS)
    iWorking = iWorking | block_metadata.method(1,:)==SAMPLE_START_METHODS(ii);
end

% Identify aliquots also by finding the start of a new sequence (value of sequence row decreases).
idxStartNewAliquotAll = find(iWorking | ([0 diff(block_metadata.sequenceRow(1,:))]<0));
aliquotLengthsAll = diff(idxStartNewAliquotAll);
aliquotLengthsAll(end+1)=length(block_deltas)-(idxStartNewAliquotAll(end)-1);

% TAKE ONLY THE ALIQUOTS WITH THE MODE COMMON NUMBER OF BLOCKS!
% At this stage, it is necessary to also include aliquots with one block
% more than the mode so as not to exclude PIS aliquots. The PIS blocks are
% not included in the final output as they are filtered out below.
MODE_ALIQUOT_LENGTH = mode(aliquotLengthsAll);
idxStartNewAliquotUse = idxStartNewAliquotAll(aliquotLengthsAll==MODE_ALIQUOT_LENGTH | aliquotLengthsAll==(MODE_ALIQUOT_LENGTH+1));
NUM_ALIQUOTS = length(idxStartNewAliquotUse);

% Identify the Blocks to Inlcude
idxBlocksToUse = idxStartNewAliquotUse + [0:3]'; % PIS blocks are excluded here: + [0:3] includes only the first through fourth blocks
block_deltas = block_deltas(:,idxBlocksToUse,:);

% Fill the Array of Delta Values
aliquot_deltas = reshape(block_deltas,MODE_BLOCK_LENGTH,MODE_ALIQUOT_LENGTH,NUM_ALIQUOTS,size(cycle_deltas,2));
aliquot_deltas = permute(aliquot_deltas,[3 4 2 1]);

for ii = 1:numel(metadata_fields)
    block_metadata.(metadata_fields{ii}) = block_metadata.(metadata_fields{ii})(:,idxBlocksToUse);
    aliquot_metadata.(metadata_fields{ii}) = squeeze(permute(reshape(block_metadata.(metadata_fields{ii}),MODE_BLOCK_LENGTH,MODE_ALIQUOT_LENGTH,NUM_ALIQUOTS,size(cycle_metadata.(metadata_fields{ii}),3)),[3 4 2 1]));
end


end