%% SPC Plotting
% Identifies the data from the South Pole Ice Core and plots various parts
% of the dataset.

%% Identify the relevant aliquots and their replicates
% Use only the aliquots that have 'SPICE' in the filename and than have a
% decimal point (indicating a bottom depth) in the ID1 string.

iSPC = contains(aliquot_metadata.filename(1,1,:),'SPICE') & contains(aliquot_metadata.ID1(1,1,:),'.');

%Split ID1 into bottom depth and replicate id
% Here I have to do some messing around to remove trailing spaces from ID1
% and to deal with samples that don't have a replicate id in ID1. I also
% assume that there are no replicates beyond rep d in the dataset.
workingStr = sprintf('%s_*',aliquot_metadata.ID1(1,1,iSPC));
workingStr(workingStr==' ') = '';
workingStr = replace(workingStr,'a_*','a*');
workingStr = replace(workingStr,'b_*','b*');
workingStr = replace(workingStr,'c_*','c*');
workingStr = replace(workingStr,'d_*','d*');

bottomDepthRep = sscanf(workingStr,'%f%c*',[2 Inf]);
rep = char(bottomDepthRep(2,:));
bottomDepth = bottomDepthRep(1,:);

iRepA = rep=='a'; iRepB = rep=='b'; iRepC = rep=='c'; iRepD = rep=='d';

% Assign a replicate id to the aliquots that don't have one
% N.B. This approach does not take into account the date that the sample 
% was run when assigning replicate id. For example, there is at least one
% instance of rep c being measured before rep a and b because the sample
% measured first was not given a replicate id and then, when it was
% remeasured a second and third time, these analyses were labelled a and b.

idxNoRep = find(~iRepA & ~iRepB & ~iRepC & ~iRepD);
for ii = 1:length(idxNoRep)
    if ~ismember(bottomDepth(idxNoRep(ii)),bottomDepth(iRepA))
        iRepA(idxNoRep(ii)) = true;
        bottomDepthRep(2,idxNoRep(ii)) = 'a'; rep(idxNoRep(ii)) = 'a';
    elseif ~ismember(bottomDepth(idxNoRep(ii)),bottomDepth(iRepB))
        iRepB(idxNoRep(ii)) = true;
        bottomDepthRep(2,idxNoRep(ii)) = 'b'; rep(idxNoRep(ii)) = 'b';
    elseif ~ismember(bottomDepth(idxNoRep(ii)),bottomDepth(iRepC))
        iRepC(idxNoRep(ii)) = true;
        bottomDepthRep(2,idxNoRep(ii)) = 'c'; rep(idxNoRep(ii)) = 'c';
    elseif ~ismember(bottomDepth(idxNoRep(ii)),bottomDepth(iRepD))
        iRepD(idxNoRep(ii)) = true;
        bottomDepthRep(2,idxNoRep(ii)) = 'd'; rep(idxNoRep(ii)) = 'd';
    else
        warning('There is a Sample with an uncertain Replicate ID')
    end
end

if sum(~iRepA & ~iRepB & ~iRepC & ~iRepD)>0
    warning('There is a Sample with an uncertain Replicate ID')
end



        
        
