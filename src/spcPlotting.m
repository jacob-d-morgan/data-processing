%% SPC Plotting
% Identifies the data from the South Pole Ice Core and plots various parts
% of the dataset.

%% Identify the relevant aliquots and their replicates
% Use only the aliquots that have 'SPICE' in the filename and than have a
% decimal point (indicating a bottom depth) in the ID1 string.

iSPC = contains(aliquot_metadata.filename(1,1,:),'SPICE') & contains(aliquot_metadata.ID1(1,1,:),'.');

workingStr = sprintf('%s_*',aliquot_metadata.ID1(1,1,iSPC));
workingStr(workingStr==' ') = '';
workingStr = replace(workingStr,'a_*','a*');
workingStr = replace(workingStr,'b_*','b*');
workingStr = replace(workingStr,'c_*','c*');
workingStr = replace(workingStr,'d_*','d*');

bottomDepthRep = sscanf(workingStr,'%f%c*',[2 Inf]);

rep = char(bottomDepthRep(2,:));
iRepA = rep=='a'; iRepB = rep=='b'; iRepC = rep=='c'; iRepD = rep=='d';
iNoRep = rep(~iRepA & ~iRepB & ~iRepC & ~iRepD);

