%% SPC Data Pull
% Identifies the data from the South Pole Ice Core from within the XP
% master dataset.

%% Identify the relevant aliquots and their replicates
% Use only the aliquots that have 'SPICE' in the filename and than have a
% decimal point (indicating a bottom depth) in the ID1 string.

iSPC = contains(aliquot_metadata.filename(1,1,:),'SPICE') & contains(aliquot_metadata.ID1(1,1,:),'.');

spc_aliquots = squeeze(aliquot_means_pisCorr_csCorr_ljaCorr(:,:,:,iSPC))';

% Split ID1 into bottom depth and replicate id
% Here I have to do some messing around to remove trailing spaces from ID1
% and to deal with samples that don't have a replicate id in ID1. I also
% assume that there are no replicates beyond rep d in the dataset.
workingStr = sprintf('%s_*',aliquot_metadata.ID1(1,1,iSPC));
workingStr(workingStr==' ') = '';
workingStr = replace(workingStr,'a_*','a*');
workingStr = replace(workingStr,'b_*','b*');
workingStr = replace(workingStr,'c_*','c*');
workingStr = replace(workingStr,'d_*','d*');
workingStr = replace(workingStr,'e_*','e*');

bottomDepthRep = (sscanf(workingStr,'%f%c*',[2 Inf])');
rep = string(char(bottomDepthRep(:,2)));
bottomDepth = bottomDepthRep(:,1);


% Make sure that there aren't any duplicate replicate ids for each sampling
% depth, i.e. make sure there aren't two samples labelled 890.17 a (hint:
% there is). Reassign replicate ids for these samples.
repIDs = ["a" "b" "c" "d" "e" "_"];
for ii=1:length(repIDs)-1
    iRep = rep==repIDs(ii);
    [~, idxUnique] = unique(bottomDepth(iRep),'stable');
    duplicateDepths = setdiff(1:sum(iRep),idxUnique);
    changeRep = rep(iRep); changeRep(duplicateDepths)=repIDs(ii+1);
    rep(iRep)=changeRep;
end

iRepA = rep=='a'; iRepB = rep=='b'; iRepC = rep=='c'; iRepD = rep=='d'; iRepE = rep=='e';


% Assign a replicate id to the aliquots that don't have one
% N.B. This approach does not take into account the date that the sample
% was run when assigning replicate id. For example, there is at least one
% instance of rep c being measured before rep a and b because the sample
% measured first was not given a replicate id and then, when it was
% remeasured a second and third time, these analyses were labelled a and b.

idxNoRep = find(~iRepA & ~iRepB & ~iRepC & ~iRepD & ~iRepE);
for ii = 1:length(idxNoRep)
    if ~ismember(bottomDepth(idxNoRep(ii)),bottomDepth(iRepA))
        iRepA(idxNoRep(ii)) = true;
        bottomDepthRep(idxNoRep(ii),2) = "a"; rep(idxNoRep(ii)) = "a";
    elseif ~ismember(bottomDepth(idxNoRep(ii)),bottomDepth(iRepB))
        iRepB(idxNoRep(ii)) = true;
        bottomDepthRep(idxNoRep(ii),2) = "b"; rep(idxNoRep(ii)) = "b";
    elseif ~ismember(bottomDepth(idxNoRep(ii)),bottomDepth(iRepC))
        iRepC(idxNoRep(ii)) = true;
        bottomDepthRep(idxNoRep(ii),2) = "c"; rep(idxNoRep(ii)) = "c";
    elseif ~ismember(bottomDepth(idxNoRep(ii)),bottomDepth(iRepD))
        iRepD(idxNoRep(ii)) = true;
        bottomDepthRep(idxNoRep(ii),2) = "d"; rep(idxNoRep(ii)) = "d";
    elseif ~ismember(bottomDepth(idxNoRep(ii)),bottomDepth(iRepE))
        iRepE(idxNoRep(ii)) = true;
        bottomDepthRep(idxNoRep(ii),2) = "e"; rep(idxNoRep(ii)) = "e";
    else
        warning('There is a Sample with an uncertain Replicate ID')
    end
end

if sum(~iRepA & ~iRepB & ~iRepC & ~iRepD & ~iRepE)>0
    warning('There is a Sample with an uncertain Replicate ID')
end



%% Assemble the replicates and replicate means into tables

repA = table();
repA.bottomDepth = bottomDepth(iRepA);
repA.rep_repA = rep(iRepA);
repA.d15N_repA = spc_aliquots(iRepA,4);
repA.d18O_repA = spc_aliquots(iRepA,5);
repA.d17O_repA = spc_aliquots(iRepA,6);
repA.d36Ar_repA = spc_aliquots(iRepA,7);
repA.d38Ar_repA = spc_aliquots(iRepA,8);
repA.dO2N2_repA = spc_aliquots(iRepA,9);
repA.dArN2_repA = spc_aliquots(iRepA,10);

repB = table();
repB.bottomDepth = bottomDepth(iRepB);
repB.rep_repB = rep(iRepB);
repB.d15N_repB = spc_aliquots(iRepB,4);
repB.d18O_repB = spc_aliquots(iRepB,5);
repB.d17O_repB = spc_aliquots(iRepB,6);
repB.d36Ar_repB = spc_aliquots(iRepB,7);
repB.d38Ar_repB = spc_aliquots(iRepB,8);
repB.dO2N2_repB = spc_aliquots(iRepB,9);
repB.dArN2_repB = spc_aliquots(iRepB,10);

repC = table();
repC.bottomDepth = bottomDepth(iRepC);
repC.rep_repC = rep(iRepC);
repC.d15N_repC = spc_aliquots(iRepC,4);
repC.d18O_repC = spc_aliquots(iRepC,5);
repC.d17O_repC = spc_aliquots(iRepC,6);
repC.d36Ar_repC = spc_aliquots(iRepC,7);
repC.d38Ar_repC = spc_aliquots(iRepC,8);
repC.dO2N2_repC = spc_aliquots(iRepC,9);
repC.dArN2_repC = spc_aliquots(iRepC,10);

repD = table();
repD.bottomDepth = bottomDepth(iRepD);
repD.rep_repD = rep(iRepD);
repD.d15N_repD = spc_aliquots(iRepD,4);
repD.d18O_repD = spc_aliquots(iRepD,5);
repD.d17O_repD = spc_aliquots(iRepD,6);
repD.d36Ar_repD = spc_aliquots(iRepD,7);
repD.d38Ar_repD = spc_aliquots(iRepD,8);
repD.dO2N2_repD = spc_aliquots(iRepD,9);
repD.dArN2_repD = spc_aliquots(iRepD,10);

repE = table();
repE.bottomDepth = bottomDepth(iRepE);
repE.rep_repE = rep(iRepE);
repE.d15N_repE = spc_aliquots(iRepE,4);
repE.d18O_repE = spc_aliquots(iRepE,5);
repE.d17O_repE = spc_aliquots(iRepE,6);
repE.d36Ar_repE = spc_aliquots(iRepE,7);
repE.d38Ar_repE = spc_aliquots(iRepE,8);
repE.dO2N2_repE = spc_aliquots(iRepE,9);
repE.dArN2_repE = spc_aliquots(iRepE,10);

spc_replicates = outerjoin(outerjoin(outerjoin(outerjoin(repA,repB,'MergeKeys',1),repC,'MergeKeys',1),repD,'MergeKeys',1),repE,'MergeKeys',1);

idxVars = find(contains(spc_replicates.Properties.VariableNames,'d15N'));
spc = table;
spc.bottomDepth = spc_replicates.bottomDepth;
spc.d15N = nanmean(spc_replicates{:,idxVars},2);
spc.d18O = nanmean(spc_replicates{:,idxVars+1},2);
spc.d17O = nanmean(spc_replicates{:,idxVars+2},2);
spc.d36Ar = nanmean(spc_replicates{:,idxVars+3},2);
spc.d38Ar = nanmean(spc_replicates{:,idxVars+4},2);
spc.dO2N2 = nanmean(spc_replicates{:,idxVars+5},2);
spc.dArN2 = nanmean(spc_replicates{:,idxVars+6},2);


%% Plot everything against depth
stackedFig(7,'XDir','reverse','Overlap',[0 20 0 20 0 20]);
stackedFigAx(); xlim([0 1800]);

stackedFigAx(1)
plot(spc_replicates.bottomDepth,[spc_replicates.d15N_repA spc_replicates.d15N_repB spc_replicates.d15N_repC spc_replicates.d15N_repD],'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d15N,'o','Color','none','MarkerFaceColor',lineCol(1))
ylabel('SPC \delta^{15}N [per mil]');

stackedFigAx(2)
plot(spc_replicates.bottomDepth,[spc_replicates.d18O_repA spc_replicates.d18O_repB spc_replicates.d18O_repC spc_replicates.d18O_repD],'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d18O,'o','Color','none','MarkerFaceColor',lineCol(2))
ylabel('SPC \delta^{18}O [per mil]');

stackedFigAx(3)
plot(spc_replicates.bottomDepth,[spc_replicates.d17O_repA spc_replicates.d17O_repB spc_replicates.d17O_repC spc_replicates.d17O_repD],'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d17O,'o','Color','none','MarkerFaceColor',lineCol(2)*0.7)
ylabel('SPC \delta^{17}O [per mil]');

stackedFigAx(4)
plot(spc_replicates.bottomDepth,[spc_replicates.d36Ar_repA spc_replicates.d36Ar_repB spc_replicates.d36Ar_repC spc_replicates.d36Ar_repD],'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d36Ar,'o','Color','none','MarkerFaceColor',lineCol(5))
ylabel('SPC \delta^{40/36}Ar [per mil]');

stackedFigAx(5)
plot(spc_replicates.bottomDepth,[spc_replicates.d38Ar_repA spc_replicates.d38Ar_repB spc_replicates.d38Ar_repC spc_replicates.d38Ar_repD],'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d38Ar,'o','Color','none','MarkerFaceColor',lineCol(5)*0.7)
ylabel('SPC \delta^{40/38}Ar [per mil]');

stackedFigAx(6)
plot(spc_replicates.bottomDepth,[spc_replicates.dO2N2_repA spc_replicates.dO2N2_repB spc_replicates.dO2N2_repC spc_replicates.dO2N2_repD],'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.dO2N2,'o','Color','none','MarkerFaceColor',lineCol(7))
ylabel('SPC \delta^O_2/N_2 [per mil]');

stackedFigAx(7)
plot(spc_replicates.bottomDepth,[spc_replicates.dArN2_repA spc_replicates.dArN2_repB spc_replicates.dArN2_repC spc_replicates.dArN2_repD],'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.dArN2,'o','Color','none','MarkerFaceColor',lineCol(4))
ylabel('SPC \delta^Ar/N_2 [per mil]');

stackedFigAx();
xlabel('Bottom Depth [m]')