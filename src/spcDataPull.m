%% SPC Data Pull
% Identifies the data from the South Pole Ice Core from within the XP
% master dataset.

%% Identify the relevant aliquots and their replicates
% Use only the aliquots that have 'SPICE' in the filename and than have a
% decimal point (indicating a bottom depth) in the ID1 string.

iSPC = contains(aliquot_metadata.filename(:,1,1),'SPICE') & contains(aliquot_metadata.ID1(:,1,1),'.');

spc_aliquot_means = nanmean(mean(aliquot_deltas_pisCorr_csCorr_ljaCorr(iSPC,:,:,:),4),3);

% Split ID1 into bottom depth and replicate id
% Here I have to do some messing around to remove trailing spaces from ID1
% and to deal with samples that don't have a replicate id in ID1. I also
% assume that there are no replicates beyond rep d in the dataset.
workingStr = sprintf('%s_*',aliquot_metadata.ID1(iSPC,1,1));
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
repA.bottomDepth_repA = bottomDepth(iRepA);
repA.rep_repA = rep(iRepA);
repA.d15N_repA = spc_aliquot_means(iRepA,1);
repA.d18O_repA = spc_aliquot_means(iRepA,2);
repA.d17O_repA = spc_aliquot_means(iRepA,3);
repA.d36Ar_repA = spc_aliquot_means(iRepA,4);
repA.d38Ar_repA = spc_aliquot_means(iRepA,5);
repA.dO2N2_repA = spc_aliquot_means(iRepA,6);
repA.dArN2_repA = spc_aliquot_means(iRepA,7);

repB = table();
repB.bottomDepth = bottomDepth(iRepB);
repB.bottomDepth_repB = bottomDepth(iRepB);
repB.rep_repB = rep(iRepB);
repB.d15N_repB = spc_aliquot_means(iRepB,1);
repB.d18O_repB = spc_aliquot_means(iRepB,2);
repB.d17O_repB = spc_aliquot_means(iRepB,3);
repB.d36Ar_repB = spc_aliquot_means(iRepB,4);
repB.d38Ar_repB = spc_aliquot_means(iRepB,5);
repB.dO2N2_repB = spc_aliquot_means(iRepB,6);
repB.dArN2_repB = spc_aliquot_means(iRepB,7);

repC = table();
repC.bottomDepth = bottomDepth(iRepC);
repC.bottomDepth_repC = bottomDepth(iRepC);
repC.rep_repC = rep(iRepC);
repC.d15N_repC = spc_aliquot_means(iRepC,1);
repC.d18O_repC = spc_aliquot_means(iRepC,2);
repC.d17O_repC = spc_aliquot_means(iRepC,3);
repC.d36Ar_repC = spc_aliquot_means(iRepC,4);
repC.d38Ar_repC = spc_aliquot_means(iRepC,5);
repC.dO2N2_repC = spc_aliquot_means(iRepC,6);
repC.dArN2_repC = spc_aliquot_means(iRepC,7);

repD = table();
repD.bottomDepth = bottomDepth(iRepD);
repD.bottomDepth_repD = bottomDepth(iRepD);
repD.rep_repD = rep(iRepD);
repD.d15N_repD = spc_aliquot_means(iRepD,1);
repD.d18O_repD = spc_aliquot_means(iRepD,2);
repD.d17O_repD = spc_aliquot_means(iRepD,3);
repD.d36Ar_repD = spc_aliquot_means(iRepD,4);
repD.d38Ar_repD = spc_aliquot_means(iRepD,5);
repD.dO2N2_repD = spc_aliquot_means(iRepD,6);
repD.dArN2_repD = spc_aliquot_means(iRepD,7);

repE = table();
repE.bottomDepth = bottomDepth(iRepE);
repE.bottomDepth_repE = bottomDepth(iRepE);
repE.rep_repE = rep(iRepE);
repE.d15N_repE = spc_aliquot_means(iRepE,1);
repE.d18O_repE = spc_aliquot_means(iRepE,2);
repE.d17O_repE = spc_aliquot_means(iRepE,3);
repE.d36Ar_repE = spc_aliquot_means(iRepE,4);
repE.d38Ar_repE = spc_aliquot_means(iRepE,5);
repE.dO2N2_repE = spc_aliquot_means(iRepE,6);
repE.dArN2_repE = spc_aliquot_means(iRepE,7);

spc_replicates = outerjoin(outerjoin(outerjoin(outerjoin(repA,repB,'MergeKeys',1),repC,'MergeKeys',1),repD,'MergeKeys',1),repE,'MergeKeys',1);
spc_replicates = mergevars(spc_replicates,{'d15N_repA','d15N_repB','d15N_repC','d15N_repD','d15N_repE'},'NewVariableName','d15N');
spc_replicates = mergevars(spc_replicates,{'d18O_repA','d18O_repB','d18O_repC','d18O_repD','d18O_repE'},'NewVariableName','d18O');
spc_replicates = mergevars(spc_replicates,{'d17O_repA','d17O_repB','d17O_repC','d17O_repD','d17O_repE'},'NewVariableName','d17O');
spc_replicates = mergevars(spc_replicates,{'d36Ar_repA','d36Ar_repB','d36Ar_repC','d36Ar_repD','d36Ar_repE'},'NewVariableName','d36Ar');
spc_replicates = mergevars(spc_replicates,{'d38Ar_repA','d38Ar_repB','d38Ar_repC','d38Ar_repD','d38Ar_repE'},'NewVariableName','d38Ar');
spc_replicates = mergevars(spc_replicates,{'dO2N2_repA','dO2N2_repB','dO2N2_repC','dO2N2_repD','dO2N2_repE'},'NewVariableName','dO2N2');
spc_replicates = mergevars(spc_replicates,{'dArN2_repA','dArN2_repB','dArN2_repC','dArN2_repD','dArN2_repE'},'NewVariableName','dArN2');
spc_replicates = mergevars(spc_replicates,{'bottomDepth_repA','bottomDepth_repB','bottomDepth_repC','bottomDepth_repD','bottomDepth_repE'},'NewVariableName','bottomDepthReps');
spc_replicates = mergevars(spc_replicates,{'rep_repA','rep_repB','rep_repC','rep_repD','rep_repE'},'NewVariableName','repReps');

spc = table;
spc.bottomDepth = spc_replicates.bottomDepth;
spc.d15N = nanmean(spc_replicates.d15N,2);
spc.d18O = nanmean(spc_replicates.d18O,2);
spc.d17O = nanmean(spc_replicates.d17O,2);
spc.d36Ar = nanmean(spc_replicates.d36Ar,2);
spc.d38Ar = nanmean(spc_replicates.d38Ar,2);
spc.dO2N2 = nanmean(spc_replicates.dO2N2,2);
spc.dArN2 = nanmean(spc_replicates.dArN2,2);

%% Outlier Detection
load spice_ageModel.mat
spc.gasAge = interp1(spice_ageModel.depth,spice_ageModel.gasAge,spc.bottomDepth,'linear','extrap');
spc_replicates.gasAge = interp1(spice_ageModel.depth,spice_ageModel.gasAge,spc_replicates.bottomDepth,'linear','extrap');
spc_replicates.gasAgeReps = interp1(spice_ageModel.depth,spice_ageModel.gasAge,spc_replicates.bottomDepthReps,'linear','extrap');

low = table; upp = table; center = table;
[~,low.d15N,upp.d15N,center.d15N]=isoutlier(spc.d15N(2:end),'movmed',1000,'SamplePoints',spc.gasAge(2:end));
[~,low.d18O,upp.d18O,center.d18O]=isoutlier(spc.d18O(2:end),'movmed',1000,'SamplePoints',spc.gasAge(2:end));
[~,low.d17O,upp.d17O,center.d17O]=isoutlier(spc.d17O(2:end),'movmed',1000,'SamplePoints',spc.gasAge(2:end));
[~,low.d36Ar,upp.d36Ar,center.d36Ar]=isoutlier(spc.d36Ar(2:end),'movmed',1000,'SamplePoints',spc.gasAge(2:end));
[~,low.d38Ar,upp.d38Ar,center.d38Ar]=isoutlier(spc.d38Ar(2:end),'movmed',1000,'SamplePoints',spc.gasAge(2:end));
[~,low.dO2N2,upp.dO2N2,center.dO2N2]=isoutlier(spc.dO2N2(2:end),'movmed',1000,'SamplePoints',spc.gasAge(2:end));
[~,low.dArN2,upp.dArN2,center.dArN2]=isoutlier(spc.dArN2(2:end),'movmed',1000,'SamplePoints',spc.gasAge(2:end));

%%%% CAUTION! These two lines are a workaround that deals with the NaN gas 
% age that occurs because the shallowest ice sample is in the firn 
% according to the age model so has a gas age of NaN. The shallowest sample
% is anomalously light in all measured ratios for both replicates, likely
% some sort of gas loss signal due to partially open bubbles
spc.gasAge(1) = 2014-1950; spc_replicates.gasAge(1) = 1950-2014; spc_replicates.gasAgeReps(1,~isnan(spc_replicates.bottomDepthReps(1,:))) = 1950-2014; % Set to -64 yr BP (2014 CE, year that drilling started)
low = [low(1,:); low]; upp = [upp(1,:); upp]; center = [center(1,:); center];


% Identify outlying replicates manually using the thresholds from isoutlier
spc_replicateOutliers = table;
spc_replicateOutliers.d15N = spc_replicates.d15N < low.d15N | spc_replicates.d15N > upp.d15N;
spc_replicateOutliers.d18O = spc_replicates.d18O < low.d18O | spc_replicates.d18O > upp.d18O;
spc_replicateOutliers.d17O = spc_replicates.d17O < low.d17O | spc_replicates.d17O > upp.d17O;
spc_replicateOutliers.d36Ar = spc_replicates.d36Ar < low.d36Ar | spc_replicates.d36Ar > upp.d36Ar;
spc_replicateOutliers.d38Ar = spc_replicates.d38Ar < low.d38Ar | spc_replicates.d38Ar > upp.d38Ar;
spc_replicateOutliers.dO2N2 = spc_replicates.dO2N2 < low.dO2N2 | spc_replicates.dO2N2 > upp.dO2N2;
spc_replicateOutliers.dArN2 = spc_replicates.dArN2 < low.dArN2 | spc_replicates.dArN2 > upp.dArN2;


%% Plot everything against depth
stackedFig(7,'XDir','reverse','Overlap',[0 20 0 20 0 20]);
stackedFigAx(); xlim([-100 56000]);

stackedFigAx(1)
H=shadedErrorBar(spc.gasAge,center.d15N,[upp.d15N-center.d15N center.d15N-low.d15N]','-r'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3;
h1=plot(spc.gasAge,spc_replicates.d15N,'.','Color',lineCol(1));
h2=plot(spc.gasAge,spc.d15N,'-','Color',lineCol(1)*0.7,'LineWidth',0.5);
h3= plot(spc_replicates.gasAgeReps(spc_replicateOutliers.d15N),spc_replicates.d15N(spc_replicateOutliers.d15N),'o','Color','none','MarkerFaceColor',lineCol(10));
legend([h1(1) h2 H.mainLine h3],{'Replicates','Replicate Means','Moving Median','Outliers'},'Location','North','Orientation','Horizontal');
ylabel('\delta^{15}N [per mil]')

stackedFigAx(2)
H=shadedErrorBar(spc.gasAge,center.d18O,[upp.d18O-center.d18O center.d18O-low.d18O]','-r'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3;
h1=plot(spc.gasAge,spc_replicates.d18O,'.','Color',lineCol(3));
h2=plot(spc.gasAge,spc.d18O,'-','Color',lineCol(3)*0.7,'LineWidth',0.5);
h3= plot(spc_replicates.gasAgeReps(spc_replicateOutliers.d18O),spc_replicates.d18O(spc_replicateOutliers.d18O),'o','Color','none','MarkerFaceColor',lineCol(10));
%legend([h1(1) h2 H.mainLine h3],{'Replicates','Replicate Means','Moving Median','Outliers'},'Location','North','Orientation','Horizontal');
set(gca,'YDir','reverse');
ylabel('\delta^{18}O [per mil]')

stackedFigAx(3)
H=shadedErrorBar(spc.gasAge,center.d17O,[upp.d17O-center.d17O center.d17O-low.d17O]','-r'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3;
h1=plot(spc.gasAge,spc_replicates.d17O,'.','Color',lineCol(3)*0.7);
h2=plot(spc.gasAge,spc.d17O,'-','Color',lineCol(3)*0.7*0.7,'LineWidth',0.5);
h3= plot(spc_replicates.gasAgeReps(spc_replicateOutliers.d17O),spc_replicates.d17O(spc_replicateOutliers.d17O),'o','Color','none','MarkerFaceColor',lineCol(10));
%legend([h1(1) h2 H.mainLine h3],{'Replicates','Replicate Means','Moving Median','Outliers'},'Location','North','Orientation','Horizontal');
set(gca,'YDir','reverse');
ylabel('\delta^{17}O [per mil]')

stackedFigAx(4)
H=shadedErrorBar(spc.gasAge,center.d36Ar,[upp.d36Ar-center.d36Ar center.d36Ar-low.d36Ar]','-r'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3;
h1=plot(spc.gasAge,spc_replicates.d36Ar,'.','Color',lineCol(5));
h2=plot(spc.gasAge,spc.d36Ar,'-','Color',lineCol(5),'LineWidth',0.5);
h3= plot(spc_replicates.gasAgeReps(spc_replicateOutliers.d36Ar),spc_replicates.d36Ar(spc_replicateOutliers.d36Ar),'o','Color','none','MarkerFaceColor',lineCol(10));
%legend([h1(1) h2 H.mainLine h3],{'Replicates','Replicate Means','Moving Median','Outliers'},'Location','North','Orientation','Horizontal');
ylabel('\delta^{40}/_{36}Ar [per mil]')

stackedFigAx(5)
H=shadedErrorBar(spc.gasAge,center.d38Ar,[upp.d38Ar-center.d38Ar center.d38Ar-low.d38Ar]','-r'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3;
h1=plot(spc.gasAge,spc_replicates.d38Ar,'.','Color',lineCol(5)*0.7);
h2=plot(spc.gasAge,spc.d38Ar,'-','Color',lineCol(5)*0.7*0.7,'LineWidth',0.5);
h3= plot(spc_replicates.gasAgeReps(spc_replicateOutliers.d38Ar),spc_replicates.d38Ar(spc_replicateOutliers.d38Ar),'o','Color','none','MarkerFaceColor',lineCol(10));
%legend([h1(1) h2 H.mainLine h3],{'Replicates','Replicate Means','Moving Median','Outliers'},'Location','North','Orientation','Horizontal');
ylabel('\delta^{40}/_{38}Ar [per mil]')

stackedFigAx(6)
H=shadedErrorBar(spc.gasAge,center.dO2N2,[upp.dO2N2-center.dO2N2 center.dO2N2-low.dO2N2]','-r'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3;
h1=plot(spc.gasAge,spc_replicates.dO2N2,'.','Color',lineCol(7));
h2=plot(spc.gasAge,spc.dO2N2,'-','Color',lineCol(7)*0.7,'LineWidth',0.5);
h3= plot(spc_replicates.gasAgeReps(spc_replicateOutliers.dO2N2),spc_replicates.dO2N2(spc_replicateOutliers.dO2N2),'o','Color','none','MarkerFaceColor',lineCol(10));
%legend([h1(1) h2 H.mainLine h3],{'Replicates','Replicate Means','Moving Median','Outliers'},'Location','North','Orientation','Horizontal');
ylabel('\delta^{40}/_{38}Ar [per mil]')

stackedFigAx(7)
H=shadedErrorBar(spc.gasAge,center.dArN2,[upp.dArN2-center.dArN2 center.dArN2-low.dArN2]','-r'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3;
h1=plot(spc.gasAge,spc_replicates.dArN2,'.','Color',lineCol(4));
h2=plot(spc.gasAge,spc.dArN2,'-','Color',lineCol(4)*0.7,'LineWidth',0.5);
h3= plot(spc_replicates.gasAgeReps(spc_replicateOutliers.dArN2),spc_replicates.dArN2(spc_replicateOutliers.dArN2),'o','Color','none','MarkerFaceColor',lineCol(10));
%legend([h1(1) h2 H.mainLine h3],{'Replicates','Replicate Means','Moving Median','Outliers'},'Location','North','Orientation','Horizontal');
ylabel('\delta^{40}/_{38}Ar [per mil]')


stackedFigAx;
xlabel('Age [yr BP]');
ticklabelformat(gca,'x','%.0f')

stackedFigReset;

%% Remove Outliers
% Set outliers to nan and recalculate depth means.

spc_replicates.d15N(spc_replicateOutliers.d15N)=nan;
spc_replicates.d18O(spc_replicateOutliers.d18O)=nan;
spc_replicates.d17O(spc_replicateOutliers.d17O)=nan;
spc_replicates.d36Ar(spc_replicateOutliers.d36Ar)=nan;
spc_replicates.d38Ar(spc_replicateOutliers.d38Ar)=nan;
spc_replicates.dO2N2(spc_replicateOutliers.dO2N2)=nan;
spc_replicates.dArN2(spc_replicateOutliers.dArN2)=nan;

spc.d15N = nanmean(spc_replicates.d15N,2);
spc.d18O = nanmean(spc_replicates.d18O,2);
spc.d17O = nanmean(spc_replicates.d17O,2);
spc.d36Ar = nanmean(spc_replicates.d36Ar,2);
spc.d38Ar = nanmean(spc_replicates.d38Ar,2);
spc.dO2N2 = nanmean(spc_replicates.dO2N2,2);
spc.dArN2 = nanmean(spc_replicates.dArN2,2);

%% Plot everything against depth

stackedFig(7,'XDir','reverse','Overlap',[0 20 0 20 0 20]);
stackedFigAx(); xlim([0 1800]);

stackedFigAx(1)
plot(spc_replicates.bottomDepth,spc_replicates.d15N,'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d15N,'o','Color','none','MarkerFaceColor',lineCol(1))
ylabel('SPC \delta^{15}N [per mil]');

stackedFigAx(2)
plot(spc_replicates.bottomDepth,spc_replicates.d18O,'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d18O,'o','Color','none','MarkerFaceColor',lineCol(3))
ylabel('SPC \delta^{18}O [per mil]');

stackedFigAx(3)
plot(spc_replicates.bottomDepth,spc_replicates.d17O,'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d17O,'o','Color','none','MarkerFaceColor',lineCol(3)*0.7)
ylabel('SPC \delta^{17}O [per mil]');

stackedFigAx(4)
plot(spc_replicates.bottomDepth,spc_replicates.d36Ar,'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d36Ar,'o','Color','none','MarkerFaceColor',lineCol(5))
ylabel('SPC \delta^{40}/_{36}Ar [per mil]');

stackedFigAx(5)
plot(spc_replicates.bottomDepth,spc_replicates.d38Ar,'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.d38Ar,'o','Color','none','MarkerFaceColor',lineCol(5)*0.7)
ylabel('SPC \delta^{40}/_{38}Ar [per mil]');

stackedFigAx(6)
plot(spc_replicates.bottomDepth,spc_replicates.dO2N2,'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.dO2N2,'o','Color','none','MarkerFaceColor',lineCol(7))
ylabel('SPC \deltaO_2/N_2 [per mil]');

stackedFigAx(7)
plot(spc_replicates.bottomDepth,spc_replicates.dArN2,'.','Color',lineCol(9));
plot(spc.bottomDepth,spc.dArN2,'o','Color','none','MarkerFaceColor',lineCol(4))
ylabel('SPC \deltaAr/N_2 [per mil]');

stackedFigAx();
xlabel('Bottom Depth [m]')
stackedFigReset