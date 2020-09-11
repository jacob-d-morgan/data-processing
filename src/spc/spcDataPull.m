%% SPC Data Pull
% Identifies the data from the South Pole Ice Core from within the XP
% master sheet for 2016, 2017, and 2018.

clearvars; cluk;
clc;

%% Import Raw Data
% Imports the 'standard' raw dataset of aliquot delta values and their
% metadata. Does not include any cycles with a gas configuration different
% to Air+ or where the 28N2 beam is either saturated or absent, or cycles
% with a non-standard (different than the mode) number of cycles or blocks.

filesToImport = [
    "XP-2018(excelExportIntensityJDM).csv"; ...
    "XP-2017(excelExportIntensityJDM).csv"; ...
    "XP-2016(excelExportIntensityJDM).csv"; ...
    ];

spcMasterSheet = makeMasterSheet(filesToImport);


%% Identify the relevant aliquots and their replicates
% Use only the aliquots that have 'SPICE' in the filename and than have a
% decimal point (indicating a bottom depth) in the ID1 string.

iSPC = contains(spcMasterSheet.metadata.filename(:,1,1),'SPICE') & contains(spcMasterSheet.metadata.ID1(:,1,1),'.');

% Split ID1 into bottom depth and replicate id
exp2match = '(?<depth>\d*.\d*)\s*(?<rep>\w*)';
tokens = regexp(spcMasterSheet.metadata.ID1(iSPC,1,1),exp2match,'names');

bottomDepth = double(cellfun(@(x) x.depth,tokens));
rep = cellfun(@(x) x.rep,tokens);

%% Assemble Tables of Data
% Arrange the data in a table using the bottom depth as a grouping variable
% to identify replicates.

% Calculate the Sample Means for SPC
delta_names = string(spcMasterSheet.deltas_corr.Properties.VariableNames);
spc_aliquotMeans = varfun(@(x) mean(x,[3 4]),(spcMasterSheet.deltas_corr(iSPC,:)));
spc_aliquotMeans.Properties.VariableNames = delta_names;

% Identify Replicate Samples
spc.metadata = table;
[grp,spc.metadata.bottomDepth] = findgroups(bottomDepth);

numReps = nan(max(grp),1);
for ii_grp = grp'
    numReps(ii_grp) = sum(grp==ii_grp);
end    

% Build Tables of Sample Replicates and Means
spc.sampleRepsNoRej = table;
spc.sampleMeansNoRej = table;
for ii_delta = delta_names
    spc.sampleRepsNoRej.(ii_delta) = nan(max(grp),max(numReps));
    for jj_grp = grp'
        temp = spc_aliquotMeans{grp==jj_grp,ii_delta};
        spc.sampleRepsNoRej.(ii_delta)(jj_grp,1:length(temp)) = temp';
    end
end

for ii_delta = delta_names
    spc.sampleMeansNoRej.(ii_delta) = splitapply(@mean,spc_aliquotMeans.(ii_delta),grp);
end


%% Outlier Detection
% Flag all replicates that fall far from the running median for rejection.
% The isoutlier() function doesn't work well with the array of replicates
% so this is done in two steps: (1) find the rejection thresholds using the
% running median and MedAD of the array of replicate means, (2) apply the
% thresholds to the array of replicate samples.

load spice_ageModel.mat
spc.metadata.gasAge = interp1(spice_ageModel.depth,spice_ageModel.gasAge,spc.metadata.bottomDepth,'linear','extrap');

reject.low = table; reject.upp = table; reject.center = table;
for ii_delta = delta_names
    [~,reject.low.(ii_delta),reject.upp.(ii_delta),reject.center.(ii_delta)]=isoutlier(spc.sampleMeansNoRej.(ii_delta)(2:end),'movmed',1000,'SamplePoints',spc.metadata.gasAge(2:end));
end

%%%% CAUTION! These two lines are a workaround that deals with the NaN gas 
% age that occurs because the shallowest ice sample is in the firn 
% according to the age model so has a gas age of NaN. The shallowest sample
% is anomalously light in all measured ratios for both replicates, likely
% some sort of gas loss signal due to partially open bubbles
spc.metadata.gasAge(1) = 2014-1950; % Set to -64 yr BP (2014 CE, year that drilling started)
reject.low = [reject.low(1,:); reject.low]; 
reject.upp = [reject.upp(1,:); reject.upp]; 
reject.center = [reject.center(1,:); reject.center];


% Identify outlying replicates manually using the thresholds from isoutlier
spc.rejections = table;
for ii_delta = delta_names
    spc.rejections.(ii_delta) = spc.sampleRepsNoRej.(ii_delta) < reject.low.(ii_delta) | spc.sampleRepsNoRej.(ii_delta) > reject.upp.(ii_delta);
end

%% Plot everything against age
stackedFig(7,'XDir','reverse','Overlap',[0 20 0 20 0 20]);

for ii = 1:length(delta_names)
    stackedFigAx(ii)
    H=shadedErrorBar(spc.metadata.gasAge,reject.center.(delta_names(ii)),[reject.upp.(delta_names(ii))-reject.center.(delta_names(ii)) reject.center.(delta_names(ii))-reject.low.(delta_names(ii))]','-r'); 
    delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3;
    h1=plot(spc.metadata.gasAge,spc.sampleRepsNoRej.(delta_names(ii)),'.','Color',lineCol(ii));
    h2=plot(spc.metadata.gasAge,spc.sampleMeansNoRej.(delta_names(ii)),'-','Color',lineCol(ii)*0.7,'LineWidth',0.5);
    h3= plot(repmat(spc.metadata.gasAge,[max(numReps),1]),spc.sampleRepsNoRej.(delta_names(ii))(:),'o','Color','none','MarkerIndices',find(spc.rejections.(delta_names(ii))),'MarkerFaceColor',lineCol(10));
    ylabel(delta_names(ii))
end

legend(stackedFigAx(7),[h1(1) h2 H.mainLine h3],{'Replicates','Replicate Means','Moving Median','Outliers'},'Location','North','Orientation','Horizontal');

stackedFigAx;
xlim([-100 56000]);
xlabel('Age [yr BP]');
ticklabelformat(gca,'x','%.0f')

stackedFigReset;

%% Remove Outliers
% Set outliers to nan and recalculate depth means.

spc.sampleReps = spc.sampleRepsNoRej;
spc.sampleMeans = spc.sampleMeansNoRej;
for ii_delta = delta_names
    spc.sampleReps.(ii_delta)(spc.rejections.(ii_delta)) = nan;
    spc.sampleMeans.(ii_delta) = nanmean(spc.sampleReps.(ii_delta),2);
end


%% Plot everything against depth

stackedFig(7,'XDir','reverse','Overlap',[0 20 0 20 0 20]);
stackedFigAx(); xlim([0 1800]);

for ii = 1:length(delta_names)
    stackedFigAx(ii)
    plot(spc.metadata.bottomDepth,spc.sampleReps.(delta_names(ii)),'.','Color',lineCol(9));
    plot(spc.metadata.bottomDepth,spc.sampleMeans.(delta_names(ii)),'.','Color',lineCol(ii),'MarkerSize',10)
    ylabel(['SPC ' delta_names(ii) '[per mil]']);
end

stackedFigAx();
xlabel('Bottom Depth [m]')
stackedFigReset
