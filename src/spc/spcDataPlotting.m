%% SPC Data Plotting
% Plots figures from the SPC XP dataset.

clearvars;
cluk; clc

%% Import Data
% Load the relevant csv files and compile the data into the SPC dataset.

filesToImport = [
    "XP-2018(excelExportIntensityJDM).csv"; ...
    "XP-2017(excelExportIntensityJDM).csv"; ...
    "XP-2016(excelExportIntensityJDM).csv"; ...
    ];

spcMasterSheet = makeMasterSheet(filesToImport);
exp2match = '^(?<depth>\d\d\d*(\.?)\d*)\s*(?<rep>\w?)$'; % Beginning; depth: two or more numeric chars, maybe with a decimal point; maybe a space; rep: maybe a single alphanumeric (or underscore) char; End.
spc = makeIceCoreDataset(spcMasterSheet,'SPICE',exp2match);

%% Outlier Detection
% Flag all replicates that fall far from the running median for rejection.
% The isoutlier() function doesn't work well with the array of replicates
% so this is done in two steps: (1) find the rejection thresholds using the
% running median and MedAD of the array of replicate means, (2) apply the
% thresholds to the array of replicate samples.

load spice_ageModel.mat

spc.metadata.gasAge = interp1(spice_ageModel.depth,spice_ageModel.gasAge,spc.metadata.bottomDepth,'linear','extrap');
spc.metadata.gasAge(1) = 1950-2014; % Set to -64 yr BP (2014 CE, year that drilling started)

spc.rejections = table;
delta_names = string(spc.sampleMeansNoRej.Properties.VariableNames);
reject.low = table; reject.upp = table; reject.cen = table;
for ii_delta = delta_names
    reject.cen.(ii_delta) = movmedian(spc.sampleMeansNoRej.(ii_delta),250,'omitnan','SamplePoints',spc.metadata.gasAge);
    dev = spc.sampleRepsNoRej.(ii_delta) - reject.cen.(ii_delta);
    dev = dev(~isnan(dev)); % Muse omit nan values for CDF to end at 1.0
    
    [CDF,edges] = histcounts(dev(:),'BinMethod','fd','Normalization','cdf');
    
    reject.low.(ii_delta) = reject.cen.(ii_delta) + edges(find(CDF>0.01,1,'first'));
    reject.upp.(ii_delta) = reject.cen.(ii_delta) + edges(find(CDF>0.99,1,'first'));
    iRej = (spc.sampleRepsNoRej.(ii_delta) > reject.upp.(ii_delta)) | (spc.sampleRepsNoRej.(ii_delta) < reject.low.(ii_delta));
    
    spc.rejections.(ii_delta) = iRej;
end

%%%% CAUTION! These two lines are a workaround that deals with the NaN gas 
% age that occurs because the shallowest ice sample is in the firn 
% according to the age model so has a gas age of NaN. The shallowest sample
% is anomalously light in all measured ratios for both replicates, likely
% some sort of gas loss signal due to partially open bubbles
% spc.metadata.gasAge(1) = 2014-1950; % Set to -64 yr BP (2014 CE, year that drilling started)
% reject.low = [reject.low(1,:); reject.low]; 
% reject.upp = [reject.upp(1,:); reject.upp]; 
% reject.center = [reject.center(1,:); reject.center];


%% Plot everything against age
stackedFig(numel(delta_names),'XDir','reverse');

for ii = 1:length(delta_names)
    stackedFigAx(ii)
    H=shadedErrorBar(spc.metadata.gasAge,reject.cen.(delta_names(ii)),[reject.upp.(delta_names(ii))-reject.cen.(delta_names(ii)) reject.cen.(delta_names(ii))-reject.low.(delta_names(ii))]','-r'); 
    delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3;
    h1=plot(spc.metadata.gasAge,spc.sampleRepsNoRej.(delta_names(ii)),'.','Color',lineCol(ii));
    h2=plot(spc.metadata.gasAge,spc.sampleMeansNoRej.(delta_names(ii)),'-','Color',lineCol(ii)*0.7,'LineWidth',0.5);
    h3= plot(repmat(spc.metadata.gasAge,[max(spc.metadata.numRepsNoRej),1]),spc.sampleRepsNoRej.(delta_names(ii))(:),'o','Color','none','MarkerIndices',find(spc.rejections.(delta_names(ii))),'MarkerFaceColor',lineCol(10));
    ylabel(delta_names(ii))
end

legend(stackedFigAx(numel(delta_names)),[h1(1) h2 H.mainLine h3],{'Replicates','Replicate Means','Moving Median','Outliers'},'Location','North','Orientation','Horizontal');

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

stackedFig(numel(delta_names),'XDir','reverse');
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

%% Make Gravitational Correction
% Correct the isotopic and elemental ratios for gravitational enrichment in
% the firn column, using the measured d15N.

massDiff = table;
massDiff.d15N = 0;
massDiff.d18O = 2;
massDiff.d17O = 1;
massDiff.d4036Ar = 4;
massDiff.d4038Ar = 2;
massDiff.dO2N2 = 4;
massDiff.dArN2 = 12;

delta_names_grav = append(delta_names(2:end),'grav');
for ii = 1:numel(delta_names_grav)
    spc.sampleReps.(delta_names_grav(ii)) = spc.sampleReps.(delta_names(ii+1)) - massDiff.(delta_names(ii+1))*spc.sampleReps.d15N;
    spc.sampleMeans.(delta_names_grav(ii)) = nanmean(spc.sampleReps.(delta_names_grav(ii)),2);
end

%% Calculate Pair Differences

spc.pairDiffs = table;  
delta_names = string(spc.sampleReps.Properties.VariableNames);
for ii = 1:numel(delta_names)
    spc.pairDiffs.(delta_names(ii)) = diff(spc.sampleReps.(delta_names(ii)),1,2);
end


%% Pair Difference Plots

iBub = spc.metadata.bottomDepth < 700;
iBCTZ = spc.metadata.bottomDepth > 700 & spc.metadata.bottomDepth < 1140;
iClath = spc.metadata.bottomDepth > 1140;

figure; hold on;
P=gmregress(spc.sampleMeans.dO2N2grav(iClath),spc.sampleMeans.dArN2grav(iClath));
plot(spc.sampleMeans.dO2N2grav,polyval(P,spc.sampleMeans.dO2N2grav),'-k')
plot(spc.sampleMeans.dO2N2grav(iClath,:),spc.sampleMeans.dArN2grav(iClath,:),'.','Color',lineCol(4));
xlabel('\deltaO_2/N_2 [per mil]'); ylabel('\deltaAr/N_2 [per mil]');
title('Ar/N_2 vs O_2/N_2 in SPC Clathrate Ice')

xx = spc.pairDiffs.dO2N2grav;
varsToPlot = ["dArN2grav" "d15N" "d18Ograv" "d17Ograv" "d4036Argrav" "d4038Argrav"];
numVars = length(varsToPlot);
depthsToPlot = [iBub iBCTZ any([iBub,iBCTZ],2) iClath]; numDepths = size(depthsToPlot,2);
cols = [lineCol(4); lineCol(1); lineCol(3); lineCol(3)*0.7; lineCol(5); lineCol(5)*0.7];
axCount = 1;

figure;
for ii = 1:numDepths
    for jj = 1:numVars
        
        x = xx(depthsToPlot(:,ii),:);
        y = spc.pairDiffs.(varsToPlot(jj))(depthsToPlot(:,ii),:);
        
        %x = xx(~isnan(xx) & ~isnan(yy));
        %y = yy(~isnan(xx) & ~isnan(yy));
        
        
        [R,P] = corrcoef([x(:) y(:)],'Rows','Complete');
        r(ii,jj) = R(1,2);
        p(ii,jj) = P(1,2);
        
        
        subplot(numDepths,numVars,axCount); hold on;
        if p(ii,jj) < 0.1
            [pFit,~,CI,pVal] = gmregress(x,y,0.001);
            plot(x,polyval(pFit,x),'-k');
            text(-8,0.75*min(get(gca,'YLim')),{['m = ' sprintf('%.2f',pFit(1)*1000) ' \pm ' sprintf('%.2f',(pFit(1)-CI(1,1))*1000) ' per meg/per mil'];['p = ' num2str(pVal)]});
        end
        
        plot(x,y,'.','Color',cols(jj,:))
        xlabel('\Delta\deltaO_2/N_2 [per mil]');
        ylabel(varsToPlot(jj));
        
        axCount=axCount+1;
    end
end
