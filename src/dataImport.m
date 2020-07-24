%% dataImport %%

% clearvars;
cluk; clc;

% set(0,'defaultFigureVisible','off'); disp('Run dataImport: Turning Figures Off');

%% Import Data
% Imports the 'standard' raw dataset of aliquot delta values and their
% metadata. Does not include any cycles with a gas configuration
% different to Air+ or where the 28N2 beam is saturated or absent, or
% cycles with a non-standard number of cycles or blocks.
% Rather than reading in the files each time, it's easier to read them in
% once using the makeRawDataset() line and then use the clearvars -EXCEPT
% line to preserve the variables read in by this line for future runs.
% Toggle the comments on the indicated line(s) to adjust whether or not the
% files are re-read.

clearvars -EXCEPT aliquot_deltas metadata aliquot_deltas_pis aliquot_metadata_pis;

filesToImport = [
    "XP-2018(excelExportIntensityJDM).csv"; ...
    "XP-2017(excelExportIntensityJDM).csv"; ...
    "XP-2016(excelExportIntensityJDM).csv"
    ];

% Generate 'Standard' Raw Dataset
% [aliquot_deltas,metadata,aliquot_deltas_pis,aliquot_metadata_pis] = makeRawDataset(filesToImport,'includePIS',true); % <-- TOGGLE MY COMMENT

aliquot_metadata = metadata.metadata;
delta_names = metadata.delta_names;
delta_labels = metadata.delta_labels;
delta_units = metadata.delta_units;

metadata_fields = string(fieldnames(aliquot_metadata))';


%% Make PIS Correction
% Correct the delta values in aliquot_deltas for the effect of pressure
% imbalance in the bellows.

% Calculate PIS Values
[calcPis,pisStats] = calculatePisValues(aliquot_deltas,aliquot_metadata,aliquot_deltas_pis,aliquot_metadata_pis);

% Make PIS Correction
calcPis(pisStats.rejections) = nan;
[aliquot_deltas_pisCorr,PIS] = makePisCorr(aliquot_deltas,aliquot_metadata.msDatetime,aliquot_metadata.pressureImbal,calcPis);

%% Plot a time-series of the PIS and related parameters
% This is useful to identify aliquots where the PIS block did no run
% correctly or where the r-squared is low, suggesting a poor determination
% of the PIS. These cases can be filtered out below.

% Import Mass Spec Events to Annotate Figures
massSpecEvents = readtable('spreadsheet_metadata.xlsx','Sheet',1);
massSpecEvents.Event = categorical(massSpecEvents.Event);

iPIS = ~isnan(aliquot_deltas_pis(:,:,1,1));
if sum(any(iPIS,2)) ~= sum(all(iPIS,2)) % Check that all delta values identify each PIS experiment
    warning('Warning: Some delta values are misisng a PIS block for one or more experiments')
    iPIS = any(iPIS,2); % If some delta values are missing a PIS block somehow, calculate the PIS for the other delta values anyway
else
    iPIS = iPIS(:,1); % Otherwise, just make iPis a vector (from a matrix) by taking the first column.
end

stackedFig(3,'RelSize',[0.4 1.7 0.9],'Overlap',[-10 -10]);
stackedFigAx
xlim(datenum(['01 Jan 2016'; '31 Dec 2018']))

% Plot R-Squared of each PIS Experiment
stackedFigAx(1)
for ii = 1:numel(delta_names)
    set(gca,'ColorOrderIndex',ii)
    plot(aliquot_metadata.msDatenum(iPIS,1,1),pisStats.rSq(iPIS,ii),'-ok','MarkerIndices',find(pisStats.rejections(iPIS,ii)),'MarkerFaceColor',lineCol(ii)); % Plot all the r-squared values with markers for rejected values 
    legH(ii)=plot(aliquot_metadata.msDatenum(~pisStats.rejections(:,ii) & iPIS,1,1),pisStats.rSq(~pisStats.rejections(:,ii) & iPIS,ii),'s-','MarkerFaceColor',lineCol(ii)); % Replot as overlay, omitting rejected aliquots
end
legend(legH,delta_names,'Orientation','Horizontal','Location','South')
ylabel('r^2');
ylim([0.7 1]);

% Plot PIS Value for each PIS Experiment
stackedFigAx(2)
for ii = 1:numel(delta_names)
    set(gca,'ColorOrderIndex',ii)
    plot(aliquot_metadata.msDatenum(iPIS,1,1),pisStats.measuredPis(iPIS,ii),'-ok','MarkerIndices',find(pisStats.rejections(iPIS,ii)),'MarkerFaceColor',lineCol(ii)); % Plot all the r-squared values with markers for rejected values 
    legH(ii)=plot(aliquot_metadata.msDatenum(~pisStats.rejections(:,ii) & iPIS,1,1),pisStats.measuredPis(~pisStats.rejections(:,ii) & iPIS,ii),'s-','MarkerFaceColor',lineCol(ii)); % Replot as overlay, omitting rejected aliquots
end
ylabel('PIS [per mil/per mil]');
ylim([-0.01 0.005]);

% Plot Pressure Imbalance for each PIS Experiment
stackedFigAx(3)
plot(aliquot_metadata.msDatenum(iPIS & any(~pisStats.rejections,2),1,1),pisStats.pImbal(iPIS & any(~pisStats.rejections,2),:),'-^','Color',lineCol(9));
text(aliquot_metadata.msDatenum(iPIS & any(~pisStats.rejections,2),1,1),pisStats.pImbal(iPIS & any(~pisStats.rejections,2)),aliquot_metadata_pis.ID1(iPIS & any(~pisStats.rejections,2),1,1));
ylabel('Pressure Imbalance [per mil]');
ylim([-600 0])

% Set Labels & Limits etc.
stackedFigAx();
xlim(datenum(['01 Jan 2016'; '31 Dec 2018']))
datetick('x','keeplimits');
drawnow;
stackedFigReset

% Plot Time Series of PIS Value Used for Each Aliquot
stackedFig(numel(delta_names))
for ii=1:numel(delta_names)
    stackedFigAx(ii)
    plot(aliquot_metadata.msDatenum(~pisStats.rejections(:,ii),1,1),pisStats.measuredPis(~pisStats.rejections(:,ii),ii),'o','Color','none','MarkerFaceColor',lineCol(ii))
    plot(aliquot_metadata.msDatenum(:,1,1),PIS(:,ii,1,1),'.','Color',lineCol(ii)*0.5)
    ylabel(delta_labels{ii})
end

% Add Date of Filament Changes/Refocusing
stackedFigAx

iPlot = massSpecEvents.Event == "New Filament" | massSpecEvents.Event == "Refocus";
plot(datenum([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)']),repmat(get(gca,'YLim')',1,sum(iPlot)),'-r');
text(datenum(massSpecEvents.EndDate(iPlot)),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color','r')

% Set Labels & Limits etc.
title('PIS Values used for Correction')
xlabel('Date')
xlim(datenum(["01-Jan-2016" "01-Jan-2019"]))
datetick('x','keeplimits')
drawnow;
stackedFigReset


%% Make Chemical Slope Correction
% Still need to:
%   1) Weed out the questionable blocks/aliquots that mess up some of the
%      chem slope experiments
%   2) Figure out how I'm going to do the CS Correction for Ar isotopes

% Calculate CS Values
iCS = contains(aliquot_metadata.ID1(:,1,1),'CS');
iCS_AddO2 = iCS & contains(aliquot_metadata.ID1(:,1,1),{'15','N'});
iCS_AddN2 = iCS & contains(aliquot_metadata.ID1(:,1,1),{'18','O'});

[CS_15N,csStats_15N] = calculateChemSlope(aliquot_deltas_pisCorr(:,6,:,:),aliquot_deltas_pisCorr(:,1,:,:),aliquot_metadata,iCS_AddO2);
[CS_ArN2,csStats_ArN2] = calculateChemSlope(aliquot_deltas_pisCorr(:,6,:,:),aliquot_deltas_pisCorr(:,7,:,:),aliquot_metadata,iCS_AddO2);
[CS_18O,csStats_18O] = calculateChemSlope((1/(aliquot_deltas_pisCorr(:,6,:,:)/1000+1)-1)*1000,aliquot_deltas_pisCorr(:,2,:,:),aliquot_metadata,iCS_AddN2);
[CS_17O,csStats_17O] = calculateChemSlope((1/(aliquot_deltas_pisCorr(:,6,:,:)/1000+1)-1)*1000,aliquot_deltas_pisCorr(:,3,:,:),aliquot_metadata,iCS_AddN2);

x_temp = [((aliquot_deltas_pisCorr(:,7,:,:)./1000+1).^-1-1)*1000 ((aliquot_deltas_pisCorr(:,6,:,:)/1000+1)./(aliquot_deltas_pisCorr(:,7,:,:)./1000+1)-1)*1000]; % predictor variables = dN2/Ar AND dO2Ar (= [q_o2n2/q_arn2 -1]*1000)
[CS_4036Ar,csStats_36Ar] = calculateChemSlope(x_temp,aliquot_deltas_pisCorr(:,4,:,:),aliquot_metadata,iCS_AddN2 | iCS_AddO2,true);
[CS_4038Ar,csStats_38Ar] = calculateChemSlope(x_temp,aliquot_deltas_pisCorr(:,5,:,:),aliquot_metadata,iCS_AddN2 | iCS_AddO2,true);


% Make CS Corrections
csValues = [{CS_15N} {CS_18O} {CS_17O} {CS_ArN2} {CS_4036Ar} {CS_4038Ar}];
csRegressors = [
    {aliquot_deltas_pisCorr(:,delta_names=='d15N',:,:);}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='d18O',:,:)}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='d17O',:,:)}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:);}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='d4036Ar',:,:)}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='d4038Ar',:,:);}, ...
    ];
csPredictors = [
    {aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:);}, ...
    {((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1).^-1-1)*1000}, ...
    {((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1).^-1-1)*1000}, ...
    {aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:);}, ...
    {[((aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1).^-1-1)*1000 ((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1)./(aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1)-1)*1000]}, ...
    {[((aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1).^-1-1)*1000 ((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1)./(aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1)-1)*1000]}, ...
    ];

csCorr = cell(size(csValues)); CS = cell(size(csValues));
for ii = 1:length(csValues)
    [csCorr{ii},CS{ii}] = makeCsCorr(csRegressors{ii},aliquot_metadata.msDatetime,csPredictors{ii},csValues{ii});
end

aliquot_deltas_pisCorr_csCorr = aliquot_deltas_pisCorr;
aliquot_deltas_pisCorr_csCorr(:,[1 2 3 7 4 5],:,:) = [csCorr{:}];

%% Plot all the Chem Slope Experiments
% Plot all the different chem slope experiments for all the different chem
% slope effects

% Make Variables for Plotting
csStats = [csStats_15N csStats_18O csStats_17O csStats_ArN2];
csXLabels = {delta_labels(delta_names=='dO2N2') ['\deltaN_2/O_2 [' char(8240) ']'] ['\deltaN_2/O_2 [' char(8240) ']'] delta_labels(delta_names=='dO2N2')};
csYLabels = {delta_labels(delta_names=='d15N') delta_labels(delta_names=='d18O') delta_labels(delta_names=='d17O') delta_labels(delta_names=='dArN2')};

% Plot the Univariate Chem Slopes
for ii = 1:length(csStats)
    figure
    idx_csExperiments = find(~isnat(csStats(ii).csDatetime));
    for jj=1:length(idx_csExperiments)
        subplot(1,length(idx_csExperiments),jj); hold on
        plot(csStats(ii).xData{idx_csExperiments(jj)},csStats(ii).yData{idx_csExperiments(jj)},'xk') % Plot the individual aliquots
        plot(csStats(ii).xData{idx_csExperiments(jj)},polyval([csStats(ii).slope(idx_csExperiments(jj),:) csStats(ii).intercept(idx_csExperiments(jj))],csStats(ii).xData{idx_csExperiments(jj)}),'-r') % Plot the fitted line
        text(max(csStats(ii).xData{idx_csExperiments(jj)}),min(csStats(ii).yData{idx_csExperiments(jj)}),compose('CS = %.2f per meg/per mil\nr^2 = %.4f',csStats(ii).slope(idx_csExperiments(jj),:)*1000,csStats(ii).corrcoef(idx_csExperiments(jj)).^2),'HorizontalAlignment','Right','VerticalAlignment','bottom')
        
        xlabel(csXLabels{ii});
        ylabel(csYLabels{ii});
        
        title(datestr(csStats(ii).csDatetime(idx_csExperiments(jj)),'yyyy-mmm-dd'))
    end
    ax = get(gcf,'Children');
    set(ax,'XLim',[min([ax.XLim]) max([ax.XLim])]);
    set(ax,'YLim',[min([ax.YLim]) max([ax.YLim])]);
    suptitle([csYLabels{ii} ' Chem Slopes'])
end

% Plot the Bivariate Chem Slopes
csStats = [csStats_36Ar csStats_38Ar];
csXLabels = {['\deltaN_2/Ar [' char(8240) ']'] ['\deltaN_2/Ar [' char(8240) ']']};
csYLabels = {['\deltaO_2/Ar [' char(8240) ']'] ['\deltaO_2/Ar [' char(8240) ']']};
csZLabels = {delta_labels(delta_names=='d4036Ar') delta_labels(delta_names=='d4038Ar')};

for ii = 1:length(csStats)
    figure
    idx_csExperiments = find(~isnat(csStats(ii).csDatetime));
    for jj=1:length(idx_csExperiments)
        subplot(1,length(idx_csExperiments),jj); hold on
        plot3(csStats(ii).xData{idx_csExperiments(jj)}(:,1),csStats(ii).xData{idx_csExperiments(jj)}(:,2),csStats(ii).yData{idx_csExperiments(jj)},'xk'); % Plot the individual aliquots
        [X1,X2]=meshgrid(linspace(min(csStats(ii).xData{idx_csExperiments(jj)}(:,1)),max(csStats(ii).xData{idx_csExperiments(jj)}(:,1)),20),linspace(min(csStats(ii).xData{idx_csExperiments(jj)}(:,2)),max(csStats(ii).xData{idx_csExperiments(jj)}(:,2)),20)); % Create a regularly spaced grid
        surf(X1,X2,X1.*csStats(ii).slope(idx_csExperiments(jj),1) + X2*csStats(ii).slope(idx_csExperiments(jj),2) + csStats(ii).intercept(idx_csExperiments(jj)));  % Plot the fitted surface on the grid
        
        text(max(csStats(ii).xData{idx_csExperiments(jj)}(:,1)),min(csStats(ii).xData{idx_csExperiments(jj)}(:,2)),min(csStats(ii).yData{idx_csExperiments(jj)}),compose('CS N_2/Ar = %.2f per meg/per mil\nCS O_2/Ar = %.2f per meg/per mil',csStats(ii).slope(idx_csExperiments(ii),1)*1000,csStats(ii).slope(idx_csExperiments(ii),2)*1000),'HorizontalAlignment','Right','VerticalAlignment','bottom');
        
        xlabel(csXLabels(ii))
        ylabel(csYLabels(ii))
        zlabel(csZLabels(ii))
        title(datestr(csStats(ii).csDatetime(idx_csExperiments(jj)),'yyyy-mmm-dd'))
        
        colormap(cbrewer('seq','Greens',20));
        view([40 45]);
    end
    pos=get(gca,'Position');
    colorbar;
    set(gca,'Position',pos);
    
    ax = get(gcf,'Children');
    set(ax(2:end),'XLim',[min([ax(2:end).XLim]) max([ax(2:end).XLim])]);
    set(ax(2:end),'YLim',[min([ax(2:end).YLim]) max([ax(2:end).YLim])]);
    set(ax(2:end),'ZLim',[min([ax(2:end).ZLim]) max([ax(2:end).ZLim])]);
    set(ax(2:end),'CLim',[min([ax(2:end).CLim]) max([ax(2:end).CLim])]);
    
    set(ax(1),'Limits',[min([ax(2:end).CLim]) max([ax(2:end).CLim])]);
    
    suptitle([csZLabels{ii} ' Chem Slopes'])
end


%% Plot the Time-Series of the CS Used for Each Aliquot

csToPlot = [CS{:}];
csNames = {'\delta^{15}N CS','\delta^{18}O CS','\delta^{17}O CS','\deltaAr/N_2 CS','\delta^{40}/_{36}Ar CS','\delta^{40}/_{36}Ar CS','\delta^{40}/_{38}Ar CS','\delta^{40}/_{38}Ar CS'};
slopes_36Ar = [csStats_36Ar.slope]; slopes_38Ar = [csStats_38Ar.slope];
csDatetimes = {[csStats_15N.csDatetime] [csStats_18O.csDatetime] [csStats_17O.csDatetime] [csStats_ArN2.csDatetime] [csStats_36Ar.csDatetime] [csStats_36Ar.csDatetime] [csStats_38Ar.csDatetime] [csStats_38Ar.csDatetime]};
csSlopes = {[csStats_15N.slope] [csStats_18O.slope] [csStats_17O.slope] [csStats_ArN2.slope] slopes_36Ar(:,1) slopes_36Ar(:,2) slopes_38Ar(:,1) slopes_38Ar(:,2)};

stackedFig(numel(csNames))
for ii=1:numel(csNames)
    stackedFigAx(ii)
    plot(datenum([csDatetimes{ii}]),[csSlopes{ii}],'ok','MarkerFaceColor',lineCol(ii))
    plot(aliquot_metadata.msDatenum(:,1,1),csToPlot(:,ii),'.','Color',lineCol(ii)*0.5)
    ylabel(csNames{ii})
end

% Add Date of Filament Changes/Refocusing
stackedFigAx

iPlot = massSpecEvents.Event == "New Filament" | massSpecEvents.Event == "Refocus";
plot(datenum([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)']),repmat(get(gca,'YLim')',1,sum(iPlot)),'-r');
text(datenum(massSpecEvents.EndDate(iPlot)),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color','r')

% Set Labels & Limits etc.
title('CS Values used for Correction')
xlabel('Date')
xlim(datenum(["01-Jan-2016" "01-Jan-2019"]))
datetick('x','keeplimits')
stackedFigReset

%% Calculate the LJA Normalization Values
% Identified the aliquots that correspond to the different batches of LJA
% measurements. Plots the distribution (box plots) of aliquot means for
% each batch and identifies and removes any outliers. Calculates the mean
% of aliquot means for each batch to use for LJA normalization.

% Identify the LJA aliquots
iLja = contains(aliquot_metadata.ID1(:,1,1),'LJA');

[calcLja,ljaStats] = calculateLjaValues(aliquot_deltas,aliquot_metadata,iLja);

% Plot the distribution of all LJA aliquot means
figure
for ii = 1:numel(delta_names)
    subplot(1,numel(delta_names),ii)
    histogram(mean(mean(aliquot_deltas_pisCorr_csCorr(iLja,ii,:,:),4),3))
    axis('square');
    xlabel('\delta [per mil]'); ylabel('Counts')
    title(delta_labels(ii))
end

% Make Box Plots


for ii = 1:numel(delta_names)
    allLjaAliquots = []; allLjaAliquotsGrp = [];
    numAliquots = []; stdevAliquots = []; labels = strings;
    for jj = 1:length(ljaStats)
       allLjaAliquots =  [allLjaAliquots; ljaStats(jj).aliquots{ii}];
       allLjaAliquotsGrp = [allLjaAliquotsGrp; repmat(jj,size(ljaStats(jj).aliquots{ii}))];
       numAliquots(jj,:) = ljaStats(jj).N(ii);
       stdevAliquots(jj,:) = ljaStats(jj).stdevs(ii);
       labels(jj) = string(datestr(ljaStats(jj).datetime(ii),'yyyy-mmm-dd'));
    end
    figure; hold on;
     hBoxPlot=boxplot(allLjaAliquots,allLjaAliquotsGrp,'Labels',labels,'Notch','off');
    %hBoxPlot=boxplot(nanmean(mean(aliquot_deltas_pisCorr_csCorr(iLja,ii,:,:),4),3),idxLjaAliquots(iLja),'Labels',labels,'Notch','off');
     plot(allLjaAliquotsGrp,allLjaAliquots,'.k');
    %plot(idxLjaAliquots(iLja),nanmean(mean(aliquot_deltas_pisCorr_csCorr(iLja,ii,:,:),4),3),'.k');
     text(1:max(allLjaAliquotsGrp),repmat(min(allLjaAliquots)*1.05,1,length(ljaStats)),compose(['N = %d\n \\sigma = %.3f' char(8240) '\nSEM = %.4f' char(8240)],numAliquots,stdevAliquots,stdevAliquots./numAliquots),'HorizontalAlignment','center');
    
    ylim('auto');
    ylabel('\delta_{LJA} [per mil]')
    title(['LJA ' delta_labels(ii)])
end


%% Make the LJA Correction

LJA = calcLja;
LJA = fillmissing(LJA,'previous',1);

aliquot_deltas_pisCorr_csCorr_ljaCorr = ((aliquot_deltas_pisCorr_csCorr/1000+1)./(LJA/1000+1)-1)*1000;


%%
set(0,'defaultFigureVisible','on');
disp('>> Script Complete: Turning Figures Back On');