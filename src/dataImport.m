%% dataImport %%

% clearvars;
cluk; clc;

set(0,'defaultFigureVisible','off'); disp('Run dataImport: Turning Figures Off');

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
[aliquot_deltas,metadata,aliquot_deltas_pis,aliquot_metadata_pis] = makeRawDataset(filesToImport,'includePIS',true); % <-- TOGGLE MY COMMENT

aliquot_metadata = metadata.metadata;
delta_names = metadata.delta_names;
delta_labels = metadata.delta_labels;
delta_units = metadata.delta_units;

metadata_fields = string(fieldnames(aliquot_metadata))';


%% Make PIS Correction
% Correct the delta values in aliquot_deltas for the effect of pressure
% imbalance in the bellows.

[aliquot_deltas_pisCorr,PIS,pisStats] = pisCorr(aliquot_deltas,aliquot_metadata,aliquot_deltas_pis,aliquot_metadata_pis);

%% Plot a time-series of the PIS and related parameters
% This is useful to identify aliquots where the PIS block did no run
% correctly or where the r-squared is low, suggesting a poor determination
% of the PIS. These cases can be filtered out below.

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

stackedFigAx();
datetick('x');
xlim(datenum(['01 Jan 2016'; '31 Dec 2018']))

stackedFigReset

stackedFig(numel(delta_names))
for ii=1:numel(delta_names)
    stackedFigAx(ii)
    plot(aliquot_metadata.msDatenum(~pisStats.rejections(:,ii),1,1),pisStats.measuredPis(~pisStats.rejections(:,ii),ii),'o','Color','none','MarkerFaceColor',lineCol(ii))
    plot(aliquot_metadata.msDatenum(:,1,1),PIS(:,ii,1,1),'.','Color',lineCol(ii)*0.5)
    ylabel(delta_names{ii})
end
stackedFigAx
title('PIS Values used for Correction')
xlabel('Date')
xlim(datenum(["01-Jan-2016" "01-Jan-2019"]))
datetick('x','keeplimits')
stackedFigReset


%% Calculate the Chemical Slopes
% Still need to:
%   1) Weed out the questionable blocks/aliquots that mess up some of the
%      chem slope experiments
%   2) Figure out how I'm going to do the CS Correction for Ar isotopes

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

% Plot all the Chem Slope Experiments
% Plot all the different chem slope experiments for all the different chem
% slope effects

% dO2/N2 Effect on d15N
figure
for ii=1:length(csStats_15N)
    subplot(1,length(csStats_15N),ii); hold on
    plot(csStats_15N(ii).xData,csStats_15N(ii).yData,'xk') % Plot the individual aliquots
    plot(csStats_15N(ii).xData,polyval([csStats_15N(ii).slope csStats_15N(ii).intercept],csStats_15N(ii).xData),'-r') % Plot the fitted line
    text(max(csStats_15N(ii).xData),min(csStats_15N(ii).yData),compose('CS = %.2f per meg/per mil\nr^2 = %.4f',csStats_15N(ii).slope*1000,csStats_15N(ii).corrcoef(1,2).^2),'HorizontalAlignment','Right','VerticalAlignment','bottom')
    axis([-10 350 -0.01 0.25]);
    xlabel('\deltaO_2/N_2 [per mil]');
    ylabel('\delta^{15}N [per mil]');
    title(['\delta^{15}N CS: ' datestr(csStats_15N(ii).datetime,'yyyy-mmm-dd')])
end

% dO2/N2 Effect on dAr/N2
figure
for ii=1:length(csStats_ArN2)
    subplot(1,length(csStats_ArN2),ii); hold on
    plot(csStats_ArN2(ii).xData,csStats_ArN2(ii).yData,'xk') % Plot the individual aliquots
    plot(csStats_ArN2(ii).xData,polyval([csStats_ArN2(ii).slope csStats_ArN2(ii).intercept],csStats_ArN2(ii).xData),'-r') % Plot the fitted line
    text(max(csStats_ArN2(ii).xData),min(csStats_ArN2(ii).yData),compose('CS = %.2f per meg/per mil\nr^2 = %.4f',csStats_ArN2(ii).slope*1000,csStats_ArN2(ii).corrcoef(1,2).^2),'HorizontalAlignment','Right','VerticalAlignment','bottom')
    axis([-10 350 -0.1 1.2]);
    xlabel('\deltaO_2/N_2 [per mil]');
    ylabel('\deltaAr/N_2 [per mil]');
    title(['\deltaAr/N_2 CS: ' datestr(csStats_ArN2(ii).datetime,'yyyy-mmm-dd')])
end

% dN2/O2 Effect on d18O
figure
for ii=1:length(csStats_18O)
    subplot(1,length(csStats_18O),ii); hold on
    plot(csStats_18O(ii).xData,csStats_18O(ii).yData,'xk') % Plot the individual aliquots
    plot(csStats_18O(ii).xData,polyval([csStats_18O(ii).slope csStats_18O(ii).intercept],csStats_18O(ii).xData),'-r') % Plot the fitted line
    text(max(csStats_18O(ii).xData),min(csStats_18O(ii).yData),compose('CS = %.2f per meg/per mil\nr^2 = %.4f',csStats_18O(ii).slope*1000,csStats_18O(ii).corrcoef(1,2).^2),'HorizontalAlignment','Right','VerticalAlignment','bottom')
    axis([-10 350 -0.1 0.1]);
    xlabel('\deltaN_2/O_2 [per mil]');
    ylabel('\delta^{18}O [per mil]');
    title(['\delta^{18}O CS: ' datestr(csStats_18O(ii).datetime,'yyyy-mmm-dd')])
end

% dN2/O2 Effect on d17O
figure
for ii=1:length(csStats_17O)
    subplot(1,length(csStats_17O),ii); hold on
    plot(csStats_17O(ii).xData,csStats_17O(ii).yData,'xk') % Plot the individual aliquots
    plot(csStats_17O(ii).xData,polyval([csStats_17O(ii).slope csStats_17O(ii).intercept],csStats_17O(ii).xData),'-r') % Plot the fitted line
    text(max(csStats_17O(ii).xData),min(csStats_17O(ii).yData),compose('CS = %.2f per meg/per mil\nr^2 = %.4f',csStats_17O(ii).slope*1000,csStats_17O(ii).corrcoef(1,2).^2),'HorizontalAlignment','Right','VerticalAlignment','bottom')
    axis([-10 350 -0.1 1]);
    xlabel('\deltaN_2/O_2 [per mil]');
    ylabel('\delta^{17}O [per mil]');
    title(['\delta^{17}O CS: ' datestr(csStats_17O(ii).datetime,'yyyy-mmm-dd')])
end

% dN2/Ar and dO2/Ar Effects on d40/36Ar
figure
for ii=1:length(csStats_36Ar)
    subplot(1,length(csStats_36Ar),ii); hold on
    plot3(csStats_36Ar(ii).xData(:,1),csStats_36Ar(ii).xData(:,2),csStats_36Ar(ii).yData,'xk'); % Plot the individual aliquots
    [X1,X2]=meshgrid(linspace(min(csStats_36Ar(ii).xData(:,1)),max(csStats_36Ar(ii).xData(:,1)),20),linspace(min(csStats_36Ar(ii).xData(:,2)),max(csStats_36Ar(ii).xData(:,2)),20)); % Create a regularly spaced grid
    surf(X1,X2,X1.*csStats_36Ar(ii).slope(1) + X2*csStats_36Ar(ii).slope(2) + csStats_36Ar(ii).intercept);  % Plot the fitted surface on the grid
    
    text(max(csStats_36Ar(ii).xData(:,1)),min(csStats_36Ar(ii).xData(:,2)),min(csStats_36Ar(ii).yData),compose('CS N_2/Ar = %.2f per meg/per mil\nCS O_2/Ar = %.2f per meg/per mil',csStats_36Ar(ii).slope(1)*1000,csStats_36Ar(ii).slope(2)*1000),'HorizontalAlignment','Right','VerticalAlignment','bottom');
    
    xlabel('\deltaN_2/Ar [per mil]')
    ylabel('\deltaO_2/Ar [per mil]')
    zlabel('\delta^{40}/_{36}Ar [per mil]')
    title(['\delta^{40}/_{36}Ar CS: ' datestr(csStats_36Ar(ii).datetime,'yyyy-mmm-dd')])
    
    colormap(cbrewer('seq','Greens',20));
    axis([-20 350 -20 350 -8 1]);
    caxis([-8 0]);
    view([40 45]);
end
pos=get(gca,'Position');
colorbar;
set(gca,'Position',pos);

% dN2/Ar and dO2/Ar Effects on d40/38Ar
figure
for ii=1:length(csStats_38Ar)
    subplot(1,length(csStats_38Ar),ii); hold on
    plot3(csStats_38Ar(ii).xData(:,1),csStats_38Ar(ii).xData(:,2),csStats_38Ar(ii).yData,'xk'); % Plot the individual aliquots
    [X1,X2]=meshgrid(linspace(min(csStats_38Ar(ii).xData(:,1)),max(csStats_38Ar(ii).xData(:,1)),20),linspace(min(csStats_38Ar(ii).xData(:,2)),max(csStats_38Ar(ii).xData(:,2)),20)); % Create a regularly spaced grid
    surf(X1,X2,X1.*csStats_38Ar(ii).slope(1) + X2*csStats_38Ar(ii).slope(2) + csStats_38Ar(ii).intercept);  % Plot the fitted surface on the grid
    
    text(max(csStats_38Ar(ii).xData(:,1)),min(csStats_38Ar(ii).xData(:,2)),min(csStats_38Ar(ii).yData),compose('CS N_2/Ar = %.2f per meg/per mil\nCS O_2/Ar = %.2f per meg/per mil',csStats_38Ar(ii).slope(1)*1000,csStats_38Ar(ii).slope(2)*1000),'HorizontalAlignment','Right','VerticalAlignment','bottom');
    
    xlabel('\deltaN_2/Ar [per mil]')
    ylabel('\deltaO_2/Ar [per mil]')
    zlabel('\delta^{40}/_{38}Ar [per mil]')
    title(['\delta^{40}/_{38}Ar CS: ' datestr(csStats_38Ar(ii).datetime,'yyyy-mmm-dd')])
    
    colormap(cbrewer('seq','Greens',20));
    axis([-50 350 -50 350 -16 1]);
    caxis([-8 0]);
    view([40 45]);
end
pos=get(gca,'Position');
colorbar;
set(gca,'Position',pos);


%% Make the CS Corrections

CS = [CS_15N CS_18O CS_17O zeros(size(CS_15N)) zeros(size(CS_15N)) zeros(size(CS_15N)) CS_ArN2];
CS = fillmissing(CS,'previous',1);

% Predictor Variables for Univariate Chem Slopes
CS_predictors = [aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:) ... % O2N2 CS on d15N
    ((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1).^-1-1)*1000 ... % N2O2 CS on d18O
    ((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1).^-1-1)*1000 ... % N2O2 CS on d17O
    zeros(size(aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:))) ... % d4036Ar CS Below
    zeros(size(aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:))) ... % d4038Ar CS Below
    zeros(size(aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:))) ... % No CS Corr for dO2N2
    aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)]; % O2N2 CS on dArN2

aliquot_deltas_pisCorr_csCorr = aliquot_deltas_pisCorr - CS.*CS_predictors;

% Argon Isotope CS Corrections
CS_4036Ar = fillmissing(CS_4036Ar,'previous',1);
aliquot_deltas_pisCorr_csCorr(:,delta_names=='d4036Ar',:,:) = aliquot_deltas_pisCorr_csCorr(:,delta_names=='d4036Ar',:,:) - (CS_4036Ar(:,1,:,:).*((aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1).^-1-1)*1000) - (CS_4036Ar(:,2,:,:).*((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1)./(aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1)-1)*1000);

CS_4038Ar = fillmissing(CS_4038Ar,'previous',1);
aliquot_deltas_pisCorr_csCorr(:,delta_names=='d4038Ar',:,:) = aliquot_deltas_pisCorr_csCorr(:,delta_names=='d4038Ar',:,:) - (CS_4038Ar(:,1,:,:).*((aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1).^-1-1)*1000) - (CS_4038Ar(:,2,:,:).*((aliquot_deltas_pisCorr(:,delta_names=='dO2N2',:,:)/1000+1)./(aliquot_deltas_pisCorr(:,delta_names=='dArN2',:,:)/1000+1)-1)*1000);


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
    title(delta_names(ii))
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
    title(['LJA ' delta_names(ii)])
end


%% Make the LJA Correction

LJA = calcLja;
LJA = fillmissing(LJA,'previous',1);

aliquot_deltas_pisCorr_csCorr_ljaCorr = ((aliquot_deltas_pisCorr_csCorr/1000+1)./(LJA/1000+1)-1)*1000;


%%
set(0,'defaultFigureVisible','on');
disp('>> Script Complete: Turning Figures Back On');