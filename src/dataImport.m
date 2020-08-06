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
    plot(aliquot_metadata.msDatenum(iPIS,1,1),pisStats.slope(iPIS,ii),'-ok','MarkerIndices',find(pisStats.rejections(iPIS,ii)),'MarkerFaceColor',lineCol(ii)); % Plot all the r-squared values with markers for rejected values 
    legH(ii)=plot(aliquot_metadata.msDatenum(~pisStats.rejections(:,ii) & iPIS,1,1),pisStats.slope(~pisStats.rejections(:,ii) & iPIS,ii),'s-','MarkerFaceColor',lineCol(ii)); % Replot as overlay, omitting rejected aliquots
end
ylabel('PIS [per mil/per mil]');
ylim([-0.01 0.005]);

% Plot Pressure Imbalance for each PIS Experiment
stackedFigAx(3)
plot(aliquot_metadata.msDatenum(iPIS & any(~pisStats.rejections,2),1,1),pisStats.pisImbal(iPIS & any(~pisStats.rejections,2),:),'-^','Color',lineCol(9));
text(aliquot_metadata.msDatenum(iPIS & any(~pisStats.rejections,2),1,1),pisStats.pisImbal(iPIS & any(~pisStats.rejections,2)),aliquot_metadata_pis.ID1(iPIS & any(~pisStats.rejections,2),1,1));
ylabel('Pressure Imbalance [per mil]');
ylim([-600 0])

% Set Labels & Limits etc.
stackedFigAx();
xlim(datenum(['01 Jan 2016'; '31 Dec 2018']))
datetick('x','keeplimits');
drawnow;
stackedFigReset

%% Plot Time Series of PIS Value Used for Each Aliquot
stackedFig(numel(delta_names))
for ii=1:numel(delta_names)
    stackedFigAx(ii)
    plot(aliquot_metadata.msDatenum(~pisStats.rejections(:,ii),1,1),pisStats.slope(~pisStats.rejections(:,ii),ii),'o','Color','none','MarkerFaceColor',lineCol(ii))
    plot(aliquot_metadata.msDatenum(:,1,1),PIS(:,ii),'.','Color',lineCol(ii)*0.5)
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

% == MANUALLY INCLUDE THE ONLY 2016-02-09 REP-0 IN BOTH CS EXPERIMENTS == %
iCS_AddN2(aliquot_metadata.msDatetime(:,1,1)=={'2016-02-08 13:27:59'}) = true;
% ======================================================================= %

[calcCS_15N,csStats_15N] = calculateChemSlope(aliquot_deltas_pisCorr(:,6,:,:),aliquot_deltas_pisCorr(:,1,:,:),aliquot_metadata,iCS_AddO2);
[calcCS_ArN2,csStats_ArN2] = calculateChemSlope(aliquot_deltas_pisCorr(:,6,:,:),aliquot_deltas_pisCorr(:,7,:,:),aliquot_metadata,iCS_AddO2);
[calcCS_18O,csStats_18O] = calculateChemSlope((1/(aliquot_deltas_pisCorr(:,6,:,:)/1000+1)-1)*1000,aliquot_deltas_pisCorr(:,2,:,:),aliquot_metadata,iCS_AddN2);
[calcCS_17O,csStats_17O] = calculateChemSlope((1/(aliquot_deltas_pisCorr(:,6,:,:)/1000+1)-1)*1000,aliquot_deltas_pisCorr(:,3,:,:),aliquot_metadata,iCS_AddN2);

x_temp = [((aliquot_deltas_pisCorr(:,7,:,:)./1000+1).^-1-1)*1000 ((aliquot_deltas_pisCorr(:,6,:,:)/1000+1)./(aliquot_deltas_pisCorr(:,7,:,:)./1000+1)-1)*1000]; % predictor variables = dN2/Ar AND dO2Ar (= [q_o2n2/q_arn2 -1]*1000)
[calcCS_4036Ar,csStats_36Ar] = calculateChemSlope(x_temp,aliquot_deltas_pisCorr(:,4,:,:),aliquot_metadata,iCS_AddN2 | iCS_AddO2,true);
[calcCS_4038Ar,csStats_38Ar] = calculateChemSlope(x_temp,aliquot_deltas_pisCorr(:,5,:,:),aliquot_metadata,iCS_AddN2 | iCS_AddO2,true);


% Make CS Corrections
csValues = [{calcCS_15N} {calcCS_18O} {calcCS_17O} {calcCS_ArN2} {calcCS_4036Ar} {calcCS_4038Ar}];
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
% slope effects.

% Make Variables for Plotting
csStats = [csStats_15N csStats_18O csStats_17O csStats_ArN2];
csXLabels = {delta_labels(delta_names=='dO2N2') ['\deltaN_2/O_2 [' char(8240) ']'] ['\deltaN_2/O_2 [' char(8240) ']'] delta_labels(delta_names=='dO2N2')};
csYLabels = {delta_labels(delta_names=='d15N') delta_labels(delta_names=='d18O') delta_labels(delta_names=='d17O') delta_labels(delta_names=='dArN2')};
csColors = [lineCol(1:3); lineCol(7)];

% Plot the Univariate Chem Slopes
for ii = 1:length(csStats)
    figure
    idx_csExperiments = find(~isnat(csStats(ii).csDatetime));
    for jj=1:length(idx_csExperiments)
        subplot(1,length(idx_csExperiments),jj); hold on
        
        iLjaRej = csStats(ii).rejections{idx_csExperiments(jj)};
        xToPlot = csStats(ii).xData{idx_csExperiments(jj)};
        yToPlot = csStats(ii).yData{idx_csExperiments(jj)};
        pToPlot = [csStats(ii).slope(idx_csExperiments(jj),:) csStats(ii).intercept(idx_csExperiments(jj))];
        strToPlot = compose('CS = %.2f per meg/per mil\nr^2 = %.4f',pToPlot(1)*1000,csStats(ii).rSq(idx_csExperiments(jj)));
        
        plot(xToPlot(iLjaRej),yToPlot(iLjaRej),'xk') % Plot the individual aliquots, without rejections
        plot(xToPlot(~iLjaRej),yToPlot(~iLjaRej),'ok','MarkerFaceColor',csColors(ii,:)) % Plot the individual aliquots, without rejections
        plot(xToPlot,polyval(pToPlot,xToPlot),'-','Color',csColors(ii,:)*0.7) % Plot the fitted line
        set(gca,'Children',flipud(get(gca,'Children')));
        
        text(max(xToPlot),min(yToPlot),strToPlot,'HorizontalAlignment','Right','VerticalAlignment','bottom')
        
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
csColors = flipud(cat(3,cbrewer('seq','Purples',20),cbrewer('seq','Oranges',20)));

for ii = 1:length(csStats)
    figure
    idx_csExperiments = find(~isnat(csStats(ii).csDatetime));
    for jj=1:length(idx_csExperiments)
        subplot(1,length(idx_csExperiments),jj); hold on
        
        xToPlot = csStats(ii).xData{idx_csExperiments(jj)}(:,1);
        yToPlot = csStats(ii).xData{idx_csExperiments(jj)}(:,2);
        zToPlot = csStats(ii).yData{idx_csExperiments(jj)};
        pToPlot = [csStats(ii).slope(idx_csExperiments(jj),:) csStats(ii).intercept(idx_csExperiments(jj))];
        strToPlot = compose('CS N_2/Ar = %.2f per meg/per mil\nCS O_2/Ar = %.2f per meg/per mil',pToPlot(1)*1000,pToPlot(2)*1000);
        [X1,X2]=meshgrid(linspace(min(xToPlot),max(xToPlot),20),linspace(min(yToPlot),max(yToPlot),20)); % Create a regularly spaced grid
        
        plot3(xToPlot,yToPlot,zToPlot,'xk'); % Plot the individual aliquots
        surf(X1,X2,X1.*pToPlot(1) + X2.*pToPlot(2) + pToPlot(3));  % Plot the fitted surface on the grid
        
        text(max(xToPlot),min(yToPlot),min(zToPlot),strToPlot,'HorizontalAlignment','Right','VerticalAlignment','bottom');
        
        xlabel(csXLabels(ii))
        ylabel(csYLabels(ii))
        zlabel(csZLabels(ii))
        
        title(datestr(csStats(ii).csDatetime(idx_csExperiments(jj)),'yyyy-mmm-dd'))
        
        colormap(csColors(:,:,ii));
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
    plot(aliquot_metadata.msDatenum(:,1,1),csToPlot(:,ii),'.','Color',lineCol(ii)*0.7)
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

%% Make the LJA Correction
% Identifies the aliquots that correspond to LJA measurements during
% identical MS conditions (i.e. same filament, std can etc.) and calculates
% the mean of the aliquots, after rejecting outliers. The mean of all the
% aliquots is used to normalize all aliquots measured under identical MS
% conditions, unless there is a trend in the LJA values. In this case, the
% extrapolated values are used.

% Calculate LJA Values
iLja = contains(aliquot_metadata.ID1(:,1,1),'LJA');
[ljaValues,ljaStats] = calculateLjaValues(aliquot_deltas,aliquot_metadata,iLja);

% Make LJA Correction
[deltas_ljaCorr,LJA] = makeLjaCorr(aliquot_deltas_pisCorr_csCorr,aliquot_metadata.msDatetime(:,1,1),ljaStats,ljaValues);


%% Plot All the LJA Aliquots and Values
% For each delta value, plot a time-series of the aliquots averaged for
% each set of LJA values. Also plot shading, indicating the standard
% deviation of the individual aliquots from their mean and the
% linear least-squares fit to the aliquots.

% Loop through each delta value
for ii = 1:numel(delta_names)
    stackedFig(2,'Overlap',30);
    hold on;
    idx_ljaExperiments = find(~isnat(ljaStats.ljaDatetime));
    for jj = 1:length(idx_ljaExperiments)
        % Plot LJA Data:
        stackedFigAx(1)
        iLjaRej = ljaStats.ljaRejections{idx_ljaExperiments(jj)}(:,ii);
        ljaSetTimeRange = [min(ljaStats.ljaAliquotSetDates{idx_ljaExperiments(jj)}) max(ljaStats.ljaAliquotSetDates{idx_ljaExperiments(jj)})];
        ljaP = [ljaStats.ljaSlope(idx_ljaExperiments(jj),ii) ljaStats.ljaIntercept(idx_ljaExperiments(jj),ii)];
        fittedLine = polyval(ljaP,datenum(ljaStats.ljaAliquotSetDates{idx_ljaExperiments(jj)}(~iLjaRej,:)));
        
        % plot mean and std dev of LJA aliquots
        shadedErrorBar(datenum(ljaSetTimeRange),repmat(ljaStats.ljaValues(idx_ljaExperiments(jj),ii),2,1),repmat(std(ljaStats.ljaAliquotSetDeltas{idx_ljaExperiments(jj)}(~iLjaRej,ii)),2,1),{'-','Color',lineCol(jj)})
        % plot fitted line
        plot(datenum(ljaStats.ljaAliquotSetDates{idx_ljaExperiments(jj)}(~iLjaRej,:)),fittedLine,':','color',lineCol(jj)*0.7);
        
        % plot included and rejected LJA aliquots
        plot(datenum(ljaStats.ljaAliquotSetDates{idx_ljaExperiments(jj)}),ljaStats.ljaAliquotSetDeltas{idx_ljaExperiments(jj)}(:,ii),'.','Color',lineCol(jj))
        plot(datenum(ljaStats.ljaAliquotSetDates{idx_ljaExperiments(jj)}(iLjaRej,:)),ljaStats.ljaAliquotSetDeltas{idx_ljaExperiments(jj)}(iLjaRej,ii),'xk')
        
        ylabel(['LJA' delta_labels(ii)])
        yl = ylim;
        
        % Plot Can Data:
        stackedFigAx(2)
        iCanRej = ljaStats.canRejections{idx_ljaExperiments(jj)}(:,ii);
        canSetTimeRange = [min(ljaStats.canAliquotSetDates{idx_ljaExperiments(jj)}) max(ljaStats.canAliquotSetDates{idx_ljaExperiments(jj)})];
        canP = [ljaStats.canSlope(idx_ljaExperiments(jj),ii) ljaStats.canIntercept(idx_ljaExperiments(jj),ii)];
        fittedLine = polyval(canP,datenum(ljaStats.canAliquotSetDates{idx_ljaExperiments(jj)}(~iCanRej,:)));
        
        % plot mean and std dev of aliquots
        shadedErrorBar(datenum(canSetTimeRange),repmat(mean(ljaStats.canAliquotSetDeltas{idx_ljaExperiments(jj)}(~iCanRej,ii)),2,1),repmat(std(ljaStats.canAliquotSetDeltas{idx_ljaExperiments(jj)}(~iCanRej,ii)),2,1),{'-','Color',lineCol(jj)})
        % plot fitted line
        plot(datenum(ljaStats.canAliquotSetDates{idx_ljaExperiments(jj)}(~iCanRej,:)),fittedLine,':','color',lineCol(jj)*0.7);
        
        % plot included and rejected aliquots
        plot(datenum(ljaStats.canAliquotSetDates{idx_ljaExperiments(jj)}),ljaStats.canAliquotSetDeltas{idx_ljaExperiments(jj)}(:,ii),'.','Color',lineCol(jj))
        plot(datenum(ljaStats.canAliquotSetDates{idx_ljaExperiments(jj)}(iCanRej,:)),ljaStats.canAliquotSetDeltas{idx_ljaExperiments(jj)}(iCanRej,ii),'xk')
        
        ylabel(['LJA' delta_labels(ii)])
        
    end
    
    stackedFigAx
    % Add Date of Filament Changes/Refocusing/New Std Cans
    iPlot = massSpecEvents.Event == "New Filament" | ...
        massSpecEvents.Event == "Refocus" | ...
        massSpecEvents.Event == "New Std Cans" | ...
        massSpecEvents.Event == "Swap Std Cans" | ...
        massSpecEvents.Event == "MS Change";
    plot(datenum([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)']),repmat(get(gca,'YLim')',1,sum(iPlot)),'-','Color',lineCol(10));
    text(datenum(massSpecEvents.EndDate(iPlot)),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color',lineCol(10));
    
    xlim(datenum([2016 2019],[01 01],[01 01]))
    datetick('x','mmm-yyyy','keeplimits')
    ylabel(delta_labels{ii});
    title(['LJA: ' delta_names{ii}]);
    
    yls = [ylim(stackedFigAx(1)) ylim(stackedFigAx(2))];
    linkaxes([stackedFigAx(1) stackedFigAx(2)],'y')
    
    stackedFigReset;
    
end


% Plot the Time-Series of the LJA Values Used for Each Aliquot
stackedFig(numel(delta_names))
for ii=1:numel(delta_names)
    stackedFigAx(ii)
    for jj = 1:length(idx_ljaExperiments)
        iLjaRej = ljaStats.ljaRejections{idx_ljaExperiments(jj)}(:,ii);
        ljaSetTimeRange = [min(ljaStats.ljaAliquotSetDates{idx_ljaExperiments(jj)}) max(ljaStats.ljaAliquotSetDates{idx_ljaExperiments(jj)})];
        shadedErrorBar(datenum(ljaSetTimeRange),repmat(ljaStats.ljaValues(idx_ljaExperiments(jj),ii),2,1),repmat(std(ljaStats.ljaAliquotSetDeltas{idx_ljaExperiments(jj)}(~iLjaRej,ii)),2,1),{'-','Color',lineCol(ii)});
    end
    
    plot(datenum(aliquot_metadata.msDatetime(:,1,1)),LJA(:,ii),'.','Color',lineCol(ii)*0.5)
    ylabel(delta_labels{ii})
end

% Add Date of Filament Changes/Refocusing/Can Swaps/MS Changes
stackedFigAx

iPlot = massSpecEvents.Event == "New Filament" | ...
        massSpecEvents.Event == "Refocus" | ...
        massSpecEvents.Event == "New Std Cans" | ...
        massSpecEvents.Event == "Swap Std Cans" | ...
        massSpecEvents.Event == "MS Change";
plot(datenum([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)']),repmat(get(gca,'YLim')',1,sum(iPlot)),'-','Color',lineCol(10));
text(datenum(massSpecEvents.EndDate(iPlot)),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color',lineCol(10))

% Set Labels & Limits etc.
title('LJA Values used for Correction')
xlabel('Date')
xlim(datenum([2016 2019],[01 01],[01 01]))
datetick('x','mmm-yyyy','keeplimits')
stackedFigReset


%%
set(0,'defaultFigureVisible','on');
disp('>> Script Complete: Turning Figures Back On');