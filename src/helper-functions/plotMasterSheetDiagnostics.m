function plotMasterSheetDiagnostics(masterSheet)
% PLOTMASTERSHEETDIAGNOSTICS plots certain key correction diagnostics from
% a data import and generating a master sheet.


%% Plot a time-series of the PIS and related parameters
% This is useful to identify aliquots where the PIS block did no run
% correctly or where the r-squared is low, suggesting a poor determination
% of the PIS. These cases can be filtered out below.

pisStats = masterSheet.correctionDiagnostics.PIS;
delta_names = string(masterSheet.corrDeltas.Properties.VariableNames);
delta_labels = string(masterSheet.corrDeltas.Properties.VariableDescriptions);
delta_units = string(masterSheet.corrDeltas.Properties.VariableUnits);

% Import Mass Spec Events to Annotate Figures
massSpecEvents = getMassSpecEvents;
massSpecEvents.Event = categorical(massSpecEvents.Event);

iPIS = ~isnat(pisStats.pisDatetime);

stackedFig(3,'RelSize',[0.4 1.7 0.9],'Overlap',[-10 -10]);
xlim(stackedFigAx,datenum([2013 2019],[01 01],[01 01]));

% Plot R-Squared of each PIS Experiment
stackedFigAx(1)
for ii = 1:numel(delta_names)
    set(gca,'ColorOrderIndex',ii)
    plot(datenum(pisStats.pisDatetime(iPIS,:)),pisStats.rSq(iPIS,ii),'-ok','MarkerIndices',find(pisStats.rejections(iPIS,ii)),'MarkerFaceColor',lineCol(ii)); % Plot all the r-squared values with markers for rejected values 
    legH(ii)=plot(datenum(pisStats.pisDatetime(iPIS & ~pisStats.rejections(:,ii),:)),pisStats.rSq(iPIS & ~pisStats.rejections(:,ii),ii),'-s','MarkerFaceColor',lineCol(ii)); % Replot as overlay, omitting rejected aliquots
end
legend(legH,delta_labels,'Orientation','Horizontal','Location','South')
ylabel('r^2');
ylim([0.7 1]);

% Plot PIS Value for each PIS Experiment
stackedFigAx(2)
for ii = 1:numel(delta_names)
    set(gca,'ColorOrderIndex',ii)
    plot(datenum(pisStats.pisDatetime(iPIS,:)),pisStats.slope(iPIS,ii),'-ok','MarkerIndices',find(pisStats.rejections(iPIS,ii)),'MarkerFaceColor',lineCol(ii)); % Plot all the r-squared values with markers for rejected values 
    legH(ii)=plot(datenum(pisStats.pisDatetime(iPIS & ~pisStats.rejections(:,ii),:)),pisStats.slope(iPIS & ~pisStats.rejections(:,ii),ii),'s-','MarkerFaceColor',lineCol(ii)); % Replot as overlay, omitting rejected aliquots
end
ylabel(['PIS [' char(8240) '/' char(8240) ']']);
ylim([-0.01 0.005]);

% Plot Pressure Imbalance for each PIS Experiment
stackedFigAx(3)
plot(datenum(pisStats.pisDatetime(iPIS & any(~pisStats.rejections,2),:)),pisStats.pisImbal(iPIS & any(~pisStats.rejections,2),:),'-^','Color',lineCol(9));
ylabel('Pressure Imbalance [per mil]');
ylim([-600 0])

% Set Labels & Limits etc.
stackedFigAx;
xlim(stackedFigAx,datenum([2013 2019],[01 01],[01 01]));
datetick('x','keeplimits');
drawnow;
stackedFigReset

%% Plot Time Series of PIS Value Used for Each Aliquot
stackedFig(numel(delta_names))
for ii=1:numel(delta_names)
    stackedFigAx(ii)
    plot(datenum(pisStats.pisDatetime(~pisStats.rejections(:,ii),:)),pisStats.slope(~pisStats.rejections(:,ii),ii),'o','Color','none','MarkerFaceColor',lineCol(ii))
    plot(masterSheet.metadata.msDatenum(:,1,1),masterSheet.correctionCoeffs.PIS{:,ii},'.','Color',lineCol(ii)*0.5)
    ylabel([delta_labels{ii} ' [' delta_units{ii} '/mV]'])
end

% Add Date of Filament Changes/Refocusing
stackedFigAx

iPlot = massSpecEvents.Event == "New Filament" | massSpecEvents.Event == "Refocus";
plot(datenum([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)']),repmat(get(gca,'YLim')',1,sum(iPlot)),'-r');
text(datenum(massSpecEvents.EndDate(iPlot)),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color','r')

% Set Labels & Limits etc.
title('PIS Values used for Correction')
xlim(stackedFigAx,datenum([2013 2019],[01 01],[01 01]));
datetick('x','keeplimits')
drawnow;
stackedFigReset


%% Plot all the Chem Slope Experiments
% Plot all the different chem slope experiments for all the different chem
% slope effects.

% Plot the Univariate Chem Slopes
csUnivar = {'d15N','d18O','d17O','dArN2'};

% Make Variables for Plotting
csStats = masterSheet.correctionDiagnostics.CS;
csXLabels = [delta_labels(delta_names=='dO2N2') '\deltaN_2/O_2 ' '\deltaN_2/O_2' delta_labels(delta_names=='dO2N2')];
csYLabels = [delta_labels(delta_names=='d15N') delta_labels(delta_names=='d18O') delta_labels(delta_names=='d17O') delta_labels(delta_names=='dArN2')];
csXLabels = append(csXLabels,[' [' delta_units{delta_names=='dO2N2'} ']']);
csYLabels = append(csYLabels,[' [' delta_units{delta_names=='dO2N2'} ']']);
csColors = [lineCol(1:3); lineCol(7)];

% Loop Through Different Isotope Ratio Chem Slopes 
for ii = 1:length(csUnivar)
    figure
    idx_csExperiments = find(~isnat(csStats.(csUnivar{ii}).csDatetime));
    for jj=1:length(idx_csExperiments)
        subplot(1,length(idx_csExperiments),jj); hold on
        
        iRej = csStats.(csUnivar{ii}).rejections{idx_csExperiments(jj)};
        xToPlot = mean(csStats.(csUnivar{ii}).xData{idx_csExperiments(jj)},[4 3]);
        yToPlot = mean(csStats.(csUnivar{ii}).yData{idx_csExperiments(jj)},[4 3]);
        pToPlot = [csStats.(csUnivar{ii}).slope(idx_csExperiments(jj),:) csStats.(csUnivar{ii}).intercept(idx_csExperiments(jj))];
        strToPlot = compose('CS = %.2f per meg/per mil\nr^2 = %.4f',pToPlot(1)*1000,csStats.(csUnivar{ii}).rSq(idx_csExperiments(jj)));
        
        plot(xToPlot(iRej),yToPlot(iRej),'xk') % Plot the individual aliquots, without rejections
        plot(xToPlot(~iRej),yToPlot(~iRej),'ok','MarkerFaceColor',csColors(ii,:)) % Plot the individual aliquots, without rejections
        plot(xToPlot,polyval(pToPlot,xToPlot),'-','Color',csColors(ii,:)*0.7) % Plot the fitted line
        set(gca,'Children',flipud(get(gca,'Children')));
        
        text(max(xToPlot),min(yToPlot),strToPlot,'HorizontalAlignment','Right','VerticalAlignment','bottom')
        
        xlabel(csXLabels{ii});
        ylabel(csYLabels{ii});
        
        title(datestr(csStats.(csUnivar{ii}).csDatetime(idx_csExperiments(jj)),'yyyy-mmm-dd'))
    end
    ax = get(gcf,'Children');
    set(ax,'XLim',[min([ax.XLim]) max([ax.XLim])]);
    set(ax,'YLim',[min([ax.YLim]) max([ax.YLim])]);
    suptitle([csYLabels{ii} ' Chem Slopes'])
end

% Plot the Bivariate Chem Slopes
csBivar = {'d4036Ar','d4038Ar'};

% Make Variables For Plotting
csXLabels = {['\deltaN_2/Ar [' char(8240) ']'] ['\deltaN_2/Ar [' char(8240) ']']};
csYLabels = {['\deltaO_2/Ar [' char(8240) ']'] ['\deltaO_2/Ar [' char(8240) ']']};
csZLabels = [delta_labels(delta_names=='d4036Ar') delta_labels(delta_names=='d4038Ar')];
csZLabels = append(csZLabels,[' [' delta_units{delta_names=='d4036Ar'} ']']);
csColors = flipud(cat(3,cbrewer('seq','Purples',20),cbrewer('seq','Oranges',20)));

% Loop Through Different Isotope Ratio Chem Slopes
for ii = 1:length(csBivar)
    figure
    idx_csExperiments = find(~isnat(csStats.(csBivar{ii}).csDatetime));
    for jj=1:length(idx_csExperiments)
        subplot(1,length(idx_csExperiments),jj); hold on
        
        xToPlot = mean(csStats.(csBivar{ii}).xData{idx_csExperiments(jj)}(:,1),[4 3]);
        yToPlot = mean(csStats.(csBivar{ii}).xData{idx_csExperiments(jj)}(:,2),[4 3]);
        zToPlot = mean(csStats.(csBivar{ii}).yData{idx_csExperiments(jj)},[4 3]);
        pToPlot = [csStats.(csBivar{ii}).slope(idx_csExperiments(jj),:) csStats.(csBivar{ii}).intercept(idx_csExperiments(jj))];
        strToPlot = compose('CS N_2/Ar = %.2f per meg/per mil\nCS O_2/Ar = %.2f per meg/per mil',pToPlot(1)*1000,pToPlot(2)*1000);
        [X1,X2]=meshgrid(linspace(min(xToPlot),max(xToPlot),20),linspace(min(yToPlot),max(yToPlot),20)); % Create a regularly spaced grid
        
        plot3(xToPlot,yToPlot,zToPlot,'xk'); % Plot the individual aliquots
        surf(X1,X2,X1.*pToPlot(1) + X2.*pToPlot(2) + pToPlot(3));  % Plot the fitted surface on the grid
        
        text(max(xToPlot),min(yToPlot),min(zToPlot),strToPlot,'HorizontalAlignment','Right','VerticalAlignment','bottom');
        
        xlabel(csXLabels(ii))
        ylabel(csYLabels(ii))
        zlabel(csZLabels(ii))
        
        title(datestr(csStats.(csBivar{ii}).csDatetime(idx_csExperiments(jj)),'yyyy-mmm-dd'))
        
        colormap(csColors(:,:,ii));
        view([40 45]);
    end
    pos=get(gca,'Position');
    cb=colorbar;
    cb.Label.String = csZLabels{ii};
    set(gca,'Position',pos);
    
    ax = get(gcf,'Children');
    set(ax(2:end),'XLim',[min([ax(2:end).XLim]) max([ax(2:end).XLim])]);
    set(ax(2:end),'YLim',[min([ax(2:end).YLim]) max([ax(2:end).YLim])]);
    set(ax(2:end),'ZLim',[min([ax(2:end).ZLim]) max([ax(2:end).ZLim])]);
    set(ax(2:end),'CLim',[min([ax(2:end).CLim]) max([ax(2:end).CLim])]);
    
    set(ax(1),'Limits',[min([ax(2:end).CLim]) max([ax(2:end).CLim])]);
    
    suptitle([csZLabels{ii}(1:end-4) ' Chem Slopes'])
end


%% Plot the Time-Series of the CS Used for Each Aliquot

csToPlot = masterSheet.correctionCoeffs.CS{:,:};
csYLabels = {'\delta^{15}N CS','\delta^{18}O CS','\delta^{17}O CS','\delta^{40}/_{36}Ar CS','\delta^{40}/_{36}Ar CS','\delta^{40}/_{38}Ar CS','\delta^{40}/_{38}Ar CS','\deltaAr/N_2 CS'};
csDatetimes = [csStats.d15N.csDatetime csStats.d18O.csDatetime csStats.d17O.csDatetime csStats.d4036Ar.csDatetime csStats.d4036Ar.csDatetime csStats.d4038Ar.csDatetime csStats.d4038Ar.csDatetime csStats.dArN2.csDatetime];
csSlopes = [csStats.d15N.slope csStats.d18O.slope csStats.d17O.slope csStats.d4036Ar.slope csStats.d4038Ar.slope csStats.dArN2.slope];

stackedFig(numel(csYLabels))
for ii=1:numel(csYLabels)
    stackedFigAx(ii)
    plot(datenum(csDatetimes(:,ii)),csSlopes(:,ii),'ok','MarkerFaceColor',lineCol(ii))
    plot(datenum(masterSheet.metadata.msDatetime(:,1,1)),csToPlot(:,ii),'.','Color',lineCol(ii)*0.7)
    ylabel([csYLabels{ii} ' [' char(8240) '/' char(8240) ']'])
end

% Add Date of Filament Changes/Refocusing
stackedFigAx

iPlot = massSpecEvents.Event == "New Filament" | massSpecEvents.Event == "Refocus";
plot(datenum([massSpecEvents.StartDate(iPlot)'; massSpecEvents.EndDate(iPlot)']),repmat(get(gca,'YLim')',1,sum(iPlot)),'-r');
text(datenum(massSpecEvents.EndDate(iPlot)),repmat(max(get(gca,'YLim')),sum(iPlot),1),massSpecEvents.Event(iPlot),'Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','Color','r')

% Set Labels & Limits etc.
title('CS Values used for Correction')
xlim(stackedFigAx,datenum([2013 2019],[01 01],[01 01]));
datetick('x','keeplimits')
stackedFigReset


%% Plot All the LJA Aliquots and Values
% For each delta value, plot a time-series of the aliquots averaged for
% each set of LJA values. Also plot shading, indicating the standard
% deviation of the individual aliquots from their mean and the
% linear least-squares fit to the aliquots.

ljaStats = masterSheet.correctionDiagnostics.LJA;

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
        
        ylabel(['LJA: ' delta_labels{ii}  ' [' delta_units{ii} ']'])
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
        
        ylabel(['Std. Can: ' delta_labels{ii}  ' [' delta_units{ii} ']'])
        
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
    
    xlim(datenum([2013 2019],[01 01],[01 01]))
    datetick('x','mmm-yyyy','keeplimits')
    ylabel(delta_labels{ii});
    title(['LJA: ' delta_labels{ii}]);
    
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
    
    plot(datenum(masterSheet.metadata.msDatetime(:,1,1)),masterSheet.correctionCoeffs.LJA{:,ii},'.','Color',lineCol(ii)*0.5)
    ylabel([delta_labels{ii} ' [' delta_units{ii} ']'])
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
xlim(datenum([2013 2019],[01 01],[01 01]))
datetick('x','mmm-yyyy','keeplimits')
stackedFigReset
