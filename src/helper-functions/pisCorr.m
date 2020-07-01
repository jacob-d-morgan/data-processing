function [aliquot_deltas_pisCorr] = pisCorr(aliquot_deltas,aliquot_metadata,aliquot_deltas_pis,aliquot_metadata_pis,varargin)
% PISCORR corrects all the delta values in ALIQUOT_DELTAS for the effect of
% pressure imbalance. The correction is determined using the delta values
% from the PIS experiment blocks and their metadata in ALIQUOT_DELTAS_PIS
% and ALIQUOT_METADATA_PIS.
%
%  ------------------------------------------------------------------------

%% Parse Inputs
narginchk(4,6);

flagPlotFigs = false;
if ~isempty(varargin)
    flagPlotFigs = varargin{2};
end

%% Calculate the PIS for each PIS experiment
% Identifies the aliquots containing PIS experiments and calculates the
% pressure imbalance sensitivity of each delta vlaue for each experiment.

iPIS = ~isnan(aliquot_deltas_pis(:,:,1,1));
if sum(any(iPIS,2)) ~= sum(all(iPIS,2)) % Check that all delta values identify each PIS experiment
    warning('Warning: Some delta values are misisng a PIS block for one or more experiments')
    iPIS = any(iPIS,2); % If some delta values are missing a PIS block somehow, calculate the PIS for the other delta values anyway
else
    iPIS = iPIS(:,1); % Otherwise, just make iPis a vector (from a matrix) by taking the first column.
end


% Calculate the PIS for each experiment for each of the delta values
calcPis = nan(size(aliquot_deltas,[1 2]));
calcPisRsq = nan(size(aliquot_deltas,[1 2]));
calcPisPval = nan(size(aliquot_deltas,[1 2]));
calcPisImbal = nan(size(aliquot_deltas,[1 2]));

for ii=find(iPIS)' % find the indices of the PIS aliquots and loop through them
    for jj=1:size(aliquot_deltas,2) % loop all through delta values, skip the first three columns as these are voltages and pressure imbalance
        
        y_temp = [squeeze(nanmean(aliquot_deltas(ii,jj,:,:),4)); squeeze(nanmean(aliquot_deltas_pis(ii,jj,1,:),4))]; % response variable = the looped delta value from the looped aliquot
        x_temp = [nanmean(aliquot_metadata.pressureImbal(ii,:,:),3)'; nanmean(aliquot_metadata_pis.pressureImbal(ii,:,:),3)']; % predictor variable = the pressure imbalance (col 3) from the looped variable
        m_temp = [ones(size(y_temp)) x_temp]\y_temp; % Calculate the PIS and Intercept
        
        [R_corr,pVal] = corrcoef(x_temp,y_temp); % Find the r-squared correlation coefficient for the PIS test
        [pImbal, idx] = max(abs(x_temp)); % Find the block with the max P Imbalance
        
        if idx ~= 5
            warning(['For the PIS experiment on ' datestr(aliquot_metadata.msDatenum(ii,1,1),'dd-mmm-yyyy HH:MM') ' the block with the largest imbalance is block ' num2str(idx) ', not block 5.'])
        end
        
        calcPis(ii,jj)=m_temp(2);
        calcPisRsq(ii,jj) = R_corr(1,2).^2;
        calcPisPval(ii,jj) = pVal(1,2);
        calcPisImbal(ii,jj) = pImbal * sign(x_temp(idx));
    end
end


%% Identify Anomalous PIS values

% Reject all PIS values where the P Imbalance is smaller than 100 mV
iPisRejections = abs(calcPisImbal)<100;
calcPis(iPisRejections) = nan;
calcPisRsq(iPisRejections) = nan;
calcPisPval(iPisRejections) = nan;
calcPisImbal(iPisRejections) = nan;

%pisRejectionMethodTesting; % Run to test different rejection criteria

% Identify PIS values that differ significantly from the running median for a given delta value
numRej = zeros(1,size(aliquot_deltas,2));
movWindow=49; % Moving median window = 7 weeks


for ii = 1:size(aliquot_deltas,2)
    x_temp = aliquot_metadata.msDatenum(iPIS,1,1);
    y_temp = calcPis(iPIS,ii,1,1);

    cen = movmedian(y_temp,movWindow,'omitnan','SamplePoints',x_temp); % Calculate moving median for given window width
    dev = y_temp-cen; % Detrend time-series by calculating deviation of each point from moving median

    [CDF,edges] = histcounts(dev,'BinMethod','fd','Normalization','cdf'); % Calculate CDF of the deviations, using Freedman-Diaconis rule for bin widths

    low = cen + edges(find(CDF>0.01,1,'first')); % Lower Bound = Median - First Percentile Deviation
    upp = cen + edges(find(CDF>0.99,1,'first')); % Upper Bound = Median + Ninety Ninth Percentile Deviation
    iRej = (y_temp > upp) | (y_temp < low);

    iPisRejections(iPIS,ii) = iPisRejections(iPIS,ii) | iRej; % assign the rejections back to the full-size variable
    numRej(ii) = sum(iRej); % tally the number of rejections for each combination of delta value and window size

end


%% Plot a time-series of the PIS and related parameters
% This is useful to identify aliquots where the PIS block did no run
% correctly or where the r-squared is low, suggesting a poor determination
% of the PIS. These cases can be filtered out below.

if flagPlotFigs
    
    stackedFig(3,'RelSize',[0.4 1.7 0.9],'Overlap',[-10 -10]);
    stackedFigAx
    xlim(datenum(['01 Jan 2016'; '31 Dec 2018']))
    
    % Plot R-Squared of each PIS Experiment
    stackedFigAx(1)
    for ii = 1:size(aliquot_deltas,2)
        set(gca,'ColorOrderIndex',ii)
        plot(aliquot_metadata.msDatenum(iPIS,1,1),calcPisRsq(iPIS,ii),'-ok','MarkerIndices',find(iPisRejections(iPIS,ii)),'MarkerFaceColor',lineCol(ii)); % Plot all the r-squared values with markers for rejected values
        legH(ii)=plot(aliquot_metadata.msDatenum(~iPisRejections(:,ii) & iPIS,1,1),calcPisRsq(~iPisRejections(:,ii) & iPIS,ii),'s-','MarkerFaceColor',lineCol(ii)); % Replot as overlay, omitting rejected aliquots
    end
    %legend(legH,delta_cols,'Orientation','Horizontal','Location','South')
    ylabel('r^2');
    ylim([0.7 1]);
    
    % Plot PIS Value for each PIS Experiment
    stackedFigAx(2)
    for ii = 1:size(aliquot_deltas,2)
        set(gca,'ColorOrderIndex',ii)
        plot(aliquot_metadata.msDatenum(iPIS,1,1),calcPis(iPIS,ii),'-ok','MarkerIndices',find(iPisRejections(iPIS,ii)),'MarkerFaceColor',lineCol(ii)); % Plot all the r-squared values with markers for rejected values
        legH(ii)=plot(aliquot_metadata.msDatenum(~iPisRejections(:,ii) & iPIS,1,1),calcPis(~iPisRejections(:,ii) & iPIS,ii),'s-','MarkerFaceColor',lineCol(ii)); % Replot as overlay, omitting rejected aliquots
    end
    ylabel('PIS [per mil/per mil]');
    ylim([-0.01 0.005]);
    
    % Plot Pressure Imbalance for each PIS Experiment
    stackedFigAx(3)
    plot(aliquot_metadata.msDatenum(iPIS & any(~iPisRejections,2),1,1),calcPisImbal(iPIS & any(~iPisRejections,2),:),'-^','Color',lineCol(9));
    text(aliquot_metadata.msDatenum(iPIS & any(~iPisRejections,2),1,1),calcPisImbal(iPIS & any(~iPisRejections,2)),aliquot_metadata_pis.ID1(iPIS & any(~iPisRejections,2),1,1));
    ylabel('Pressure Imbalance [per mil]');
    ylim([-600 0])
    
    stackedFigAx();
    datetick('x');
    xlim(datenum(['01 Jan 2016'; '31 Dec 2018']))
    
    stackedFigReset
    
end

%% Make PIS Correction
PIS = calcPis;
PIS(iPisRejections) = nan;
PIS = repmat(PIS,[1 1 size(aliquot_deltas,[3 4])]);
PIS = fillmissing(PIS,'previous',1);

aliquot_deltas_pisCorr = aliquot_deltas - permute(aliquot_metadata.pressureImbal,[1 4 2 3]).*PIS;

if flagPlotFigs
    stackedFig(size(aliquot_deltas,2))
    for ii=1:size(aliquot_deltas,2)
        stackedFigAx(ii)
        plot(aliquot_metadata.msDatenum(~iPisRejections(:,ii),1,1),calcPis(~iPisRejections(:,ii),ii),'o','Color','none','MarkerFaceColor',lineCol(ii))
        plot(aliquot_metadata.msDatenum(:,1,1),PIS(:,ii,1,1),'.','Color',lineCol(ii)*0.5)
        %ylabel(delta_cols{ii})
    end
    
    stackedFigAx
    title('PIS Values used for Correction')
    xlabel('Date')
    xlim(datenum(["01-Jan-2016" "01-Jan-2019"]))
    datetick('x','keeplimits')
    stackedFigReset
end
end % end pisCorr