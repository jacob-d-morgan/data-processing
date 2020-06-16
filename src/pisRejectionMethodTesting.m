%% PIS Outlier Rejection Method Testing
% Tests and plots several different methods for rejecting 'bad' PIS values.
% I always reject any PIS experiment where the P Imbalance is less than
% 100 mV but I need a more specific criteria to reject 'bad' PIS values for
% each delta value.

%% Significance of Correlation Coefficient
% Reject those with a p-value greater than 0.05.
%
% This method does poorly because the r-squared value depends mostly on the
% scatter of the four non-PIS blocks (see GH #16 for details).

iNotSig = calcPisPval > 0.05;

stackedFig(numel(delta_cols))
for ii=1:numel(delta_cols)
    stackedFigAx(ii)
    x = aliquot_metadata.msDatenum(:,1,1);
    y = calcPis(:,ii,1,1);
    
    plot(x(~iNotSig(:,ii,1,1)),y(~iNotSig(:,ii,1,1)),'o','Color','none','MarkerFaceColor',lineCol(ii))
    plot(x(iNotSig(:,ii,1,1)),y(iNotSig(:,ii,1,1)),'xk')
    ylabel(delta_cols{ii})
end

stackedFigAx;
xlim(datenum(["01-Jan-2016" "01-Jan-2019"]));
datetick('x','dd-mm-yyyy','keeplimits');
title('Reject PIS Values with p < 0.05');
stackedFigReset;

numRejections = sum(iNotSig);

% calcPis(iNotSig) = nan;
% calcPisRsq(iNotSig) = nan;
% calcPisImbal(iNotSig(:,1,1,1)) = nan;

%% Moving Median +/- 3 Moving MedAD
% Reject PIS values that are more than three median absolute deviations
% away from a moving median of the time series. Test various values for the
% width of the moving median window.
%
% This method seems to do okay but the distance of the rejection threshold
% from the moving median varies wildly due to the irregular spacing of the
% PIS experiments. This leads to both false-positive rejections in places. 
% It is mostly insensitive to window size for windows larger than ~6 weeks.

movWindow = 7:7:7*4*4; % Define a range of window sizes from 1 to 16 weeks
numRej = nan(length(movWindow),numel(delta_cols));

% for ii =1:length(movWindow)
for ii = find(movWindow==42) % Works best for the Moving MedAD method
    stackedFig(numel(delta_cols));
    for jj = 1:numel(delta_cols)
        stackedFigAx(jj)
        x = aliquot_metadata.msDatenum(~isnan(calcPisImbal),1,1);
        y = calcPis(~isnan(calcPisImbal),jj,1,1);
        
        [iRej,low,upp,cen] = isoutlier(y,'movmed',movWindow(ii),'SamplePoints',x); % Reject points three moving MedAD away from the moving median

        H=shadedErrorBar(x,cen,[upp-cen cen-low],'-k'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3; H.mainLine.LineWidth = 1; % Plot the rejection criteria
        plot(x(~iRej),y(~iRej),'o','Color','none','MarkerFaceColor',lineCol(jj)) % Plot the non-rejected data
        errorbar(x(iRej),y(iRej),[],[],repmat(movWindow(ii)/2,sum(iRej),1),repmat(movWindow(ii)/2,sum(iRej),1),'xk','LineWidth',2) % Plot the rejected data and indicate the window size with x error bars
        ylabel(delta_cols{jj})
        
        numRej(ii,jj) = sum(iRej(:,:,1,1)); % tally the number of rejections for each combination of delta value and window size
    end
    
    stackedFigAx;
    xlim(datenum(["01-Jan-2016" "01-Jan-2019"]));
    datetick('x','dd-mmm','KeepLimits')
    title(['Mov. Med. ' char(177) ' 3 Mov. MedADN, Window = ' num2str(movWindow(ii)) ' days'])
    stackedFigReset;

end

% Plot the sensitivity of number of rejections to window size
figure
plot(movWindow,numRej)
xlabel('Window Size [days]');
ylabel('Number of Rejections')
legend(delta_cols)


%% Moving Median +/- 3 MedAD of the Entire Time-Series
% Reject PIS values that are more than three median absolute deviations
% away from a moving median of the time series. In this case, the median
% deviation is that of the entire dataset, not just the local moving MedAD.
%
% This mehtod does okay in general. However, the MedAD of the entire time
% series is not a very relevant parameter to use for rejection. The scatter
% in each time series depends mostly on the changes associated with 
% re-tuning rather than scatter between experiments. This means that the
% distance of the rejections threshold is often too large, resulting in
% false negatives in places. I should use a more relevant metric by
% detrending the time-series. 
% It is extremely insensitive to window size for almost all window sizes.

movWindow = 7:7:7*4*4; % Define a range of window sizes from 1 to 16 weeks
numRej = nan(length(movWindow),numel(delta_cols));

% for ii =1:length(movWindow)
for ii = find(movWindow==49) % Works best for the entire time-series method
    stackedFig(numel(delta_cols));
    for jj = 1:numel(delta_cols)
        stackedFigAx(jj)
        x = aliquot_metadata.msDatenum(~isnan(calcPisImbal),1,1);
        y = calcPis(~isnan(calcPisImbal),jj,1,1);
        
        cen = movmedian(y,movWindow(ii),'omitnan','SamplePoints',x);
        MedAD = mad(y,1); % Calculate MedAD of entire dataset (flag = 1 calculates Median AD c.f. Mean AD)
        low = cen - 3*MedAD*-1/(sqrt(2)*erfcinv(3/2)); % Scaling factor used by isoutlier() makes MAD ~ Std dev (for normally dist data, which these are not).
        upp = cen + 3*MedAD*-1/(sqrt(2)*erfcinv(3/2));
        iRej = (y > upp) | (y < low);
        
        H=shadedErrorBar(x,cen,[upp-cen cen-low],'-k'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3; H.mainLine.LineWidth = 1; % Plot the rejection criteria
        plot(x(~iRej),y(~iRej),'o','Color','none','MarkerFaceColor',lineCol(jj)) % Plot the non-rejected data
        errorbar(x(iRej),y(iRej),[],[],repmat(movWindow(ii)/2,sum(iRej),1),repmat(movWindow(ii)/2,sum(iRej),1),'xk','LineWidth',2) % Plot the rejected data and indicate the window size with x error bars
        ylabel(delta_cols{jj})
        
        numRej(ii,jj) = sum(iRej(:,:,1,1)); % tally the number of rejections for each combination of delta value and window size
    end
    
    stackedFigAx;
    xlim(datenum(["01-Jan-2016" "01-Jan-2019"]));
    datetick('x','dd-mmm','KeepLimits')
    title(['Mov Med. ' char(177) '3 MedADN, Window = ' num2str(movWindow(ii)) ' days'])
    stackedFigReset;

end

% Plot the sensitivity of number of rejections to window size
figure
plot(movWindow,numRej)
xlabel('Window Size [days]');
ylabel('Number of Rejections')
legend(delta_cols)


%% Moving Median +/- 3 MedAD of the Entire Detrended Time-Series
% Removes the points with the most extreme deviations from the running
% median, i.e. the points with deviations outside of the median deviation
% (which ought to be, and is ~0) +/- 3 MedAD.
%
% This method uses a more relevant criterion for rejection. However, it is
% very rejection happy, it probably has too many false-positives. I think
% this is because the distribution of the rejections is weird so the +/- 3
% MedADN ends up rejecting more than 1% of the data. Instead, I could just
% calculate these percentiles directly.
% It is still sensitive to window size, even for windows of ~16 weeks


movWindow = 7:7:7*4*4; % Define a range of window sizes from 1 to 16 weeks
dev = nan(sum(~isnan(calcPisImbal)),numel(delta_cols),length(movWindow));
numRej = nan(length(movWindow),numel(delta_cols));

% for ii =1:length(movWindow)
for ii = movWindow==98 % Works best for the detrended method
    stackedFig(numel(delta_cols));
    for jj = 1:numel(delta_cols)
        stackedFigAx(jj)
        x = aliquot_metadata.msDatenum(~isnan(calcPisImbal),1,1);
        y = calcPis(~isnan(calcPisImbal),jj,1,1);
        
        cen = movmedian(y,movWindow(ii),'omitnan','SamplePoints',x);
        dev(:,jj,ii) = y-cen; % Detrend time-series by calculating deviation of each point from moving median
        MedAD = mad(y-cen,1); % flag = 1 calculates Median AD c.f. Mean AD
        low = cen - 3*MedAD*-1/(sqrt(2)*erfcinv(3/2)); % scaling factor used by isoutlier() makes MAD ~ StdDev for normally dist data.
        upp = cen + 3*MedAD*-1/(sqrt(2)*erfcinv(3/2));
        iRej = (y > upp) | (y < low);
        
        H=shadedErrorBar(x,cen,[upp-cen cen-low],'-k'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3; H.mainLine.LineWidth = 1; % Plot the rejection criteria
        plot(x(~iRej),y(~iRej),'o','Color','none','MarkerFaceColor',lineCol(jj)) % Plot the non-rejected data
        errorbar(x(iRej),y(iRej),[],[],repmat(movWindow(ii)/2,sum(iRej),1),repmat(movWindow(ii)/2,sum(iRej),1),'xk','LineWidth',2) % Plot the rejected data and indicate the window size with x error bars
        ylabel(delta_cols{jj})
        
        numRej(ii,jj) = sum(iRej(:,:,1,1)); % tally the number of rejections for each combination of delta value and window size
                
    end
    
    stackedFigAx;
    xlim(datenum(["01-Jan-2016" "01-Jan-2019"]));
    datetick('x','dd-mmm','KeepLimits')
    title(['Window = ' num2str(movWindow(ii)) ' days'])
    stackedFigReset;
    
end

% Plot the sensitivity of number of rejections to window size
figure
plot(movWindow,numRej)
xlabel('Window Size [days]');
ylabel('Number of Rejections')
legend(delta_cols)


%% Moving Median +/- 45.5% of the PDF of the Detrended Time-Series
% Removes the points with the most extreme deviations from the running
% median, i.e. the points with deviations outside of the median deviation
% (which ought to be, and is ~0) +/- 3 MedAD.
%
% This method seems to do well. It has asymmetrical rejection thresholds
% but I don't see this as a problem. It doesn't seem to have any big
% problem with false positives or false negatives and it is insensitive to
% window size for windows larger than ~5 weeks. It is, however, sensitive
% to the number and width of the bins (especially d40/38Ar, which has two
% large outliers). The Freedman-Diaconis method seems to do the best job as
% it is based on the IQR, which is not very sensitive to outliers, compared
% to other methods that use the std dev. I also anticipate this method
% improving once I add more data.

    
movWindow = 7:7:7*4*4; % Define a range of window sizes from 1 to 16 weeks
dev = nan(sum(~isnan(calcPisImbal)),numel(delta_cols),length(movWindow));
numRej = nan(length(movWindow),numel(delta_cols));
CDF = cell(length(movWindow),numel(delta_cols));
edges = cell(length(movWindow),numel(delta_cols));

% for ii =1:length(movWindow)
for ii = find(movWindow==49) % Works best for the CDF method
    stackedFig(numel(delta_cols));
	stackedFigAx;
    xlim(datenum(["01-Jan-2016" "01-Jan-2019"]));
    for jj = 1:numel(delta_cols)
        stackedFigAx(jj)
        x = aliquot_metadata.msDatenum(~isnan(calcPisImbal),1,1);
        y = calcPis(~isnan(calcPisImbal),jj,1,1);
        
        cen = movmedian(y,movWindow(ii),'omitnan','SamplePoints',x); % Calculate moving median for given window width
        dev(:,jj,ii) = y-cen; % Detrend time-series by calculating deviation of each point from moving median
        
        [CDF{ii,jj},edges{ii,jj}] = histcounts(dev(:,jj,ii),'BinMethod','fd','Normalization','cdf'); % Calculate CDF of the deviations, using Freedman-Diaconis rule for bin widths
        
        low = cen + edges{ii,jj}(find(CDF{ii,jj}>0.01,1,'first')); % Lower Bound = Median - First Percentile Deviation
        upp = cen + edges{ii,jj}(find(CDF{ii,jj}>0.99,1,'first')); % Upper Bound = Median + Ninety Ninth Percentile Deviation
        iRej = (y > upp) | (y < low);
        
        H=shadedErrorBar(x,cen,[upp-cen cen-low],'-k'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3; H.mainLine.LineWidth = 1; % Plot the rejection criteria
        plot(x(~iRej),y(~iRej),'o','Color','none','MarkerFaceColor',lineCol(jj)) % Plot the non-rejected data
        drawnow; ylim('manual')
        errorbar(x(iRej),y(iRej),[],[],repmat(movWindow(ii)/2,sum(iRej),1),repmat(movWindow(ii)/2,sum(iRej),1),'xk','LineWidth',2) % Plot the rejected data and indicate the window size with x error bars
        ylabel(delta_cols{jj})
        
        numRej(ii,jj) = sum(iRej(:,:,1,1)); % tally the number of rejections for each combination of delta value and window size
                
    end
    
    stackedFigAx;
    xlim(datenum(["01-Jan-2016" "01-Jan-2019"]));
    datetick('x','dd-mmm','KeepLimits')
    title(['Window = ' num2str(movWindow(ii)) ' days'])
    stackedFigReset;

end

% Plot the sensitivity of number of rejections to window size
figure
plot(movWindow,numRej)
xlabel('Window Size [days]');
ylabel('Number of Rejections')
legend(delta_cols)
