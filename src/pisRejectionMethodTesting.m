%% PIS Outlier Rejection Method Testing
% Tests and plots several different methods for rejecting 'bad' PIS values.
% I always reject any PIS experiment where the P Imbalance is less than
% 100 mV but I need a more specific criteria to reject 'bad' PIS values for
% each delta value.

%% Significance of Correlation Coefficient
% Reject those with a p-value greater than 0.05
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

movWindow = 7:7:7*4*4;
numRej = nan(length(movWindow),numel(delta_cols));
for ii =1:length(movWindow)
    stackedFig(numel(delta_cols));
    for jj = 1:numel(delta_cols)
        stackedFigAx(jj)
        x = aliquot_metadata.msDatenum(~isnan(calcPisImbal),1,1);
        y = calcPis(~isnan(calcPisImbal),jj,1,1);
        
        [iRej,low,upp,cen] = isoutlier(y,'movmed',movWindow(ii),'SamplePoints',x);

        H=shadedErrorBar(x,cen,[upp-cen cen-low],'-k'); delete(H.edge); H.patch.FaceColor = lineCol(9)*1.3; H.mainLine.LineWidth = 1;
        plot(x(~iRej),y(~iRej),'o','Color','none','MarkerFaceColor',lineCol(jj))
        errorbar(x(iRej),y(iRej),[],[],repmat(movWindow(ii)/2,sum(iRej),1),repmat(movWindow(ii)/2,sum(iRej),1),'xk','LineWidth',2)
        ylabel(delta_cols{jj})
        
        numRej(ii,jj) = sum(iRej(:,:,1,1));
    end
    
    stackedFigAx;
    xlim(datenum(["01-Jan-2016" "01-Jan-2019"]));
    datetick('x','dd-mmm','KeepLimits')
    title(['Window = ' num2str(movWindow(ii)) ' days'])

end

figure
plot(movWindow,numRej)
xlabel('Window Size [days]');
ylabel('Number of Rejections')
legend(delta_cols)