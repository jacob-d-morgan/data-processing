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