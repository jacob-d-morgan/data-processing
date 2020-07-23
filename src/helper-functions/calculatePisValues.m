function [pisValues,pisStats] = calculatePisValues(aliquot_deltas,aliquot_metadata,aliquot_deltas_pis,aliquot_metadata_pis)
% CALCULATEPISVALUES calculates the pressure imbalance sensitivity
%   Calculates the pressure imbalance sensitivity for each delta value and
%   for each pressure imbalance sensitivity experiment in
%   ALIQUOT_DELTAS_PIS. The PIS is calculated by fitting a line through the
%   PIS experiment blocks in ALIQUOT_DELTAS_PIS and the non-PIS blocks in
%   ALIQUOT_DELTAS and their respective pressure imbalances in
%   ALIQUOT_METADATA_PIS and ALIQUOT_METADATA.
% 
%   CALCULATEPISVALUES(ALIQUOT_DELTAS,ALIQUOT_METADATA,ALIQUOT_DELTAS_PIS,ALIQUOT_METADATA_PIS)
%   gives the PIS for each delta value, for each PIS experiment.
%
%   [PISVALUES,PISSTATS] = CALCULATEPISVALUES(...) outputs the statistics of the
%   fit for each PIS experiment, including:
% 
%   measuredPis - the calculated PIS value for each PIS experiment
%   rejections - the PIS values recommended for rejection (see below)
%   rSq - the r-squared value for each PIS experiment linear fit
%   pVal - the p-value of the correlation coefficient for each PIS
%   pImbal - the pressure imbalance of the PIS block in each PIS experiment
%   
%   Rejections:
%   A PIS Value is recommended for rejection if one of two criteria is met:
%       (1) The pressure imbalance is <100 mV (PIS for all delta values
%       recommened for rejection).
%       (2) The PIS value falls far from the running median of PIS values
%       for that delta value (PIS for specific delta value recommended for
%       rejection).
%  ------------------------------------------------------------------------

%% Parse Inputs
% Check for right number of in/outputs.

narginchk(4,4);
nargoutchk(0,2)


%% Calculate the PIS for each PIS experiment
% Identifies the aliquots containing PIS experiments and calculates the
% pressure imbalance sensitivity of each delta vlaue for each experiment.

iPIS = ~isnan(aliquot_deltas_pis(:,:,1,1));
if sum(any(iPIS,2)) ~= sum(all(iPIS,2)) % Check that all delta values identify each PIS experiment
    warning('Warning: Some delta values are misisng a PIS block for one or more experiments')
    iPIS = any(iPIS,2); % If some delta values are missing a PIS block somehow, calculate the PIS for the other delta values anyway
else
    iPIS = iPIS(:,1); % Otherwise, just make iPIS a vector by taking the first column.
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
% Use a running median to remove anomalous PIS values. This method of
% identifying outliers was chosen after testing several different methods.
% Run pisRejectionMethodTesting.m to compare the different approaches.

%pisRejectionMethodTesting; % Run to test different rejection criteria

% Reject all PIS values where the P Imbalance is smaller than 100 mV
iPisRejections = abs(calcPisImbal)<100;

% Identify PIS values that differ significantly from the running median for a given delta value
numRej = zeros(1,size(aliquot_deltas,2));
movWindow=49; % Moving median window = 7 weeks

% Determine rejections for each delta value
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

%% Assign Outputs
% Assign final outputs, including all the calculated PIS values in the
% pisStats variable and the calculated PIS values without the rejected
% values in the pisValues variable.

pisStats.measuredPis = calcPis;
pisStats.rejections = iPisRejections;
pisStats.rSq = calcPisRsq;
pisStats.pVal = calcPisPval;
pisStats.pImbal = calcPisImbal;

pisValues = calcPis;
    
end % end pisCorr