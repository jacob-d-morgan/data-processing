function [aliquot_deltas_pisCorr,varargout] = pisCorr(aliquot_deltas,aliquot_metadata,aliquot_deltas_pis,aliquot_metadata_pis,varargin)
% PISCORR corrects all the delta values in ALIQUOT_DELTAS for the effect of
% pressure imbalance. The correction is determined using the delta values
% from the PIS experiment blocks and their metadata in ALIQUOT_DELTAS_PIS
% and ALIQUOT_METADATA_PIS.
%
% PIS_CORR_DELTAS = pisCorr(ALIQUOT_DELTAS,ALIQUOT_METADATA,ALIQUOT_DELTAS_PIS,ALIQUOT_METADATA_PIS)
% 
% PISCORR can also optionally output the PIS values used to correct
% ALIQUOT_DELTAS and a structure of statistics relating to the fit for each
% PIS experiment.
%
% [...,PIS_VALUES] = pisCorr(...) outputs an array of PIS values used to
% correct the values in ALIQUOT_DELTAS.
%
% [...,PIS_VALUES,PIS_STATS] = pisCorr(...) outputs an array of PIS values 
% used to correct the values in ALIQUOT_DELTAS and the statistics of the
% fit for each PIS experiment, including:
% 
%   pisCalcVal - the calculated PIS value for each PIS experiment
%   pisCalcRsq - the r-squared value for each PIS experiment linear fit
%   pisCalcPval - the p-value of the correlation coefficient for each PIS experiment
%   pisCaclImbal - the pressure imbalance of the PIS block for each PIS experiment
%   
%  ------------------------------------------------------------------------

%% Parse Inputs
% Check for right number of in/outputs.

narginchk(4,6);
nargoutchk(0,3)


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
% Use a running median to remove anomalous PIS values. This method of
% identifying outliers was chosen after testing several different methods.
% Run pisRejectionMethodTesting.m to compare the different approaches.

%pisRejectionMethodTesting; % Run to test different rejection criteria

% Reject all PIS values where the P Imbalance is smaller than 100 mV
iPisRejections = abs(calcPisImbal)<100;
calcPis(iPisRejections) = nan;
calcPisRsq(iPisRejections) = nan;
calcPisPval(iPisRejections) = nan;
calcPisImbal(iPisRejections) = nan;

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

%% Make PIS Correction
% Use the calculated and filtered Pressure Imbalance Sensitivities to make
% the PIS correction.

PIS = calcPis;
PIS(iPisRejections) = nan;
PIS = repmat(PIS,[1 1 size(aliquot_deltas,[3 4])]);
PIS = fillmissing(PIS,'previous',1);

aliquot_deltas_pisCorr = aliquot_deltas - permute(aliquot_metadata.pressureImbal,[1 4 2 3]).*PIS;

%% Parse Outputs
% Assign optional outputs if requested.

varargout = {};

pisStats.measuredPis = calcPis;
pisStats.rejections = iPisRejections;
pisStats.rSq = calcPisRsq;
pisStats.pVal = calcPisPval;
pisStats.pImbal = calcPisImbal;

if nargout==2
    varargout = {PIS};
elseif nargout==3
    varargout = {PIS pisStats};
end
    
end % end pisCorr