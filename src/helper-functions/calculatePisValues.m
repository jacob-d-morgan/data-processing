function [pisValues,pisStats] = calculatePisValues(aliquot_deltas,aliquot_metadata,aliquot_deltas_pis,aliquot_metadata_pis)
% CALCULATEPISVALUES Calculate the Pressure Imbalance Sensitivities
%   Calculates the pressure imbalance sensitivity for each delta value, for
%   each pressure imbalance sensitivity experiment in ALIQUOT_DELTAS_PIS.
%   The PIS is calculated by fitting a straight line to the delta values as
%   a function of their pressure imbalance,for each aliquot run as a PIS
%   experiment.
% 
%   CALCULATEPISVALUES(ALIQUOT_DELTAS,ALIQUOT_METADATA,ALIQUOT_DELTAS_PIS,ALIQUOT_METADATA_PIS)
%   returns the PIS for each delta value, for each PIS experiment. The
%   delta values and pressure imbalance of the PIS blocks are passed to the
%   function in the variables ALIQUOT_DELTAS_PIS and ALIQUOT_METADATA_PIS,
%   and for the non-PIS blocks in ALIQUOT_DELTAS and ALIQUOT_METADATA.
%
%   [PISVALUES,PISSTATS] = CALCULATEPISVALUES(...) also outputs a structure
%   of statistics and other information for each PIS experiment, including:
% 
%   pisDatetime - the datetime of the first block in the PIS experiment
%   pressureImbal - the predictor data used to calculate the PIS
%   deltas - the regressor data used to calculate the PIS
%   slope - the slope of the fit, i.e. the PIS
%   intercept - the intercept of the fit
%   rSq - the r-squared value of the linear fit
%   pVal - the p-value of the correlation coefficient
%   pisImbal - the pressure imbalance of the PIS block
%   rejections - the PIS values recommended for rejection (see below)
%   
%   Rejections:
%   A PIS Value is recommended for rejection if one of two criteria is met:
%       (1) The magnitiude of the largest pressure imbalance is <100 mV
%       (PIS for all delta values recommened for rejection).
%       (2) The PIS values has a deviation from the running median more
%       extreme than the the 1st and 99th percentile (PIS for specific
%       delta value recommended for rejection).
%
%  ------------------------------------------------------------------------

%% Parse Inputs
% Check for right number of in/outputs.

narginchk(4,4);
nargoutchk(0,2)


%% Calculate the PIS for each PIS experiment
% Identifies the aliquots containing PIS experiments and calculates the
% pressure imbalance sensitivity of each delta value, for each experiment.

iPIS = ~isnan(aliquot_deltas_pis(:,:,1,1));

% Check that all delta values were measured in each PIS experiment
if sum(any(iPIS,2)) ~= sum(all(iPIS,2))
    % Warn if some delta values are missing a PIS block for some reason...
    warning('Warning: Some delta values are misisng a PIS block for one or more experiments')
end
iPIS = any(iPIS,2); % ...but calculate a PIS for these experiments anyway


% Calculate the PIS for each experiment, for each of the delta values
calcPis = nan(size(aliquot_deltas,[1 2]));

% Pre-allocate variables to be filled in the loop
stats = struct();
stats.pisDatetime = NaT(size(aliquot_deltas,1),1);
stats.pressureImbal = nan(size(aliquot_metadata.pressureImbal) + [0 1 0]);
stats.deltas = nan(size(aliquot_deltas) + [0 0 1 0]);
stats.slope = nan(size(aliquot_deltas,[1 2]));
stats.intercept = nan(size(aliquot_deltas,[1 2]));
stats.rSq = nan(size(aliquot_deltas,[1 2]));
stats.pVal = nan(size(aliquot_deltas,[1 2]));
stats.pisImbal = nan(size(aliquot_deltas,1),1);

% Loop Through the Indices of the PIS Aliquots
for ii=find(iPIS)'
    
    % Find the Pressure Imbalances for the PIS Experiment
    imbal = cat(2,aliquot_metadata.pressureImbal(ii,:,:),aliquot_metadata_pis.pressureImbal(ii,:,:));
    x_temp = mean(imbal,3)'; % predictor variable = pressure imbalance
    
    % Warn if the Largest Imbalance is NOT the PIS Block
    [pImbal, idx] = max(abs(x_temp));
    if idx ~= 5
        warning(['For the PIS experiment on ' datestr(aliquot_metadata.msDatenum(ii,1,1),'dd-mmm-yyyy HH:MM') ' the block with the largest imbalance is block ' num2str(idx) ', not block 5.'])
    end
    
    % Place Pressure Imbalance Data in Final Output
    stats.pisDatetime(ii) = aliquot_metadata.msDatetime(ii);
    stats.pressureImbal(ii,:,:) = imbal;
    stats.pisImbal(ii) = pImbal * sign(x_temp(idx));
    
    % Loop Through the Delta Values
    for jj=1:size(aliquot_deltas,2)
        
        delta = cat(3,aliquot_deltas(ii,jj,:,:),aliquot_deltas_pis(ii,jj,1,:));
        y_temp = squeeze(mean(delta,4)); % response variable = delta value
        
        m_temp = [x_temp ones(size(x_temp))]\y_temp; % calculate the PIS and Intercept
        [R_corr,pVal] = corrcoef(x_temp,y_temp); % calculate the correlation coefficient for the PIS test
        
        % Place Fit Data in Final Outputs
        calcPis(ii,jj)=m_temp(1);
        
        stats.deltas(ii,jj,:,:) = delta;
        stats.slope(ii,jj) = m_temp(1);
        stats.intercept(ii,jj) = m_temp(2);
        stats.rSq(ii,jj) = R_corr(1,2).^2;
        stats.pVal(ii,jj) = pVal(1,2);
    end
end


%% Identify Anomalous PIS values
% Use a running median to remove anomalous PIS values. This method of
% identifying outliers was chosen after testing several different methods.
% Run pisRejectionMethodTesting.m to compare the different approaches.

%pisRejectionMethodTesting; % Run to test different rejection criteria

% Pre-allocate Logical Variable to Log Rejections
stats.rejections = false(size(stats.slope));

% Reject all PIS values where the P Imbalance is smaller than 100 mV
stats.rejections(abs(stats.pisImbal)<100,:) = true;

% Identify PIS values that differ significantly from the running median for a given delta value
movWindow=49; % Moving median window = 7 weeks

for ii = 1:size(aliquot_deltas,2) % loop through the delta values
    x_temp = aliquot_metadata.msDatenum(iPIS,1,1);
    y_temp = calcPis(iPIS,ii,1,1);

    cen = movmedian(y_temp,movWindow,'omitnan','SamplePoints',x_temp); % Calculate moving median for given window width
    dev = y_temp-cen; % Detrend time-series by calculating deviation of each point from moving median

    [CDF,edges] = histcounts(dev,'BinMethod','fd','Normalization','cdf'); % Calculate CDF of the deviations, using Freedman-Diaconis rule for bin widths

    low = cen + edges(find(CDF>0.01,1,'first')); % Lower Bound = Median - First Percentile Deviation
    upp = cen + edges(find(CDF>0.99,1,'first')); % Upper Bound = Median + Ninety Ninth Percentile Deviation
    iRej = (y_temp > upp) | (y_temp < low);

    stats.rejections(iPIS,ii) = stats.rejections(iPIS,ii) | iRej; % assign the rejections back to the full-size variable

end

%% Assign Outputs
% Assign final output variables.

pisValues = calcPis;
pisStats = stats;
    
end % end pisCorr