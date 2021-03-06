function [pisValues,pisStats] = calculatePisValues(rawDataset,rawPisDataset)
% CALCULATEPISVALUES Calculate the Pressure Imbalance Sensitivities
%   Calculates the pressure imbalance sensitivity for each delta value, for
%   each pressure imbalance sensitivity experiment in RAWPISDATASET. The
%   PIS is calculated by fitting a straight line to the mean delta value of
%   each block as a function of the pressure imbalance of that block, for
%   each aliquot run as a PIS experiment.
% 
%   CALCULATEPISVALUES(RAWDATASET,RAWPISDATASET) returns the PIS for each
%   delta value, for each PIS experiment. The delta values and pressure
%   imbalance are passed to the function in the structure fields 'metadata'
%   and 'deltas' of RAWPISDATASET and RAWDATASET for the PIS and non-PIS
%   blocks respectively.
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

narginchk(2,2);
nargoutchk(0,2)


%% Calculate the PIS for each PIS experiment
% Identifies the aliquots containing PIS experiments and calculates the
% pressure imbalance sensitivity of each delta value, for each experiment.

iPIS = ~isnan(rawPisDataset.deltas{:,:}(:,:,1,1));

% Check that all delta values were measured in each PIS experiment
if sum(any(iPIS,2)) ~= sum(all(iPIS,2))
    % Warn if some delta values are missing a PIS block for some reason...
    warning('Warning: Some delta values are misisng a PIS block for one or more experiments')
end
iPIS = any(iPIS,2); % ...but calculate a PIS for these experiments anyway


% Calculate the PIS for each experiment, for each of the delta values
calcPis = nan(size(rawDataset.deltas{:,:},[1 2]));

% Pre-allocate variables to be filled in the loop
stats = struct();
stats.pisDatetime = NaT(size(rawDataset.deltas{:,:},1),1);
stats.pressureImbal = nan(size(rawDataset.metadata.pressureImbal) + [0 0 1 0]);
stats.deltas = nan(size(rawDataset.deltas{:,:}) + [0 0 1 0]);
stats.slope = nan(size(rawDataset.deltas{:,:},[1 2]));
stats.intercept = nan(size(rawDataset.deltas{:,:},[1 2]));
stats.rSq = nan(size(rawDataset.deltas{:,:},[1 2]));
stats.pVal = nan(size(rawDataset.deltas{:,:},[1 2]));
stats.pisImbal = nan(size(rawDataset.deltas{:,:},1),1);

% Loop Through the Indices of the PIS Aliquots
for ii=find(iPIS)'
    
    % Find the Pressure Imbalances for the PIS Experiment
    imbal = cat(3,rawDataset.metadata.pressureImbal(ii,:,:,:),rawPisDataset.metadata.pressureImbal(ii,:,:,:));
    x_temp = squeeze(mean(imbal,4)); % predictor variable = block mean pressure imbalance
    
    % Warn if the Largest Imbalance is NOT the PIS Block
    [pImbal, idx] = max(abs(x_temp));
    if idx ~= 5
        warning(['For the PIS experiment on ' datestr(rawDataset.metadata.msDatenum(ii,1,1),'dd-mmm-yyyy HH:MM') ' the block with the largest imbalance is block ' num2str(idx) ', not block 5.'])
    end
    
    % Place Pressure Imbalance Data in Final Output
    stats.pisDatetime(ii) = rawDataset.metadata.msDatetime(ii);
    stats.pressureImbal(ii,:,:,:) = imbal;
    stats.pisImbal(ii) = pImbal * sign(x_temp(idx));
    
    % Loop Through the Delta Values
    for jj=1:size(rawDataset.deltas{:,:},2)
        
        delta = cat(3,rawDataset.deltas{:,:}(ii,jj,:,:),rawPisDataset.deltas{:,:}(ii,jj,1,:));
        y_temp = squeeze(mean(delta,4)); % response variable = block mean delta value
        
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

% Reject PIS values with Largest Deviations from Running Median
movWindow=49; % Moving median window = 7 weeks

for ii = 1:size(rawDataset.deltas{:,:},2) % loop through the delta values
    x_temp = rawDataset.metadata.msDatenum(iPIS,1,1);
    y_temp = calcPis(iPIS,ii,1,1);

    % Calculate Deviations from Running Median
    cen = movmedian(y_temp,movWindow,'omitnan','SamplePoints',x_temp);
    dev = y_temp-cen;

    % Calculate CDF of the deviations, using Freedman-Diaconis rule for bin widths
    [CDF,edges] = histcounts(dev,'BinMethod','fd','Normalization','cdf');

    low = cen + edges(find(CDF>0.01,1,'first')); % Lower Bound = Median - 1st Percentile Deviation
    upp = cen + edges(find(CDF>0.99,1,'first')); % Upper Bound = Median + 99th Percentile Deviation
    iRej = (y_temp > upp) | (y_temp < low);

    % Assign the Rejections Back to the Full-size Variable
    stats.rejections(iPIS,ii) = stats.rejections(iPIS,ii) | iRej;

end

%% Assign Outputs
% Assign final output variables.

pisValues = calcPis;
pisStats = stats;
    
end % end calculatePisValues