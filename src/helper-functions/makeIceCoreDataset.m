function iceCore = makeIceCoreDataset(masterSheet,iIceCore,bottomDepths)
% MAKEICECOREDATASET compiles sample means from the  master sheet into a
% dataset of replicate measurements from an ice core.
%   ICECORE = MAKEICECOREDATASET(MASTERSHEET,iICECORE) compiles sample
%   aliquots in the master sheet identified by iICECORE into a dataset of
%   replicate measurements, sorted by depth.
%
%   MASTERSHEET must be a table of data generated by makeMasterSheet.
%
%   iICECORE must be a logical vector with the same length as the height of
%   the master sheet tables, which identifies the target ice core.
%
%   BOTTOMDEPTHS must be a numerical array of sample bottom depths.
%
%
% -------------------------------------------------------------------------

%% Parse Inputs
% Ensure that there are the right number of inputs and that they are each
% of the correct class.

narginchk(3,3)

if ~isstruct(masterSheet)
    error('The first input must be a structure, a "Master Sheet" generated by makeMasterSheet')
end

if ~islogical(iIceCore)
    error('The second input must be a logical array.')
end

if ~isnumeric(bottomDepths)
    error('The third input must be a numeric array.')
end


%% Assemble Tables of Data
% Arrange the data in a table using the bottom depth as a grouping variable
% to identify replicates.

% Calculate the Sample Means
aliquotMeans = varfun(@(x) mean(x,[3 4]),(masterSheet.corrDeltas(iIceCore,:)));
aliquotMeans.Properties = masterSheet.corrDeltas.Properties;

% Identify Replicate Samples
iceCore.metadata = table;
[grp,iceCore.metadata.bottomDepth] = findgroups(bottomDepths);

iceCore.metadata.numRepsNoRej = nan(max(grp),1);
for ii_grp = grp'
    iceCore.metadata.numRepsNoRej(ii_grp) = sum(grp==ii_grp);
end    

% Build Tables of Sample Replicates and Means
iceCore.sampleRepsNoRej = table;
iceCore.sampleMeansNoRej = table;
for ii_delta = string(masterSheet.corrDeltas.Properties.VariableNames)
    iceCore.sampleRepsNoRej.(ii_delta) = nan(max(grp),max(iceCore.metadata.numRepsNoRej));
    for jj_grp = grp'   
        temp = aliquotMeans{grp==jj_grp,ii_delta};
        iceCore.sampleRepsNoRej.(ii_delta)(jj_grp,1:length(temp)) = temp';
    end
end

for ii_delta = string(masterSheet.corrDeltas.Properties.VariableNames)
    iceCore.sampleMeansNoRej.(ii_delta) = splitapply(@mean,aliquotMeans.(ii_delta),grp);
end

% Assign Properties
iceCore.sampleRepsNoRej.Properties = masterSheet.corrDeltas.Properties;
iceCore.sampleMeansNoRej.Properties = masterSheet.corrDeltas.Properties;

iceCore.sampleRepsNoRej.Properties.Description = "Delta values for each sample aliquot analyzed on the XP. Delta values are corrected for Pressure Imbalance and Chemical Slope effects and are normalized to LJA.";
iceCore.sampleMeansNoRej.Properties.Description = "Delta values for each set of replicate sample aliquots from the same bottom depth analyzed on the XP. Delta values are corrected for Pressure Imbalance and Chemical Slope effects and are normalized to LJA.";

end %end makeIceCoreDataset