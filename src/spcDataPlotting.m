%% SPC Data Plotting
% Uses the datasets assembled by dataImport.m and spcDataPull.m to plot
% figures from the SPC XP dataset.

%% Make Gravitational Correction
spc.d18Ograv = spc.d18O - 2*spc.d15N;
spc.d17Ograv = spc.d17O - 1*spc.d15N;
spc.d36Argrav = spc.d36Ar - 4*spc.d15N;
spc.d38Argrav = spc.d38Ar - 2*spc.d15N;
spc.dO2N2grav = spc.dO2N2 - 4*spc.d15N;
spc.dArN2grav = spc.dArN2 - 12*spc.d15N;

spc_replicates.d18Ograv = spc_replicates.d18O - 2*spc_replicates.d15N;
spc_replicates.d17Ograv = spc_replicates.d17O - 1*spc_replicates.d15N;
spc_replicates.d36Argrav = spc_replicates.d36Ar - 4*spc_replicates.d15N;
spc_replicates.d38Argrav = spc_replicates.d38Ar - 2*spc_replicates.d15N;
spc_replicates.dO2N2grav = spc_replicates.dO2N2 - 4*spc_replicates.d15N;
spc_replicates.dArN2grav = spc_replicates.dArN2 - 12*spc_replicates.d15N;

%% Calculate Pair Differences
pairDiffs.d15N = spc_replicates.d15N(:,1)-spc_replicates.d15N(:,2);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),1)-spc_replicates.d15N(isnan(pairDiffs.d15N),3);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),1)-spc_replicates.d15N(isnan(pairDiffs.d15N),4);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),2)-spc_replicates.d15N(isnan(pairDiffs.d15N),3);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),2)-spc_replicates.d15N(isnan(pairDiffs.d15N),4);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),3)-spc_replicates.d15N(isnan(pairDiffs.d15N),4);

pairDiffs.d18O = spc_replicates.d18O(:,1)-spc_replicates.d18O(:,2);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),1)-spc_replicates.d18O(isnan(pairDiffs.d18O),3);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),1)-spc_replicates.d18O(isnan(pairDiffs.d18O),4);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),2)-spc_replicates.d18O(isnan(pairDiffs.d18O),3);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),2)-spc_replicates.d18O(isnan(pairDiffs.d18O),4);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),3)-spc_replicates.d18O(isnan(pairDiffs.d18O),4);

pairDiffs.d17O = spc_replicates.d17O(:,1)-spc_replicates.d17O(:,2);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),1)-spc_replicates.d17O(isnan(pairDiffs.d17O),3);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),1)-spc_replicates.d17O(isnan(pairDiffs.d17O),4);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),2)-spc_replicates.d17O(isnan(pairDiffs.d17O),3);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),2)-spc_replicates.d17O(isnan(pairDiffs.d17O),4);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),3)-spc_replicates.d17O(isnan(pairDiffs.d17O),4);

pairDiffs.d36Ar = spc_replicates.d36Ar(:,1)-spc_replicates.d36Ar(:,2);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),1)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),3);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),1)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),4);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),2)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),3);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),2)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),4);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),3)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),4);

pairDiffs.d38Ar = spc_replicates.d38Ar(:,1)-spc_replicates.d38Ar(:,2);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),1)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),3);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),1)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),4);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),2)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),3);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),2)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),4);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),3)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),4);

pairDiffs.dO2N2 = spc_replicates.dO2N2(:,1)-spc_replicates.dO2N2(:,2);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),1)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),3);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),1)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),4);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),2)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),3);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),2)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),4);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),3)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),4);

pairDiffs.dArN2 = spc_replicates.dArN2(:,1)-spc_replicates.dArN2(:,2);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),1)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),3);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),1)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),4);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),2)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),3);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),2)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),4);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),3)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),4);



