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