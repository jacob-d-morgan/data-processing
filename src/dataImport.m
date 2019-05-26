%% dataImport %%

% Import Data and Fix Variable Types
%xp2018 = readtable('XP-2018(excelExportIntensityJDM).csv');
xp2018.TimeCode = datenum(xp2018.TimeCode);
xp2018.isRef = logical(xp2018.isRef);


% Sort by Time Code
idxTimeCode = find(cellfun(@(varName) strcmp('TimeCode',varName),xp2018.Properties.VariableNames)==1);
xp2018 = sortrows(xp2018,idxTimeCode);


% Split into Sample and Reference Measurements

workingTable = xp2018_int(1:10000,:);

