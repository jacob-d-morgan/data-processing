%% dataImport %%

% Import Data and Fix Variable Types
%xp2018 = readtable('XP-2018(excelExportIntensityJDM).csv');
xp2018.TimeCode = datenum(xp2018.TimeCode);
xp2018.IsRef__ = logical(xp2018.IsRef__);

% Sort by Time Code
idxTimeCode = find(cellfun(@(varName) strcmp('TimeCode',varName),xp2018.Properties.VariableNames)==1);
xp2018 = sortrows(xp2018,idxTimeCode);

% Interpolate ST Voltages and Calculate Delta Values
intSA = xp2018(~xp2018.IsRef__,1:9);
intST = array2table(interp1(find(xp2018.IsRef__), ...
                xp2018{xp2018.IsRef__,1:9}, ...
                find(~xp2018.IsRef__)));
intST.Properties.VariableNames = intSA.Properties.VariableNames;

delta = table();
delta.d15N = ((intSA.rIntensity29./intSA.rIntensity28)./(intST.rIntensity29./intST.rIntensity28) - 1)*1000;
delta.d18O = ((intSA.rIntensity34./intSA.rIntensity32)./(intST.rIntensity34./intST.rIntensity32) - 1)*1000;
delta.d17O = ((intSA.rIntensity33./intSA.rIntensity32)./(intST.rIntensity33./intST.rIntensity32) - 1)*1000;
delta.d4036Ar = ((intSA.rIntensity40./intSA.rIntensity36)./(intST.rIntensity40./intST.rIntensity36) - 1)*1000;
delta.d4038Ar = ((intSA.rIntensity40./intSA.rIntensity38)./(intST.rIntensity40./intST.rIntensity38) - 1)*1000;
delta.dO2N2 = ((intSA.rIntensity32./intSA.rIntensity28)./(intST.rIntensity32./intST.rIntensity28) - 1)*1000;
delta.dArN2 = ((intSA.rIntensity40./intSA.rIntensity28)./(intST.rIntensity40./intST.rIntensity28) - 1)*1000;

delta.pressure = (intSA.rIntensity28./intST.rIntensity28 - 1)*1000;


