function importedData = csvReadIsodat(filename)
% CSVREADISODAT   Import an .csv file created by reprocessing ISODAT files
%   IMPORTEDDATA = CSVREADISODAT(FILENAME) reads the file(s) in FILENAME
%   and outputs them in the table, IMPORTEDDATA. Files are read using
%   READTABLE, the built-in MATLAB function.


% -------------------------------------------------------------------------

% Check Input
if isstring(filename)
    filesToRead = filename;
elseif ischar(filename)
    filesToRead = string(filename);
elseif iscellstr(filename)
    filesToRead = string(filename);
else
    error('INVALID INPUT: Input variable must be a string array, char array, or cell array of char arrays');
end


% Check Import Options for Each File
opts = cell(1,length(filesToRead));
numVars = nan(1,length(filesToRead));
varNames = cell(1,length(filesToRead));
varTypes = cell(1,length(filesToRead));

for ii = 1:length(filesToRead)
    opts{ii} = detectImportOptions(filesToRead(ii),'ReadVariableNames',true,'Delimiter',',');
    numVars(ii) = length(opts{ii}.VariableNames);
    varNames{ii} = string(opts{ii}.VariableNames);
    varTypes{ii} = string(opts{ii}.VariableTypes);
end


if length(unique(numVars)) > 1
    error(join(['Invalid Input: Input files appear to have a different number of variables. ', ...
        'detectImportOpts suggests that the input files have ', ...
        join(string(numVars(1:end-1)),', ') ', and ' string(numVars(end)) ' variables repectively.' ],''));
end

if numVars(1) ~= 38
    warning(['Expected the inuput files to contain 38 variables but the files contains ' num2str(numVars(1)), ...
        '. Has there been a change in the reprocessing template?']);
end

if ~isequal(varNames{:})
    error('Invalid Input: Input files appear to have different variable names or a different order of variables.')
end

if ~isequal(varTypes{:})
    warning('Input files appear to have different variable types or a different order of variable types. Variable types of the first file will be used.')
end


% Set My Default Import Options (Based off Auto-detected Options)
optsToUse = opts{1};
optsToUse = setvartype(optsToUse,contains(optsToUse.VariableNames,'rIntensity'),'double'); % Import all the intensities as doubles, even if they are all NaN
optsToUse = setvartype(optsToUse,string(optsToUse.VariableTypes)=='char','string'); % Import all text as string arrays
optsToUse = setvartype(optsToUse,'TimeCode','datetime');
optsToUse = setvaropts(optsToUse,'TimeCode','InputFormat','yyyy/MM/dd HH:mm:ss','DatetimeFormat','yyyy-MM-dd HH:mm:ss');
optsToUse = setvartype(optsToUse,'Date','datetime');
optsToUse = setvaropts(optsToUse,'Date','InputFormat','MM/dd/yy','DatetimeFormat','yyyy-MM-dd HH:mm:ss');
optsToUse = setvartype(optsToUse,{'Identifier1','Identifier2'},'string');

optsToUse.VariableNames(string(optsToUse.VariableNames)=='Row')={'SequenceRow'}; % Otherwise this inteferes with table.Row, which is already a built-in MATLAB command
optsToUse.VariableNames(string(optsToUse.VariableNames)=='TimeCode')={'datetime'};

% Read In Files
importedData = table;
for ii = 1:length(filesToRead)
    disp(compose('Reading File %i: Working...',ii))
    readFile = readtable(filesToRead(ii),optsToUse);
    importedData = [importedData; readFile];
    
    missingIsRef = isnan(importedData.IsRef__);
    if any((missingIsRef))
        importedData(isnan(importedData.IsRef__),:) = [];
        warning(['Missing IsRef Data: Some integrations in file ' num2str(ii) ' are not recorded as either Sample or Reference. Omitting these ' num2str(sum(missingIsRef)) ' rows.'])
    end
    
    disp(compose('Reading File %i: Complete',ii))
end

importedData = sortrows(importedData,'datetime');
importedData.datenum = datenum(importedData.datetime);
importedData.IsRef__ = logical(importedData.IsRef__);

end