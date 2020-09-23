function massSpecEvents = getMassSpecEvents
% GETMASSSPECEVENTS loads a table of data on Mass Spec events
%   Reads data from the first sheet of spreadsheet_metadata.xlsx into a
%   table containing the date of important events that happened to the mass
%   spectrometer. Examples include filament changes, refocussing, or a 
%   change in the standard gas.

warning('off','MATLAB:table:ModifiedAndSavedVarnames');
massSpecEvents = readtable('spreadsheet_metadata.xlsx','Sheet',1);
warning('on','MATLAB:table:ModifiedAndSavedVarnames');