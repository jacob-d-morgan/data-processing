function [cycle_deltas,cycle_metadata] = calcCycles(importedData,isRef)
% CALCCYCLES converts imported ISODAT data to cycle deltas and metadata
%   [CYCLEDELTAS,CYCLEMETADATA] = CALCCYCLES(IMPORTEDDATA,ISREF) calculates
%   delta values from the table of imported data. 
%   ISREF is the variable that indicates whether a given line is a sample
%   (false) or reference (true) measurement. It can be either (1) a logical
%   variable with one column and the same number of rows as importedData,
%   (2) the name of a variable in importedData, (3) the index of a variable
%   in importedData, or (4) a logical index with one row and the same
%   number of columns as imported data.
%   CYCLEMETADATA is the metadata from the sample rows.



end