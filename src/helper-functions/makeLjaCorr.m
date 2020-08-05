function [deltasLjaCorr,LJA] = makeLjaCorr(deltasRaw,aliquotDates,ljaValues)
% MAKELJACORR normalizes delta values to measurements of La Jolla Air
%   For an array of delta values, DELTAS_RAW, measured against a standard
%   can, each delta value is normalized to La Jolla Air using measurements
%   of LJA against the same standard can.
%   
%   LJAVALUES should be an M-by-N matrix where M and N are the size of
%   DELTAS_RAW along the first two dimensions. Missing (NaN) values are
%   generally filled with the mean of all LJA measurements against the same
%   standard can, using the same focussing and filament, as logged in
%   spreadsheet_metadata.xlsx. An excpetion to this is if there is a trend
%   in the values of LJA against the standard can through time. If the
%   trend is larger than the scatter in the LJA measurements, LJAVALUES is
%   a linear interpolation/extrapolation of the measured values over the
%   correction period. If no LJA measurements were made during a given
%   correction period, the next available LJA value is used. The LJA values
%   used for each aliquot are given in the second output, LJA.
%
%   e.g. [DELTAS_CSCORR,CSS] = MAKELJACORR(DELTAS_RAW,ALIQUOTDATES,LJAVALUES) 
%
% -------------------------------------------------------------------------