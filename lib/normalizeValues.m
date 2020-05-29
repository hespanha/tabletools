function [normalizedValues,normalizationFactor]=normalizeValues(raw,rowSum)
% Scale=normalizeValues(raw,rowSum)
% Divide each row of a matrix by a normalization factor so that 
% the rows of the normalized matrix sum to a desired value.
%
% Inputs:
% raw    [nRows,nCols] - Matrix to normalize
% rowSum [nRows,1]     - Desired row sums. When a single scalar is given
%                        all rows should addup to the given value
% Output:
% normalizedValues    [nRows,nCols] - scaled matrix
% normalizationFactor [nRows,1]     - column vector of values by which each row was divided.
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

if size(rowSum,1)==1
    rowSum=repmat(rowSum,size(raw,1),1);
end
normalizationFactor=sum(raw,2)./rowSum;
normalizedValues=raw./repmat(normalizationFactor,1,size(raw,2));
