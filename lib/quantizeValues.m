function [classes,breakpoints]=quantizeValues(values,categories,breakpoints)
% [classes,breakpoints]=quantizeValues(values,categories,breakpoints)
% 
% Quantizes a vector into nC discrete classes. NaN values are placed in the 1st class.
%
% Inputs:
%   values      [nValues,1] - vector of values to be quantized
%   categories  [nC,m]      - categories (each a row m-vector)
%   breakpoints [nC-1,1]    - threshold values specifying the break points. Specifically,
%                               category 1 when value <= breakpoints(1)
%                               category 2 when breakpoints(1) < value <= breakpoints(2)
%                                ...
%                               category nC-1 when breakpoints(nC-2) < value <= breakpoints(nC-1)
%                               category nC when value > breakpoints(nC-1)
%                             Alternative, it can be either of the following strings:
%                               'naturalBreaks' - breakpoints are automatically computed by finding
%                                                 "natural" classes (uses k-means clustering)
%                               'equalBins'     - breakpoints are automatically computed to achieve
%                                                 an equal-sized bins between the minimum and maximum
% Outputs:
%   classes     [nValues,m] - class for each point
%   breakpoints [nC-1,1]    - breakpoints used for classification
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

nC=size(categories,1);
nValues=size(values,1);
goodValues=all(~isnan(values),2);
nValuesNnan=sum(goodValues);

if ischar(breakpoints)
    switch breakpoints
      case {'naturalBreaks'}
        mn=min(values);mx=max(values);
        breakpoints=(mn+(mx-mn)*(1:nC)'/(nC+1));
        % fprintf('Before: ')
        % disp(breakpoints)
        if 1
            s=warning('off','stats:kmeans:EmptyCluster');
            fprintf('quantizeValues: # values=%d, # valuesNnan=%d, # categories=%d\n',nValues,nValuesNnan,nC);
            if nValuesNnan<nC
                fprintf('quantizeValues: not enough values for natural breaks (# values=%d, # valuesNnan=%d, # categories=%d)\n',nValues,nValuesNnan,nC);
                % not enough values
                c=unique(values(goodValues,:));
                c(end+1:nC,1)=inf;
            else
                [idx,c]=kmeans(values,nC,'start',breakpoints,'emptyaction','singleton',...
                               'replicates',1,'MaxIter',500,'OnlinePhase','off');
            end
            warning(s);
        else
            [idx,c]=litekmeans(values',nC);
        end
        c=sort(c);
        breakpoints=(c(1:end-1)+c(2:end))/2;
        % fprintf('After: ')
        % disp(breakpoints)
      case {'equalBins'}
        mn=min(values);mx=max(values);
        breakpoints=mn+(mx-mn)*(1:nC-1)'/nC;
      otherwise 
        error('quantizeValues: unkown method to compute breakpoints ''%s''\n',breakpoints);
    end
end
ids=1+sum(repmat(values,1,length(breakpoints))>repmat(breakpoints',nValues,1),2);
classes=nan(length(ids),size(categories,2));
classes(~isnan(ids),:)=categories(ids(~isnan(ids)),:);

