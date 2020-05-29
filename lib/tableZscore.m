function tbl=tableZscore(tbl,ignore)
% tableOut=tableZscore(tableIn)
%
% Returns in tableOut a version of tableIn, where all numeric columns were 
% replaced by their zscores (i.e., normalized to have zero means and unit
% standard deviations -- see help zscore). 
%
% or
%
% tableOut=tableZscore(tableIn,ignore)
%
% Does not change the columns with indices in the vector ''ignore''
%
% NaN values are ignored for the computation of means standard deviation

if nargin<2
    ignore=[];
end

for i=1:size(tbl,2)
    x=tbl{:,i};
    if islogical(x)
        x=double(x);
    end
    if ~ismember(i,ignore) && isnumeric(x)
        %tbl{:,i}=zscore(tbl{:,i});
        tbl{:,i}=(x-repmat(nanmean(x),size(x,1),1))./repmat(nanstd(x),size(x,1),1);
    end
end
