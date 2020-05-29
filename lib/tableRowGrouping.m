function varargout=tableRowGrouping(baseTable,query,groupingVariables,varargin)
% [output variable 1, output variable 2, ...]=tableRowGrouping(...
%        baseTable,...
%        query,
%        {'grouping variable 1','grouping variable 2',... },...
%        {'column a','column b',...}, function 1,...
%        {'column c','column d',...}, function 2,...
%        ... )
%
% First selects a subset of the rows of the table based on a query.
% A query is a function whose input is the table and that returns
% which rows should be keep. E.g., To find the rows whose length
% of the field 'name' is smaller than 20, use:
%   query = @(tbl) length(tbl.name)<20
% See 'help tableSelect'' for further details on the query format.
%
% Then groups the rows of a table based on the values of a given set of
% 'grouping variables', in the sense that each row of the ''grouped''
% table corresponds to a unique value of the 'grouping variables'.
%
% Finally, returns a number of output variables with as many rows as
% the number of rows of the grouped variable, which are obtained by
% providing to given functions the values of specific columns of the
% original table that correspond to each specific row of the grouped
% table.
%
% Closely related to rowfun and varfun, but returns arrays and not a
% table.
%
% Copyright (C) 2013-2014 Stacy & Joao Hespanha

fprintf('  tableRowGrouping: baseTable(%dx%d)... ',size(baseTable));
t0=clock;

if ~isempty(query)
    k=tableSelect(baseTable,query);
    fprintf('  selecting rows in baseTable(%dx%d): ',size(baseTable));
    t1=clock;
    baseTable=baseTable(k,:);
    fprintf('  done baseTable(%dx%d) (%.3f sec)\n',size(baseTable),etime(clock,t1));
end

groupingValues=baseTable{:,groupingVariables};
[values,~,idx]=unique(groupingValues,'rows');
nFuns=(length(varargin)-2)/2;
inputVars=cell(nFuns,1);
funs=cell(nFuns,1);
varargout=cell(nFuns,1);
if mod(length(varargin),2)
    error('tableRowGrouping: number of input arguments must be odd (%d instead)\n',length(varargin)+3);
end
for i=1:2:length(varargin)
    % get input variable indices
    variables=varargin{i}; 
    if ~iscell(variables)
        variables={variables};
    end
    [is,inputVars{(i+1)/2}]=ismember(variables,baseTable.Properties.VariableNames);
    if ~all(is)
        disp('table variable names:');
        disp(baseTable.Properties.VariableNames')
        disp('variable names not found:')
        disp(variables(~is)')
        error('tableRowGrouping: unknown table variable\n');
    end
    % get function
    func=varargin{i+1};
    if ischar(func)
        switch (func)
          case{'first'}
            func=@(x)x(1,:);
          case{'last'}
            func=@(x)x(end,:);
          case {'count'}
            func=@(x)size(x,1);
          case{'sum'}
            func=@(x)sum(x,1);
          case{'mean'}
            func=@(x)mean(x,1);
          case{'median'}
            func=@(x)median(x,1);
          case{'sum_nonan'}
            func=@(x)sum(x(all(~isnan(x),2),:),1);
          case{'mean_nonan'}
            func=@(x)mean(x(all(~isnan(x),2),:),1);
          case{'median_nonan'}
            func=@(x)median(x(all(~isnan(x),2),:),1);
          otherwise,
            func=str2func(func);
        end
    end
    funs{(i+1)/2}=func;
end
for i=1:2:length(varargin)
    func=funs{(i+1)/2};
    testy=func(baseTable{1,inputVars{(i+1)/2}});
    if isnumeric(testy)
        varargout{(i+1)/2}=nan(length(values),0);
    else
        varargout{(i+1)/2}=cell(length(values),0);
    end
end

for j=1:length(values)
    rows=find(idx==j);
    for i=1:length(funs)
        x=baseTable{rows,inputVars{i}};
        y=funs{i}(x);
        varargout{i}(j,1:length(y))=y;
        %fprintf('  %s: %dx%d->%dx%d\n',func2str(funs{i}),size(x),size(y));
    end
end
fprintf('done %d groups (%.3f sec)\n',length(values),etime(clock,t0));



%
%    [group,grpNames,grpRowLoc] = table2gidx(a,groupVars); % leave out categories not present in data
%    n = length(grpNames);

%    function [group,gnames,gloc] = grp2idx(var,varName,reduce)

%    [gnames,gloc,group] = unique(var);
%    gnames = cellstr(gnames);



    