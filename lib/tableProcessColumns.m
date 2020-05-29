function tbl=tableProcessColumns(tbl,columnNames,script)
% tbl=tableProcessColumns(tbl,columnNames,script)
% 
% Returns a new table for which the columns in columnNames have been
% processed through the given script. 
%
% When absent script is 
%     @(x)str2double(x)
% which converts a string to a double.
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

if nargin<3
    script=@(x)str2double(x);
end

if ~iscell(columnNames)
    columnNames={columnNames};
end

if ischar(script)
    script=str2func(script);
end

%fprintf('tableProcessColumns: ');
%t0=clock;

for i=1:length(columnNames)
    fprintf('%s ',columnNames{i});
    tbl.(columnNames{i})=script(tbl{:,columnNames{i}});
end

%fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t0));

