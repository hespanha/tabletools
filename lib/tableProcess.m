function varargout=tableProcess(varargin);
% To get help, type tableProcess('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This script processes a table class object by executing:'
        '1) performing a regular expression substitution on a given set of columns'
        '2) remove duplicate table rows based on a given set of columns'
        '3) sorting the table rows based on a given set of columns'
        '4) processing a given set of columns using a given script'
        '5) selecting a set of rows to keep based on a query'
        '6) creates a new column based on a given expression'
%        '7) randomly subsamples the rows of the table'
        'The operations are performed in this order.'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','tbl',...
    'Description', {
        'Table to be processed.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Substitution
declareParameter(...
    'VariableName','regexprepColumns',...
    'DefaultValue','',...
    'Description', {
        '1) For regular expression substitution:'
        'Columns on which to perform a regular expression substitution.'
        'Multiple columns can be specified by setting this variable to'
        'be a cell array with multiple column names.'
                   });
declareParameter(...
    'VariableName','expression',...
    'DefaultValue','',...
    'Description', {
        '1) For regular expression substitution:'
        'Regular expression used in the regular expression substitution.'
        'See ''help regexprep.'''
                   });
declareParameter(...
    'VariableName','repl',...
    'DefaultValue','',...
    'Description', {
        '1) For regular expression substitution:'
        'Replacement expression used in the regular expression substitution.'
        'See ''help regexprep.'''
                   });

%% Unique
declareParameter(...
    'VariableName','uniqueColumns',...
    'DefaultValue','',...
    'Description', {
        '2) For removing duplicate rows:'
        'Columns that are considered to determine if two rows are equal.'
        'Sorting can involve multiple columns by setting this variable'
        'to be a cell array with multiple column names.'
                   });

%% Sort
declareParameter(...
    'VariableName','sortColumns',...
    'DefaultValue','',...
    'Description', {
        '3) For row sorting:'
        'Columns by which table should be sorted.'
        'Sorting can involve multiple columns by setting this variable'
        'to be a cell array with multiple column names.'
                   });
declareParameter(...
    'VariableName','sortMode',...
    'DefaultValue','ascend',...
    'Description', {
        '3) For row sorting:'
        'Determines the direction of sorting. May be either of the strings'
        '''ascend'' or ''descend'' or a cell array with the same size as'
        '''sortColumns'' to specify different directions for the different columns.'
                   });

%% Script
declareParameter(...
    'VariableName','scriptColumns',...
    'DefaultValue','',...
    'Description', {
        '4) For running script on columns:'
        'Column to process using a given script.'
        'Multiple columns can be specified by setting this variable to'
        'be a cell array with multiple column names.'
        'See ''help tableProcessColumns''.'
        });
declareParameter(...
    'VariableName','processScript',...
    'DefaultValue','',...
    'Description', {
        '4) For running script on columns:'
        'Script to be used in processing the columns specified by ''columns2process.'''
        'E.g., convert a string to a double, use:'
        '    @(str)str2double(str)'
        'and to convert a string (US format) to a date, use'
        '    @(str)extractDateUS(str)'
        'See ''help tableProcessColumns''.'
        });

%% Remove rows
declareParameter(...
    'VariableName','query',...
    'DefaultValue','',...
    'Description', {
        '5) For removing rows:'
        'Query specifing a subset of rows to keep.'
        'A query is a function whose input is the table and that returns'
        'which rows should be keep. E.g., To find the rows whose length'
        'of the field ''name'' is smaller than 20, use:'
        '  query = @(tbl) length(tbl.name)<20',
        ' '
        'The query is executed after sorting and processing the columns.'
        'See ''help tableSelect'' for further details on the query format.'
                   });

%% Add column
declareParameter(...
    'VariableName','newColumnName',...
    'DefaultValue','',...
    'Description', {
        '6) For adding columns:'
        'Name of the column to be created'
        });
declareParameter(...
    'VariableName','newColumnExpression',...
    'DefaultValue','',...
    'Description', {
        '6) For adding columns:'
        'Expression used to compute values to the new column.'
        'The expression should be a function whose input is'
        'the table and that returns the value of the new column.'
        'E.g., to create a column that is true when the table variable'
        '''name'' is empty one could use:'
        '     @(tbl)cellfun(''isempty'',tbl.name)'
        });
    
%%% Checks
declareParameter(...
    'VariableName','skipCheckClasses',...
    'DefaultValue',false,...
    'AdmissibleValues',{false,true},...
    'Description', {
        'When true, skips checking whether table columns are of a consistent types.'
        });

% $$$ %% Random sample
% $$$ declareParameter(...
% $$$     'VariableName','sampleSize',...
% $$$     'DefaultValue',inf,...
% $$$     'Description', {
% $$$         '7) For randomly sampling the rows:'
% $$$         'Maximum number of rows to keep in the random selection.'
% $$$         'Rows are randomly selected (without repetitions).'
% $$$         });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','tbl',...
    'Description', {
        'Table that was processed.'
                   });
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,params]=setParameters(nargout,varargin);
if stopNow
    return 
end

%verboseLevel=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('tableProcess:\n');
t0=clock;

%% regexprep
if ~isempty(regexprepColumns)
    fprintf('  regexprep...');
    t1=clock;
    tbl=tableRegexprep(tbl,regexprepColumns,expression,repl);
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
end

%% sort
if ~isempty(uniqueColumns)
    fprintf('  removing duplicates...');
    t1=clock;
    [~,ia,ic]=unique(tbl(:,uniqueColumns),'stable');
    tbl=tbl(ia,:);
    fprintf('done removed %d/%d rows (%d rows, %d columns, %.3f sec)\n',...
            length(ic)-length(ia),length(ic),size(tbl,1),size(tbl,2),etime(clock,t1));
end

%% sort
if ~isempty(sortColumns)
    fprintf('  sorting...');
    t1=clock;
    tbl=sortrows(tbl,sortColumns,sortMode);
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
end

%% process columns
if ~isempty(scriptColumns)
    fprintf('  script...');
    t1=clock;
    tbl=tableProcessColumns(tbl,scriptColumns,processScript);
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
end

%% query
if ~isempty(query)
    fprintf('  query...');
    t1=clock;
    k=tableSelect(tbl,query);
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
    fprintf('  removing %d columns...',size(tbl,1)-length(k));
    t1=clock;
    tbl=tbl(k,:);
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
end

%% add column
if ~isempty(newColumnName)
    fprintf('  new column...');
    t1=clock;
    tbl.(newColumnName)=newColumnExpression(tbl);
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
end

% $$$ %% sample rows
% $$$ if ~isinf(sampleSize)
% $$$     fprintf('  subsample...');
% $$$     t1=clock;
% $$$     sampleSize=min(sampleSize,size(tbl,1));
% $$$     k=randperm(size(tbl,1));
% $$$     k=sort(k(1:sampleSize));
% $$$     tbl=tbl(k,:);
% $$$     fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
% $$$ end

if ~skipCheckClasses
    %% Check types
    tbl=tableCheckClasses(tbl);
end

fprintf('done tableProcess (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);

end

