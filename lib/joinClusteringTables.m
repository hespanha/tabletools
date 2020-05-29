function varargout=joinClusteringTables(varargin);
% To get help, type joinClusteringTables('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This script joins'
        '1) a base table'
        '2) multiple tables with the features that have been used'
        '   for clustering the rows of the base table'
        '3) multiple tables with row-clustering results'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','baseTable',...
    'Description', {
        'Base table with the items that have been clustered.'
                   });
    
declareParameter(...
    'VariableName','variableNames',...
    'DefaultValue',[],...
    'Description', {
        'Names of the columns from the ''baseTable'' to be included in'
        'the joint table. When empty (i.e.,[]) all columns are included.'
                   });

declareParameter(...
    'VariableName','query',...
    'DefaultValue','',...
    'Description', {
        'Query specifing a subset of rows of the base table to keep.'
        'A query is a function whose input is the table and that returns'
        'which rows should be keep. E.g., To find the rows whose length'
        'of the field ''name'' is smaller than 20, use:'
        '  query = @(tbl) length(tbl.name)<20',
        ' '
        'See ''help tableSelect'' for further details on the query format.'
                   });

declareParameter(...
    'VariableName','featuresTables',...
    'DefaultValue',{},...
    'Description', {
        'Cell array of tables with the features that have been used'
        'for clustering the rows of the ''baseTable''.'
                   });

declareParameter(...
    'VariableName','featuresNames',...
    'DefaultValue',{},...
    'Description', {
        'Array of string with names for the ''featuresTables'''
        'This array should have the same length as ''featuresTables'''
        'These names will be used as prefixes to the names of the'
        'columns of the ''featuresTables'' when inserted into the joint table.'
                   });

declareParameter(...
    'VariableName','clusteringTables',...
    'DefaultValue',{},...
    'Description', {
        'Cell array of tablescluster objects with the clustering of the rows'
        'of the ''baseTable''.'
                   });

declareParameter(...
    'VariableName','clusteringNames',...
    'DefaultValue',{},...
    'Description', {
        'Array of string with names for the ''clusteringTables'''
        'This array should have the same length as ''clusteringTables'''
        'These names will be used as prefixes to the names of the'
        'columns of the ''clusteringTables'' when inserted into the joint table.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','baseTable',...
    'Description', {
        'Joint table.'
                   });
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters and inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,params]=setParameters(nargout,varargin);
if stopNow
    return 
end

%verboseLevel=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~iscell(featuresTables)
    featuresTables={featuresTables};
end
if ~iscell(featuresNames)
    featuresNames={featuresNames};
end
if ~iscell(clusteringTables)
    clusteringTables={clusteringTables};
end
if ~iscell(clusteringNames)
    clusteringNames={clusteringNames};
end

fprintf('joinClusteringTables:\n');
t0=clock();

if ~isempty(query)
    k=tableSelect(baseTable,query);
    fprintf('  selecting rows in baseTable(%dx%d): ',size(baseTable));
    t1=clock;
    baseTable=baseTable(k,:);
    fprintf('  done baseTable(%dx%d) (%.3f sec)\n',size(baseTable),etime(clock,t1));
end

if ~isempty(variableNames)
    if ~ismember('primaryKey',variableNames)
        variableNames={'primaryKey',variableNames{:}};
    end

    fprintf('  selecting %d columns in baseTable(%dx%d): ',length(variableNames),size(baseTable));
    fprintf('%s ',variableNames{:});
    t1=clock;
    VNs=intersect(variableNames,baseTable.Properties.VariableNames);
    if length(VNs)<length(variableNames)
        ignored=setdiff(variableNames,VNs);
        fprintf('\n    WARNING: %d fields do not exist: ',length(ignored));
        fprintf('%s ',ignored{:});
        fprintf('\n');
    end
    baseTable=baseTable(:,VNs);
    fprintf('  done baseTable(%dx%d) (%.3f sec)\n',size(baseTable),etime(clock,t1));
end

% Join tables with features
for i=1:length(featuresTables)
    t=featuresTables{i};
    %  Rename columns to avoid conflicts:
    fprintf('  renaming variables in %s(%dx%d)...',featuresNames{i},size(t));
    t1=clock;
    t.Properties.VariableNames=cellfun(@(x)sprintf('%s_%s',featuresNames{i},x),...
                                           t.Properties.VariableNames,'UniformOutput',0);
    fprintf('  done %s(%dx%d) (%.3f sec)\n',featuresNames{i},size(t),etime(clock,t1));
    % innerjoin
    fprintf('  innerjoin %s(%dx%d) to baseTable(%dx%d)...',featuresNames{i},size(t),size(baseTable));
    t1=clock;
    baseTable=innerjoin(baseTable,t,'LeftKeys','primaryKey','RightKeys',sprintf('%s_primaryKey',featuresNames{i}));
    fprintf('  done baseTable(%dx%d) (%.3f sec)\n',size(baseTable),etime(clock,t1));
end

% Joing tables with clustering
for i=1:length(clusteringTables)
    t=clusteringTables{i};

    %  Rename columns to avoid conflicts:
    fprintf('  renaming variables in %s.tbl(%dx%d)...',clusteringNames{i},size(t.tbl));
    t1=clock;
    t.tbl.Properties.VariableNames=cellfun(@(x)sprintf('%s_%s',clusteringNames{i},x),...
                                           t.tbl.Properties.VariableNames,'UniformOutput',0);
    fprintf('  done %s.tbl(%dx%d) (%.3f sec)\n',clusteringNames{i},size(t.tbl),etime(clock,t1));
    fprintf('  renaming variables in %s.clusters(%dx%d)...',clusteringNames{i},size(t.clusters));
    t1=clock;
    t.clusters.Properties.VariableNames=cellfun(@(x)sprintf('%s_%s',clusteringNames{i},x),...
                                           t.clusters.Properties.VariableNames,'UniformOutput',0);
    fprintf('  done %s.clusters(%dx%d) (%.3f sec)\n',clusteringNames{i},size(t.clusters),etime(clock,t1));
    % innerjoin
    fprintf('  innerjoin %s(%dx%d) to baseTable(%dx%d)...',clusteringNames{i},size(t.tbl),size(baseTable));
    t1=clock;
    baseTable=innerjoin(baseTable,t.tbl,'leftKeys','primaryKey','rightKeys',sprintf('%s_primaryKey',clusteringNames{i}));
    fprintf('  done baseTable(%dx%d) (%.3f sec)\n',size(baseTable),etime(clock,t1));
    % outerjoin: with left join, so that clusters without documents are ignored
    fprintf('  outerjoin %s(%dx%d) to baseTable(%dx%d)...',clusteringNames{i},size(t.clusters),size(baseTable));
    t1=clock;
    baseTable=outerjoin(baseTable,t.clusters,...
                'leftKeys',sprintf('%s_id',clusteringNames{i}),...
                'rightKeys',sprintf('%s_id',clusteringNames{i}),...
                'type','left','MergeKeys',true);
    fprintf('  done baseTable(%dx%d) (%.3f sec)\n',size(baseTable),etime(clock,t1));
end
fprintf('  done joinClusteringTables (%.3f sec)\n',etime(clock,t0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);

