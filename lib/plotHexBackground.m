function varargout=plotHexBackground(varargin);
% To get help, type plotHexBackground('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha


declareParameter(...
    'Help', {
        'This script plots an hexagonal SOM grid.'
        'The color of the grid hexagons can reflect additional clustering information'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','tbl',...
    'Description', {
        'Input table containing the features and SOM clustering data.'
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
    'VariableName','somPositions',...
    'DefaultValue','som_position',...
    'Description', {
        'String containing the name of the table column with the SOM cell positions.'
                   });

declareParameter(...
    'VariableName','cluster_id',...
    'DefaultValue','',...
    'Description', {
        'String containing the name of the table column with the cluster'
        'ids (from 1 to nClusters) used to color the hexagonal grid.'
        'The color each SOM hexagon is determined by which cluster contains'
        'more rows of the table in that cluster.'
                   });
    
declareParameter(...
    'VariableName','colors',...
    'DefaultValue',[.8,.8,.8],...
    'Description', {
        'Variable defining the color used for the hexagonal grid:'
        '[1x3] arrays        - color to be used for every hexagon.'
        '[nClustersx3] array - colors to be used for the hexagons of each cluster.'
        '                      The variable ''clustering'' must be defined.'
                   });

declareParameter(...
    'VariableName','nTicks',...
    'DefaultValue',8,...
    'AdmissibleValues',{0,1,2,3,4,5,6,7,8,9,10},...
    'Description', {
        'Number of (equaly-spaced) latitude/longitude labels printed in the'
        'margins of the SOM.'
        'The labels will be single capital letters for longitude,'
        'and signle numerical digits for latitute.'
                   });

declareParameter(...
    'VariableName','fontSize',...
    'DefaultValue',10,...
    'Description', {
        'Variable defining the font size for the latitude/longitude labels.'
                   });

declareParameter(...
    'VariableName','lineWidth',...
    'DefaultValue',.04,...
    'Description', {
        'Line width for the hexagons'
                   });


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,parameters__]=setParameters(nargout,varargin);
if stopNow
    if nargout>0
        varargout{1}=parameters__;
    end
    return 
end

%verboseLevel=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('plotHexBackground: starting...');
t0=clock;

som.centers=callAndStore(...
    @(query,grouping,processing,fun)tableRowGrouping(tbl,query,grouping,processing,fun),...
    query,...           % Row selection query
    {somPositions},...  % Grouping Variables
    {somPositions},...  % Processed Variables
    'first');          % function

if size(colors,1)>1 && isempty(cluster_id)
    error('plotHexBackground: %d(>1) colors given but ''cluster_id'' has not been specified\n',size(colors,1));
end

if ~isempty(cluster_id)
    nClusters=max(tbl.(cluster_id));

    if nClusters ~= size(colors,1)
        error('plotHexBackground: mismatch between number of clusters in tbl.%s (%d) and number of rows in ''colors''\n',cluster_id,nClusters,size(colors,1));
    end

    som.kmeans_hist=callAndStore(...
        @(query,grouping,processing,fun)tableRowGrouping(tbl,query,grouping,processing,fun),...
        query,...                  % Row selection query
        {somPositions},...         % Grouping Variables
        {cluster_id},...           % Processed Variables
        @(x)hist(x,1:nClusters)'); % function
    [~,som.cluster]=max(som.kmeans_hist,[],2);

    colors=colors(som.cluster,:);
end


%% background hexs color - cluster id
bandCharts('centers',som.centers,...
           'radii',{.5/cos(pi/6),-lineWidth},'colors',{colors,[1,1,1]},...
           'arc0',-pi/2,...
           'pointsPerCircle',7,'patchProperties',{'lineStyle','none'});

%% axis ticks
axisTicks=[(min(som.centers(:,1))-1/nTicks*(max(som.centers(:,1))-min(som.centers(:,1))))*ones(nTicks,1),linspace(min(som.centers(:,2)),max(som.centers(:,2)),nTicks)';
           linspace(min(som.centers(:,1)),max(som.centers(:,1)),nTicks)',(min(som.centers(:,2))-1/nTicks/1.5*(max(som.centers(:,2))-min(som.centers(:,2))))*ones(nTicks,1)];
textTicks=cellstr(char(['1'+(0:nTicks-1)';
                    'A'+(0:nTicks-1)']));

labelChart('centers',axisTicks,'labels',textTicks,'fontSize',fontSize,'textProperties',{'Color','k'});

fprintf('done plotHexBackground (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,parameters__);
