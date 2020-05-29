function varargout=plotClusterLabels(varargin);
% To get help, type plotClusterLabels('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha


declareParameter(...
    'Help', {
        'This script takes a table "tbl" as input and prints a text label per'
        'table row, with the following attributes determined by the'
        'corresponding table entries:'
        '  label position - Determined by tbl.(position)'
        '  label text     - Determined by tbl.(label)'
        '  label size     - Determined by tbl.(sizeValues)'
        '                   The convertions of values in tbl.(sizeValues)'
        '                   to the actual font size of the label can take two forms:'
        '                   . tbl.(sizeValues) contains integer indices that are'
        '                     used to select the sizes from the "fontSizes" array.'
        '                   . tbl.(sizeValues) contains (potentially non-integer)'
        '                     values that are quantized by a given quantization method'
        '                     and the quantized value produces the integer indices'
        '                     that are used to select the sizes from the "fontSizes" array.'
        ' '
        'When the field "cluster_id" is provided, each label corresponds to'
        'one cluster of rows, instead of a single row.'
        ' '
        'When working with clusters, "merging functions" must be provided'
        'to determine how the position/text values of the different'
        'rows in the same cluster are combined to determine the position/text'
        'of the label corresponding to the cluster.'
        'Typically, merging functions sum, average, or count the values.'
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
    'VariableName','positions',...
    'DefaultValue','som_position',...
    'Description', {
        'String containing the name of the table column that define the'
        'position of the labels.'
                   });

declareParameter(...
    'VariableName','labels',...
    'DefaultValue','som_name',...
    'Description', {
        'String containing the name of the table column with the text of'
        'the labels.'
                   });

declareParameter(...
    'VariableName','sizeValues',...
    'DefaultValue','',...
    'Description', {
        'String containing the name of the table column with the values'
        'that define the labels'' font size.'
        'When empty, the array "fontSizes" should have a single row, which'
        'specified the font size of every label.'
                   });

declareParameter(...
    'VariableName','sizeQuantize',...
    'DefaultValue','naturalBreaks',...
    'Description', {
        'When a non-empty string, the values of tbl.(sizeValues) are quantized'
        'using the function quantizeValues(), through the method specified'
        'by this string.'
        'Typically, equal to "naturalBreaks" or "equalBins".'
        ' '
        'See "help quantizeValues".'
                   });
    
declareParameter(...
    'VariableName','fontSizes',...
    'DefaultValue',8,...
    'Description', {
        'Variable defining the font sizes for the labels:'
        '[1x1] arrays     - font size (in pts) for every label'
        '[nSizesx1] array - the font size (in pts) of the i-th row of this'
        '                   table is used, when tbl.(colorValues) (or its'
        '                   quantized value) is equal to the integer i.'
        'Any label with font size 0 will be omited.'
                   });

declareParameter(...
    'VariableName','cluster_id',...
    'DefaultValue','',...
    'Description', {
        'String containing the name of the table column that defines the clusters.'
        ' '
        'Rows with the same value in the field "cluster_id" are assumed to'
        'be in the same cluster. The field "cluster_id" can be of any type'
        '(provided that it can be handled by the "unique" matlab function).'
                   });
    
declareParameter(...
    'VariableName','positionMerging',...
    'DefaultValue','median',...
    'Description', {
        'When "cluster_id" is not empty, this parameter specifies the'
        'function used to merge the values of all the rows of tbl.(positions)'
        'associated with the same cluster.'
        'Typical values include "median", "mean", "first", or a'
        'function handle.'
        ' '
        'See "help tableRowGrouping".'
                   });

declareParameter(...
    'VariableName','labelMerging',...
    'DefaultValue','first',...
    'Description', {
        'When "cluster_id" is not empty, this parameter specifies the'
        'function used to merge the values of all the rows of tbl.(radiusValues)'
        'associated with the same cluster.'
        'Typical values include "first", or a function handle.'
        ' '
        'See "help tableRowGrouping".'
                   });

declareParameter(...
    'VariableName','sizeMerging',...
    'DefaultValue','median',...
    'Description', {
        'When "cluster_id" is not empty, this parameter specifies the'
        'function used to merge the values of all the rows of tbl.(sizeValues)'
        'associated with the same cluster.'
        'Typical values include "median", "mean", "first", or a'
        'function handle.'
        ' '
        'See "help tableRowGrouping".'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','cluster_ids',...
    'Description', {
        'Arrays with the cluster ids.'
                   });
    
declareOutput(...
    'VariableName','labels',...
    'Description',{
        'String array with the cluster labels.'
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

fprintf('plotClusterLabels: starting...');
t0=clock;

fprintf('\nATTENTION: callAndStore() assumes the table has not changed.\n');
fprintf('           When the table changes do "clearAndStore()"\n');

%% Get centers
if isempty(cluster_id)
    centers=callAndStore(...
        @(query)tbl.(positions)(tableSelect(tbl,query),:),query);
    labels=callAndStore(...
        @(query)tbl.(labels)(tableSelect(tbl,query),:),query);
    cluster_ids=[];
else
    centers=callAndStore(...
        @(query,grouping,processing,fun)tableRowGrouping(tbl,query,grouping,processing,fun),...
        query,...         % Row selection query
        cluster_id,...    % Grouping Variables
        positions,...     % Processed Variables
        positionMerging); % Merging function
    
    cluster_ids=callAndStore(...
        @(query,grouping,processing,fun)tableRowGrouping(tbl,query,grouping,processing,fun),...
        query,...         % Row selection query
        cluster_id,...    % Grouping Variables
        cluster_id,...    % Processed Variables
        'first');         % Merging function
    
    labels=callAndStore(...
        @(query,grouping,processing,fun)tableRowGrouping(tbl,query,grouping,processing,fun),...
        query,...      % Document selection query
        cluster_id,... % Grouping Variables
        labels,...     % Processed Variables
        labelMerging); % Merging function
end

%% Get font sizes
if isempty(sizeValues)
    if size(fontSizes,1)~=1
        error('\nplotClusterDisk: empty "sizeValues" but fontSizes has more than one row (%dx%d)\n',size(fontSizes));
    end
    sizeBreakpoints=[];
    sizeValues=fontSizes;
else
    % get values from table
    if isempty(cluster_id)
        sizeValues=callAndStore(...
            @(query)tbl.(sizeValues)(tableSelect(tbl,query),:),query);
    else
        sizeValues=callAndStore(...
            @(query,grouping,processing,fun)tableRowGrouping(tbl,query,grouping,processing,fun),...
            query,...        % Document selection query
            cluster_id,...   % Grouping Variables
            sizeValues,...   % Processed Variables
            sizeMerging);   % function
    end
    % Quantize size
    if isempty(sizeQuantize)
        if any(sizeValues~=round(sizeValues))
            error('\nplotClusterDisk: empty "sizeQuantize" but some size values are not integer-valued.\n');
        end
        if max(sizeValues)>size(fontSizes,1)
            error('\nplotClusterDisk: empty "sizeQuantize" but some size values are larger than the size of "fontSizes" (%dx%d)\n',size(fontSizes));
        end
        if min(sizeValues)<1
            error('\nplotClusterDisk: empty "sizeQuantize" but some size values are smaller than 1\n');
        end

        sizeValues=fontSizes(sizeValues,:);
        sizeBreakpoints=[];
    else
        [sizeValues,sizeBreakpoints]=quantizeValues(sizeValues,fontSizes,sizeQuantize);
    end
end

k=sizeValues==0;
centers(k,:)=[];
labels(k)=[];
sizeValues(k)=[];
if ~isempty(cluster_ids)
    cluster_ids(k,:)=[];
end

%% som cell names
labelChart('centers',centers,'labels',labels,'fontSize',sizeValues,...
           'wrapDelimiter',' ','linesPerLabel',3,'textProperties',{'Color','k'});

fprintf('done plotClusterLabels (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,parameters__);
