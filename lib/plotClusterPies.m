function varargout=plotClusterPies(varargin);
% To get help, type plotClusterPies('help')
%
% Copyright (C) 2014 Joao Hespanha

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

declareParameter(...
    'Help', {
        'This script takes a table "tbl" as input and plots one pie chart'
        'per table row, with the following pie attributes determined by the'
        'corresponding table entries:'
        '  center position - Determined by tbl.(positions)'
        '  pie radius - Determined by tbl.(radiusValues)'
        '               The convertion of values in tbl.(radiusValues)'
        '               to the actual radius of a pie generally involves'
        '               some form of scaling, determined by "radiusScale".'
        '  pie slice angles - Determined by tbl.(pieValues)'
        ' '
        'When the field "cluster_id" is provided, each pie chart corresponds'
        'to one cluster of rows, instead of a single row.'
        'Rows with the same value in the field "cluster_id" are assumed to'
        'be in the same cluster. The field "cluster_id" can be of any type'
        '(provided that it can be handled by the "unique" matlab function).'
        ' '
        'When working with clusters, "merging functions" must be provided'
        'to determine how the position/radius values of the different'
        'rows in the same cluster are combined to determine the position/radius'
        'of the pie chart corresponding to the cluster.'
        'Typically, merging functions sum, average, or count the values.'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','tbl',...
    'Description', {
        'Input table containing the clustered features.'
                   });

declareParameter(...
    'VariableName','query',...
    'DefaultValue','',...
    'Description', {
        'Query specifing a subset of rows of the base table to keep.'
        'A query is a function whose input is the table and that returns'
        'which rows should be keep. E.g., To find the rows whose length'
        'of the field "name" is smaller than 20, use:'
        '  query = @(tbl) length(tbl.name)<20',
        ' '
        'See "help tableSelect" for further details on the query format.'
                   });

declareParameter(...
    'VariableName','positions',...
    'DefaultValue','som_position',...
    'Description', {
        'String containing the name of the table column that define the'
        'position of the pies'' centers.'
                   });

declareParameter(...
    'VariableName','radiusValues',...
    'Description', {
        'String containing the name of the table column with the values'
        'that define the pies'' radii.'
        'Instead, it can be also be a numerical value, in which case this'
        'value is the radius for every pie.'
                   });

declareParameter(...
    'VariableName','radiusScale',...
    'DefaultValue',{'perceptualDisks',1},...
    'Description', {
        'Scaling applied to values of tbl.(radiusValues)'
        '(possibly merged using "radiusMerging")'
        'to determine the actual radius of pies.'
        'Can either be the structure returned by the function scaleValues(),'
        'or a cell containing {method,range} to be passed to the'
        'function scaleValues().'
        ' '
        'See "help scaleValues".'
                   });
    
declareParameter(...
    'VariableName','pieValues',...
    'Description', {
        'String containing the name of the table column with the values'
        'that define the angles of the pie slices.'
                   });

declareParameter(...
    'VariableName','colorMap',...
    'DefaultValue',[.8,.8,.8],...
    'Description', {
        'Variable defining the colors used for the pie slices:'
        '[nSlicesx3] array - colors used for the pie slices.'
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
    'VariableName','radiusMerging',...
    'DefaultValue','median',...
    'Description', {
        'When "cluster_id" is not empty, this parameter specifies the'
        'function used to merge the values of all the rows of tbl.(radiusValues)'
        'associated with the same cluster.'
        'Typical values include "median", "mean", "first", or a'
        'function handle.'
        ' '
        'See "help tableRowGrouping".'
                   });

declareParameter(...
    'VariableName','pieMerging',...
    'DefaultValue','median',...
    'Description', {
        'When "cluster_id" is not empty, this parameter specifies the'
        'function used to merge the values of all the rows of tbl.(pieValues)'
        'associated with the same cluster.'
        'Typical values include "median", "mean", "first", or a'
        'function handle.'
        ' '
        'See "help tableRowGrouping".'
                   });

declareParameter(...
    'VariableName','trimRadius',...
    'DefaultValue',.03,...
    'Description', {
        'Radius of the trim drawn outside each pie.'
        'When zero, no trim is drawn.'
                   });

declareParameter(...
    'VariableName','trimColor',...
    'DefaultValue',[1,1,1],...
    'Description', {
        '[1x3] vector with the RGB colors of the outer trim.'
                   });

declareParameter(...
    'VariableName','solid',...
    'DefaultValue',true,...
    'AdmissibleValues',{false,true},...
    'Description', {
        'When true, draws a solid pie with an outer trim.'
        'When false, draws a pie ring with width defined by "trimRadius" and the'
        'color of the ring''s interior is defined by "trimColor"'
                   });

declareParameter(...
    'VariableName','pointsPerCircle',...
    'DefaultValue',72,...
    'Description', {
        'Scalar with the number of points used to draw the curved side of an whole.'
        'circle. A large number will result in a smoother curve, but slower drawing.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','radiusScale',...
    'Description',{
        'Scaling used for the radius, as returned by the function scaleValues().'
        ' '
        'See "help scaleValues".'
        });

declareOutput(...
    'VariableName','radiusValuesUnscaled',...
    'Description',{
        'Values of the radius actually ploted (prior to scaling).'
        'These values may be useful to subsequently plot legends.'
        });


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,parameters__]=setParameters(nargout,varargin);
if stopNow
    if nargout>0
        varargout{1}=parameters__;
    end
    return 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('plotClusterPies: starting...');
t0=clock;

fprintf('\nATTENTION: callAndStore() assumes the table has not changed.\n');
fprintf('           When the table changes do "clearAndStore()"\n');

%% Get centers
if isempty(cluster_id)
    centers=callAndStore(...
        @(query)tbl.(positions)(tableSelect(tbl,query),:),query);
else
    centers=callAndStore(...
        @(query,grouping,processing,fun)tableRowGrouping(tbl,query,grouping,processing,fun),...
        query,...         % Row selection query
        cluster_id,...    % Grouping Variables
        positions,...     % Processed Variables
        positionMerging); % Merging function
end

%% Get radius
if ~isnumeric(radiusValues)
    % get values from table
    if isempty(cluster_id)
        radiusValues=callAndStore(...
            @(query)tbl.(radiusValues)(tableSelect(tbl,query),:),query);
    else
        radiusValues=callAndStore(...
            @(query,grouping,processing,fun)tableRowGrouping(tbl,query,grouping,processing,fun),...
            query,...        % Document selection query
            cluster_id,...   % Grouping Variables
            radiusValues,... % Processed Variables
            radiusMerging);  % Merging function
    end
end
radiusValuesUnscaled=radiusValues;
% Scale radius
if iscell(radiusScale)
    [radiusValues,radiusScale]=scaleValues(radiusValues,radiusScale{:});
else
    [radiusValues,radiusScale]=scaleValues(radiusValues,radiusScale);
end

%% Get angles for the pie slices
if isempty(cluster_id)
    pieValues=callAndStore(...
        @(query)tbl.(pieValues)(tableSelect(tbl,query),:),query);
else
    pieValues=callAndStore(...
        @(query,grouping,processing,fun)tableRowGrouping(tbl,query,grouping,processing,fun),...
        query,...        % Document selection query
        cluster_id,...   % Grouping Variables
        pieValues,...    % Processed Variables
        pieMerging);     % function
    pieValues=normalizeValues(pieValues,2*pi);
end

%% Plot
if solid
    if trimRadius>0
        bandCharts('centers',centers,...
                   'radii',{radiusValues+trimRadius,-trimRadius},...
                   'arcs',pieValues,...
                   'colors',{trimColor,colorMap},...
                   'pointsPerCircle',pointsPerCircle,'patchProperties',{'lineStyle','none'});
    else
        bandCharts('centers',centers,...
                   'radii',radiusValues,...
                   'arcs',pieValues,...
                   'colors',colorMap,...
                   'pointsPerCircle',pointsPerCircle,'patchProperties',{'lineStyle','none'});
    end
else
    bandCharts('centers',centers,...
               'radii',{radiusValues,-trimRadius},...
               'arcs',pieValues,...
               'colors',{colorMap,trimColor},...
               'pointsPerCircle',pointsPerCircle,'patchProperties',{'lineStyle','none'});
end

fprintf('done plotClusterPies (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,parameters__);

