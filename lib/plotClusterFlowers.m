function varargout=plotClusterFlowers(varargin);
% To get help, type plotClusterTimeSeries('help')
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
        'This script takes a table "tbl" as input and plots one flower chart'
        'per table row, with the following flower chart attributes determined'
        'by the corresponding table entries:'
        '  center position - Determined by tbl.(positions)'
        '  flower petal radius - Determined by tbl.(radiusValues)'
        '               The convertion of values in tbl.(radiusValues)'
        '               to the actual radius of a petal generally involves'
        '               some form of scaling, determined by "radiusScale".'
        ' '
        'A "flower chart" looks somewhat like a pie chart, with all the "slices"'
        '(called petals) with equal angles, but different radii.'
        ' '
        'When the field "cluster_id" is provided, each flower chart corresponds'
        'to one cluster of rows, instead of a single row.'
        'Rows with the same value in the field "cluster_id" are assumed to'
        'be in the same cluster. The field "cluster_id" can be of any type'
        '(provided that it can be handled by the "unique" matlab function).'
        ' '
        'When working with clusters, "merging functions" must be provided'
        'to determine how the position/radius values of the different'
        'rows in the same cluster are combined to determine the position/radius'
        'of the flower chart corresponding to the cluster.'
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
        'position of the flowers'' centers.'
                   });

declareParameter(...
    'VariableName','radiusValues',...
    'Description', {
        'String containing the name of the table column with the values'
        'that define the petals'' radii.'
                   });

declareParameter(...
    'VariableName','radiusScale',...
    'DefaultValue',{'perceptualDisks',1},...
    'Description', {
        'Scaling applied to values of tbl.(radiusValues)'
        '(possibly merged using "radiusMerging")'
        'to determine the actual radii of petals.'
        'Can either be the structure returned by the function scaleValues(),'
        'or a cell containing {method,range} to be passed to the'
        'function scaleValues().'
        ' '
        'See "help scaleValues".'
                   });
    
declareParameter(...
    'VariableName','radiusNormalize',...
    'DefaultValue',false,...
    'AdmissibleValues',{false,true},...
    'Description', {
        'When true, the values of the petals radii are normalized so that'
        'all 1st petals add up to the same value, all 2nd petals add up'
        'to the same value, etc.'
    });

declareParameter(...
    'VariableName','sharpnessValues',...
    'DefaultValue',0,...
    'Description', {
        'Scalar value that specifies the sharpness of the pie slices,'
        'or string containing the name of the table column with the values'
        'that define the petals'' sharpness.'
        'Positive sharpness values permit the construction of pie charts that'
        'look like flowers:'
        '  0  - regular pie slice'
        '  .1 - pie slice that looks like a fat flower petal'
        '  10 - pie slice that looks like a long and skinny flower petal'
        'Intermediate values are also allowed, producing different types of petals.'
                   });

declareParameter(...
    'VariableName','sharpnessScale',...
    'DefaultValue',{'linear',[]},...
    'Description', {
        'Scaling applied to values of tbl.(sharpnessValues)'
        '(possibly merged using "sharpnessMerging")'
        'to determine the actual radii of petals.'
        'Can either be the structure returned by the function scaleValues(),'
        'or a cell containing {method,range} to be passed to the'
        'function scaleValues().'
        ' '
        'See "help scaleValues".'
                   });
    
declareParameter(...
    'VariableName','colorMap',...
    'DefaultValue',[.8,.8,.8],...
    'Description', {
        'Variable defining the colors used for the flower petals:'
        '[nPetalsx3] array - colors used for the flower petals.'
        ' '
        'When the number of columns returned by radiusValues/radiusMergeing'
        'is larger than the rows of "colorsMap", the excess columns are ignored.'
        'This is useful to ignore the last value returned by histc()'
        'which is typically zero.'
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
    'VariableName','sharpnessMerging',...
    'DefaultValue','median',...
    'Description', {
        'When "cluster_id" is not empty, this parameter specifies the'
        'function used to merge the values of all the rows of tbl.(sharpnessValues)'
        'associated with the same cluster.'
        'Typical values include "median", "mean", "first", or a'
        'function handle.'
        ' '
        'See "help tableRowGrouping".'
                   });

declareParameter(...
    'VariableName','arc0',...
    'DefaultValue',-pi/2,...
    'Description', {
        '1x1 scalar with the angle (in radians) from which arcs are measured.'
        'A zero value for arc0 corresponds to arcs starting from the 3 o''clock,'
        'and a value of -pi/2 corresponds to arcs starting from the 12 o''clock'
        'horizontal line. For normal axis, angles are measure counter clockwise'
        'and for ''axis ij'' they are measured clockwise.'
                   });

declareParameter(...
    'VariableName','trimRadius',...
    'DefaultValue',.03,...
    'Description', {
        'Radius of the trim drawn outside each flower char.'
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
        'When true, draws a solid flower chart with an outer trim.'
        'When false, draws a flower chart ring with width defined by "trimRadius"'
        'and the color of the ring''s interior is defined by "trimColor"'
                   });

declareParameter(...
    'VariableName','halfSlices',...
    'AdmissibleValues',{0,1,-1},...
    'DefaultValue',0,...
    'Description', {
        'Scalar value that specifies the type of pie slice/petal';
        '   0  - regular pie slice/petal'
        '  +1  - clockwise half pie slice/petal'
        '  -1  - anti-clockwise half pie slice/petal'
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
        'Scaling used for the values of the radius, as returned by the function scaleValues().'
        ' '
        'See "help scaleValues".'
        });

declareOutput(...
    'VariableName','radiusValuesUnscaled',...
    'Description',{
        'Values of the radius actually ploted (prior to scaling).'
        'These values may be useful to subsequently plot legends.'
        });

declareOutput(...
    'VariableName','radiusValuesScaled',...
    'Description',{
        'Values of the radius actually ploted (after scaling).'
        'These values may be useful for debuging.'
        });

declareOutput(...
    'VariableName','sharpnessScale',...
    'Description',{
        'Scaling used for the values of the sharpness, as returned by the function scaleValues().'
        ' '
        'See "help scaleValues".'
        });

declareOutput(...
    'VariableName','sharpnessValuesUnscaled',...
    'Description',{
        'Values of the sharpness actually ploted (prior to scaling).'
        'These values may be useful to subsequently plot legends.'
        });

declareOutput(...
    'VariableName','sharpnessValuesScaled',...
    'Description',{
        'Values of the sharpness actually ploted (after to scaling).'
        'These values may be useful for debuging.'
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

fprintf('plotClusterFlower: starting...');
t0=clock;

fprintf('\nATTENTION: callAndStore() assumes the table has not changed.\n');
fprintf('           When the table changes do "clearAndStore()"\n');

%% Get centers & radius
if isempty(cluster_id)
    t=tbl(:,positions);
    fun=@(query,positions)t.(positions)(tableSelect(t,query),:);
    centers=callAndStore(fun,query,positions);
    t=tbl(:,radiusValues);    
    fun=@(query,radiusValues)t.(radiusValues)(tableSelect(t,query),:);
    radiusValues=callAndStore(fun,query,radiusValues);
    if ischar(sharpnessValues)
        t=tbl(:,sharpnessValues);    
        fun=@(query,sharpnessValues)t.(sharpnessValues)(tableSelect(t,query),:);
        sharpnessValues=callAndStore(fun,query,sharpnessValues);
    end
else
    t=tbl(:,{cluster_id,positions});
    fun=@(query,grouping,processing,mfun)...
        tableRowGrouping(t,query,grouping,processing,mfun);
    centers=callAndStore(...
        fun,...
        query,...         % Row selection query
        cluster_id,...    % Grouping Variables
        positions,...     % Processed Variables
        positionMerging); % Merging function;
    if iscell(radiusValues)
        t=tbl(:,[cluster_id,radiusValues]);
    else
        t=tbl(:,{cluster_id,radiusValues});
    end
    fun= @(query,grouping,processing,mfun)...
         tableRowGrouping(t,query,grouping,processing,mfun);
    radiusValues=callAndStore(...
        fun,...
        query,...        % Document selection query
        cluster_id,...   % Grouping Variables
        radiusValues,... % Processed Variables
        radiusMerging);  % Merging function
    if ischar(sharpnessValues)
        t=tbl(:,{cluster_id,sharpnessValues});
        fun=@(query,grouping,processing,mfun)...
            tableRowGrouping(t,query,grouping,processing,mfun);
        sharpnessValues=callAndStore(...
            fun,...
            query,...        % Document selection query
            cluster_id,...   % Grouping Variables
            sharpnessValues,... % Processed Variables
            sharpnessMerging);  % Merging function
    end
end

%% Scale radius
radiusValuesUnscaled=radiusValues;
if radiusNormalize
    radiusValues=normalizeValues(radiusValues',1)';
end
if iscell(radiusScale)
    [radiusValues,radiusScale]=scaleValues(radiusValues,radiusScale{:});
else
    [radiusValues,radiusScale]=scaleValues(radiusValues,radiusScale);
end
if radiusNormalize
    radiusScales=nan;
end

%% Scale sharpness
sharpnessValuesUnscaled=sharpnessValues;
if iscell(sharpnessScale)
    [sharpnessValues,sharpnessScale]=scaleValues(sharpnessValues,sharpnessScale{:});
else
    [sharpnessValues,sharpnessScale]=scaleValues(sharpnessValues,sharpnessScale);
end


%% Get colors

if size(colorMap,1)<size(radiusValues,2)
    fprintf('plotClusterFlower: removing last %d values of the time radius, because # colors %d < # values %d\n',size(radiusValues,2)-size(colorMap,1),size(colorMap,1),size(radiusValues,2));
    radiusValues=radiusValues(:,1:size(colorMap,1));
elseif size(colorMap,1)>size(radiusValues,2)
    error('plotClusterFlower: number of colors (%d) larger than the number of values in the radius (%d)\n',size(colorMap,1),size(radiusValues,2));
end

%% Get arcs
nSlices=size(colorMap,1);
arcs=normalizeValues(ones(1,nSlices),2*pi);

%% Plot
if solid
    if trimRadius>0
        bandCharts('centers',centers,...
                   'radii',{radiusValues+trimRadius,-trimRadius},...
                   'arc0',arc0,...
                   'arcs',arcs,...
                   'slicesSharpness',{sharpnessValues,sharpnessValues},...
                   'halfSlices',halfSlices,...
                   'colors',{trimColor,colorMap},...
                   'pointsPerCircle',pointsPerCircle,'patchProperties',{'lineStyle','none'});
    else
        bandCharts('centers',centers,...
                   'radii',radiusValues,...
                   'arc0',arc0,...
                   'arcs',arcs,...
                   'slicesSharpness',{sharpnessValues,sharpnessValues},...
                   'halfSlices',halfSlices,...
                   'colors',colorMap,...
                   'pointsPerCircle',pointsPerCircle,'patchProperties',{'lineStyle','none'});
    end
else
    bandCharts('centers',centers,...
               'radii',{radiusValues,-trimRadius},...
               'arc0',arc0,...
               'arcs',arcs,...
               'slicesSharpness',sharpnessValues,...
               'halfSlices',halfSlices,...
               'colors',{colorMap,trimColor},...
               'pointsPerCircle',pointsPerCircle,'patchProperties',{'lineStyle','none'});
end

radiusValuesScaled=radiusValues;
sharpnessValuesScaled=sharpnessValues;

fprintf('done plotClusterFlower (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,parameters__);

