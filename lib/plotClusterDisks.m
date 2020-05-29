function varargout=plotClusterDisks(varargin);
% To get help, type plotClusterDisks('help')
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
        'This script takes a table "tbl" as input and plots one disk per'
        'table row, with the following disk attributes determined by the'
        'corresponding table entries:'
        '  center position - Determined by tbl.(position)'
        '  disk radius - Determined by tbl.(radiusValues)'
        '                The convertion of values in tbl.(radiusValues)'
        '                to the actual radius of a disk generally involves'
        '                some form of scaling, determined by "radiusScale".'
        '  disk color  - Determined by tbl.(colorValues)'
        '                The convertions of values in tbl.(colorValues)'
        '                to the actual RGB color of the disk can take two forms:'
        '                . tbl.(colorValues) contains integer indices that are'
        '                  used to select the colors from the given colormap.'
        '                . tbl.(colorValues) contains (potentially non-integer)'
        '                  values that are quantized by a given quantization method'
        '                  and the quantized value produces the integer indices'
        '                  that are used to select the colors from the given colormap.'
        ' '
        'When the field "cluster_id" is provided, each disk corresponds to'
        'one cluster of rows, instead of a single row.'
        ' '
        'When working with clusters, "merging functions" must be provided'
        'to determine how the position/radius/color values of the different'
        'rows in the same cluster are combined to determine the position/radius/color'
        'of the disk corresponding to the cluster.'
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
        'position of the disks'' centers.'
                   });

declareParameter(...
    'VariableName','noiseRadius',...
    'DefaultValue',0,...
    'Description', {
        'When nonzero, adds to the positions of the disks'' centers random value'
        'that is uniformly distributed in a circle of the given radius.'
        'This parameter is useful when the positions are given by SOM, where all'
        'the table rows corresponding to the the same cell are at the same position.'
        'For an hexagonal SOM use:'
        '. noiseRadius=.5/cos(pi/6) for cloud in circle fully containing the hexagon.'
        '. noiseRadius=.5 for cloud in circle fully contained inside the hexagon.'
                   });

declareParameter(...
    'VariableName','radiusValues',...
    'Description', {
        'String containing the name of the table column with the values'
        'that define the disks'' radii.'
        'Instead, it can be also be a numerical value, in which case this'
        'value is the radius for every disk.'
        ' '
        'In general, disks with zero radius are not draw at all. However,'
        'when *all* the disk radii turn out to be zero, a "cloud" of points is'
        'ploted (each point is actually a very small pentagon).'
                   });

declareParameter(...
    'VariableName','radiusScale',...
    'DefaultValue',{'perceptualDisks',1},...
    'Description', {
        'Scaling applied to values of tbl.(radiusValues)'
        '(possibly merged using "radiusMerging")'
        'to determine the actual radii of the disks.'
        'Can either be the structure returned by the function scaleValues(),'
        'or a cell containing {method,range} to be passed to the'
        'function scaleValues().'
        ' '
        'See "help scaleValues".'
                   });
    
declareParameter(...
    'VariableName','colorValues',...
    'DefaultValue','',...
    'Description', {
        'String containing the name of the table column with the values'
        'that define the disks'' colors.'
        'When empty, the "colorMap" should have a single row, which'
        'specified the color of every disk.'
                   });

declareParameter(...
    'VariableName','colorQuantize',...
    'DefaultValue','',...
    'Description', {
        'When a non-empty string, the values of tbl.(colorValues) are quantized'
        'using the function quantizeValues(), through the method specified'
        'by this string.'
        'Typically, equal to "naturalBreaks" or "equalBins".'
        ' '
        'See "help quantizeValues".'
                   });
    
declareParameter(...
    'VariableName','colorMap',...
    'DefaultValue',[.8,.8,.8],...
    'Description', {
        'Variable defining all the colors used for the disks.'
        '[1x3] arrays      - color used for every disk.'
        '[nColorsx3] array - the i-th row of this table is used, when'
        '                    tbl.(colorValues) (or its quantized value)'
        '                    is equal to the integer i.'
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
        'function used to merge the values of all the rows of tbl.(position)'
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
    'VariableName','colorMerging',...
    'DefaultValue','median',...
    'Description', {
        'When "cluster_id" is not empty, this parameter specifies the'
        'function used to merge the values of all the rows of tbl.(colorValues)'
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
        'Radius of the trim drawn outside each disk.'
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
        'When true, draws a solid disk with an outer trim.'
        'When false, draws a ring with width defined by "trimRadius" and the'
        'color of the ring''s interior is defined by "trimColor"'
        
                   });

declareParameter(...
    'VariableName','pointsPerCircle',...
    'DefaultValue',72,...
    'Description', {
        'Scalar with the number of points used to draw the curved side of an whole.'
        'circle. A large number will result in a smoother curve, but slower drawing.'
        ' '
        'When *all* the disk radii turn out to be zero, this value is ignored'
        'and replaced by 6 (leading to small pentagons).'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','radiusScale',...
    'Description',{
        'Scaling used for the radius, as returned by the function scaleValues().'
        'Can be used for subsequent calls to'
        '   plotCLusterDisks(), plotClusterPies(), plotClusterFlowers()'
        'to maintain the same scale.'
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
    'VariableName','colorBreakpoints',...
    'Description',{
        'Breakpoints used for the color quantization, as returned by quantizeValues().'
        ' '
        'See "help scaleValues".'
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

fprintf('plotClusterDisks: starting...');
t0=clock;

fprintf('\nATTENTION: callAndStore() assumes the table has not changed.\n');
fprintf('           When the table changes do "clearAndStore()"\n');

%% Get centers
if isempty(cluster_id)
    if isempty(query)
        t=tbl(:,{positions});
    else
        t=tbl;  % very wasteful, but query may require other fields
    end
    fun=@(query,positions)t.(positions)(tableSelect(t,query),:);
    centers=callAndStore(fun,query,positions);
    if noiseRadius>0 
        % make noise repeatable for each row
        s = RandStream('mt19937ar','Seed',1);
        r=rand(s,size(centers,1),1);
        th=2*pi*rand(s,size(centers,1),1);
        centers=centers+noiseRadius*[sqrt(r).*cos(th),sqrt(r).*sin(th)];
    end
else
    if isempty(query)
        t=tbl(:,{cluster_id,positions});
    else
        t=tbl;  % very wasteful, but query may require other fields
    end
    fun=@(query,grouping,processing,mfun)tableRowGrouping(t,query,grouping,processing,mfun);
    centers=callAndStore(...
        fun,...
        query,...         % Row selection query
        cluster_id,...    % Grouping Variables
        positions,...     % Processed Variables
        positionMerging); % Merging function

    if noiseRadius>0
        % make noise repeatable for each row
        s = RandStream('mt19937ar','Seed',1);
        r=rand(s,size(centers,1),1);
        th=2*pi*rand(size(centers,1),1);
        centers=centers+noiseRadius*[sqrt(r).*cos(th),sqrt(r).*sin(th)];
    end
end

%% Get radius
if isnumeric(radiusValues)
    radiusValuesUnscaled=radiusValues;
else
    % get values from table
    if isempty(cluster_id)
        if isempty(query)
            t=tbl(:,radiusValues);
        else
            t=tbl;  % very wasteful, but query may require other fields
        end
        fun=@(query,radiusValues)t.(radiusValues)(tableSelect(t,query),:);
        radiusValues=callAndStore(fun,query,radiusValues);
    else
        if isempty(query)
            if iscell(radiusValues)
                t=tbl(:,[cluster_id,radiusValues]);
            else
                t=tbl(:,{cluster_id,radiusValues});
            end
        else
            t=tbl;  % very wasteful, but query may require other fields
        end
        fun=@(query,grouping,processing,mfun)tableRowGrouping(t,query,grouping,processing,mfun);
        radiusValues=callAndStore(...
            fun,...
            query,...        % Document selection query
            cluster_id,...   % Grouping Variables
            radiusValues,... % Processed Variables
            radiusMerging);  % Merging function
    end
    radiusValuesUnscaled=radiusValues;
end
% Scale radius
if iscell(radiusScale)
    [radiusValues,radiusScale]=scaleValues(radiusValues,radiusScale{:});
else
    [radiusValues,radiusScale]=scaleValues(radiusValues,radiusScale);
end

%% Get colors
if isempty(colorValues)
    if size(colorMap,1)~=1
        error('\nplotClusterDisk: empty "colorValues" but colorMap has more than one row (%dx%d)\n',size(colorMap));
    end
    colorBreakpoints=[];
    colorValues=colorMap;
else
    % get values from table
    if isempty(cluster_id)
        if isempty(query)
            t=tbl(:,colorValues);
        else
            t=tbl;  % very wasteful, but query may require other fields
        end
        fun=@(query,colorValues)t.(colorValues)(tableSelect(t,query),:);
        colorValues=callAndStore(fun,query,colorValues);
    else
        if isempty(query)
            t=tbl(:,{cluster_id,colorValues});
        else
            t=tbl;  % very wasteful, but query may require other fields
        end
        fun=@(query,grouping,processing,mfun)tableRowGrouping(t,query,grouping,processing,mfun);
        colorValues=callAndStore(...
            fun,...
            query,...        % Document selection query
            cluster_id,...   % Grouping Variables
            colorValues,...   % Processed Variables
            colorMerging);   % function
    end
    % Quantize color
    if isempty(colorQuantize)
        if any(colorValues~=round(colorValues))
            error('\nplotClusterDisk: empty "colorQuantize" but some color values are not integer-valued.\n');
        end
        if max(colorValues)>size(colorMap,1)
            error('\nplotClusterDisk: empty "colorQuantize" but largest color value %d is larger than the size of "colorMap" (%dx%d)\n',max(colorValues),size(colorMap));
        end
        if min(colorValues)<1
            error('\nplotClusterDisk: empty "colorQuantize" but some color values are smaller than 1\n');
        end
        colorValues=colorMap(colorValues,:);
        colorBreakpoints=[];
    else
        [colorValues,colorBreakpoints]=quantizeValues(colorValues,colorMap,colorQuantize);
    end
end

%% Plot
if all(radiusValuesUnscaled==0)
    radiusValues=min(max(centers)-min(centers))/500;
    bandCharts('centers',centers,...
               'radii',radiusValues,...
               'colors',colorValues,...
               'pointsPerCircle',6,'patchProperties',{'lineStyle','none'});

else
    if solid
        if trimRadius>0
            bandCharts('centers',centers,...
                       'radii',{radiusValues+trimRadius,-trimRadius},...
                       'slicesSharpness',{0,0},...
                       'colors',{trimColor,colorValues},...
                       'pointsPerCircle',pointsPerCircle,'patchProperties',{'lineStyle','none'});
        else
            bandCharts('centers',centers,...
                       'radii',radiusValues,...
                       'slicesSharpness',0,...
                       'colors',colorValues,...
                       'pointsPerCircle',pointsPerCircle,'patchProperties',{'lineStyle','none'});
        end
    else
        bandCharts('centers',centers,...
                   'radii',{radiusValues,-trimRadius},...
                   'slicesSharpness',{0,0},...
                   'colors',{colorValues,trimColor},...
                   'pointsPerCircle',pointsPerCircle,'patchProperties',{'lineStyle','none'});
    end
end

radiusValuesScaled=radiusValues;

fprintf('done plotClusterDisks (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,parameters__);
