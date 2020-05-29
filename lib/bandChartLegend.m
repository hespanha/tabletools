function varargout=bandChartLegend(varargin);
% To get help, type bandCharts('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'Plots multiple "band charts" as a single patch object.'
        ' '
        'A "band chart" consists of a set of concentric pie charts'
        '(each called a "band") with the same center point and the'
        'same angles for each pie slices. However, different bands'
        'may have different radii and different colors.'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','chartAxis',...
    'Description', {
        'Axis of the chart for which the legend is being created.'
                   });
declareParameter(...
    'VariableName','legendAxis',...
    'Description', {
        'Axis where the legend should be created.'
                   });
declareParameter(...
    'VariableName','legendType',...
    'AdmissibleValues',{'colored-disks','colored-slices','sized-disks','sized-slices'},...
    'Description', {
        'Type of legend.'
                   });
declareParameter(...
    'VariableName','orientation',...
    'DefaultValue','vertical',...
    'AdmissibleValues',{'horizontal','vertical'},...
    'Description', {
        'Orientation along which the different elements will be drawn.'
                   });
declareParameter(...
    'VariableName','colors',...
    'DefaultValue',[0,0,0],...
    'Description', {
        'Matrix with the RGB colors for the legend elements. Should be:'
        '1) n x 3 for "colored-disks" or "colored-slices" legends'
        '         (n elements will be drawn, each with a different color)'
        '2) 1 x 3 for "sized-disks" or "sized-slices" legends'
        '         (all elements will have the same color)'
                   });
declareParameter(...
    'VariableName','radii',...
    'DefaultValue',[],...
    'Description', {
        'Vector with the radii for the legend elements. Should be:'
        '   n x 1 for "sized-disks" or "sized-slices" legends'
        '         (n elements will be drawn, each with a different radius)'
        'Ignored for "colored-disks" or "colored-slices" legends'
        '(all elements will have an automatically selected radius)'
                   });
declareParameter(...
    'VariableName','labels',...
    'Description', {
        'Can be one of three options:'
        '1) n x 1 cell array with the labels for the different elements.'
        '2) (n-1) x 1 numerical array of breakpoints, from which labels'
        '             will be constructed'
        '3) string "auto" This case should only occur with for the'
        'legendTypes "sized-disks" and "sized-slices". In this case, '
        'kmeans() is used to compute "nElements" representative radii based'
        'on the ones provided in "radii." These representatives are then used'
        'to select the radii of the legend elements and the corresponding'
        'legends.'
                   });
declareParameter(...
    'VariableName','nElements',...
    'DefaultValue',3,...
    'Description', {
        'Number of legend elements when labels="auto".'
                   });
declareParameter(...
    'VariableName','scale',...
    'DefaultValue',struct('a',1,'b',0,'p',1),...
    'Description', {
        'Structure with the scaling parameters used by scaleValues()'
        'to scale the values provided in radii when labels="auto".'
                   });
declareParameter(...
    'VariableName','arc0',...
    'DefaultValue',0,...
    'Description', {
        '1x1 scalar with the angle (in radians) from which arcs are measured.'
        'A zero value for arc0 corresponds to arcs starting from the 3 o''clock'
        'horizontal line. For normal axis, angles are measure counter clockwise'
        'and for "axis ij" they are measured clockwise.'
        'This parameter is only used for "colored-disks" or "sized-disks" legends.'
                   });
declareParameter(...
    'VariableName','pointsPerCircle',...
    'DefaultValue',36,...
    'Description', {
        'Scalar with the number of points used to draw the curved side of an whole.'
        'circle. A large number will result in a smoother curve, but slower drawing.'
                   });
declareParameter(...
    'VariableName','patchProperties',...
    'DefaultValue',{},...
    'Description', {
        'Cell array of property/value pairs to be passed to the patch command.'
        'These could include:'
        '  "LineStyle", "none" or "-" or "--" or ":" or "-."'
        '  "LineWidth", <value in pts>'
        'To see all properties use set(patch)'
                   });
declareParameter(...
    'VariableName','textProperties',...
    'DefaultValue',{},...
    'Description', {
        'Cell array of property/value pairs to be passed to the text command.'
        'These could include:'
        '  "FontSize", <value in pts>'
        '  "EdgeColor", [red,green,blue] or "k", "b", etc.'
        '  "LineWidth", <value in pts>'
        'To see all properties use set(line)'
                   });
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','hpatch',...
    'Description', {
        'handle of the patch object'
                   });
declareOutput(...
    'VariableName','htext',...
    'Description', {
        'handle of the label object'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters and inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,parameters]=setParameters(nargout,varargin);
if stopNow
    return 
end

%verboseLevel=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Duplicates the original axis to preserve sizes
%set(chartAxis,'Units','in','camerapositionmode','manual');
% fields={'Units','Projection',...
%         'Box',...
%         'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle',...
%         'XDir','YDir','ZDir','XLim','YLim','ZLim','CLim','XScale','YScale','ZScale',...
%         'Position','XAxisLocation','YAxisLocation','Clipping','Visible'};
fields={'Units','Projection',...
        'Box',...
        'XDir','YDir','ZDir','XLim','YLim','ZLim','CLim','XScale','YScale','ZScale',...
        'Position','XAxisLocation','YAxisLocation','Clipping','Visible'};
for i=1:length(fields)
    set(legendAxis,fields{i},get(chartAxis,fields{i}));
end
view(2);

if size(radii,1)<nElements
    nElements=size(radii,1);
end

if ischar(labels) && strcmp(labels,'auto') && ismember(legendType,{'sized-disks','sized-slices'})
    mn=min(radii);mx=max(radii);
    breakpoints=(mn+(mx-mn)*(1:nElements)/(nElements+1));
    % fprintf('Before: ')
    % disp(breakpoints)
    s=warning('off','stats:kmeans:EmptyCluster');
    [idx,c]=kmeans(radii,nElements,'start',breakpoints','emptyaction','singleton',...
                   'replicates',1,'MaxIter',500,'OnlinePhase','off');
    warning(s);
    %radii=scaleValues(sort(c(:)),scale);
    radii=scaleValues(linspace(min(c),max(c),nElements)',scale);
    labels=cellstr(num2str(c,'%.0f'));
else
    radii=scaleValues(radii,scale);
end

nElements = length(labels);

if isnumeric(labels)
    ls=labels;
    labels=cell(nElements+1,1);
    labels{1}=sprintf('<= %.2g',ls(1));
    for i=2:nElements
        labels{i}=sprintf('(%.2g,%.2g]',ls(i-1),ls(i));
    end
    nElements=nElements+1;
    labels{nElements}=sprintf('> %.2g',ls(end));
end


%% Compute locations for objects 
limits=axis;

tol=1.2;

switch (orientation)
  case 'horizontal'
    % horizontally spaced
    minX=limits(1);
    maxX=limits(2);
    minY=(limits(3)+limits(4))/2;
    maxY=(limits(3)+limits(4))/2;
    maxRadius=(maxX-minX)/nElements/2/tol;
    
    dX=0;
    dY=maxRadius;
    dR=0;
    if strcmp(get(legendAxis,'Ydir'),'reverse')
        dR=pi;
        dY=-dY;
    end
    % if get(legendAxis,'CameraUpVector')*[0;1;0]>0
    %     dY=-dY;
    % end
    textProp={'HorizontalAlignment','center','VerticalAlignment','top'};
  case 'vertical'
    % vertically spaced
    minX=(limits(1)+limits(2))/2;
    maxX=(limits(1)+limits(2))/2;
    minY=limits(3);
    maxY=limits(4);
    maxRadius=(maxY-minY)/nElements/2/tol;

    dX=maxRadius;
    dY=0;
    dR=0;
    if strcmp(get(legendAxis,'Xdir'),'reverse')
        fprintf('reversed xaxis\n')
        dR=pi;
        dX=-dX;
    end
    if get(legendAxis,'CameraUpVector')*[0;1;0]<0
        fprintf('camera up [%d,%d,%d]\n',get(legendAxis,'CameraUpVector'))
        dR=pi;
        dX=-dX;
    end
    textProp={'HorizontalAlignment','left','VerticalAlignment','middle'};
  otherwise
    error('bandChartLegend: unkown orientation "%s"\n',orientation);
end

if ismember(legendType,{'colored-disks','colored-slices'})
    radii=maxRadius*ones(nElements,1);
end

% Compute centers
centersX=tol*radii(1)*ones(nElements,1);
for i=2:nElements
    centersX(i)=centersX(i-1)+tol*(radii(i-1)+radii(i));
end
centersY=centersX;
centersX=minX+centersX*(maxX-minX)/(centersX(end)+tol*radii(end));
centersY=minY+centersY*(maxY-minY)/(centersY(end)+tol*radii(end));

switch (legendType)
  case {'colored-disks'},
    %% Colored disks;
    hpatch=bandCharts('centers',[centersX,centersY],...
                      'radii',radii,'colors',colors,...
                      'arc0',pi/2,'pointsPerCircle',pointsPerCircle,'patchProperties',patchProperties);
    htext=labelChart('centers',[centersX+1.5*dX,centersY+1.5*dY],...
                     'labels',labels,...
                     'textProperties',[textProp,textProperties]);
    
  case 'colored-slices',
    %% Colored slices
    hpatch=bandCharts('centers',[centersX,centersY],...
                      'arc0',pi/4+dR,'arcs',3*pi/2,...
                      'radii',radii,'colors',colors,...
                      'pointsPerCircle',pointsPerCircle,'patchProperties',patchProperties);
    htext=labelChart('centers',[centersX+1.5*dX,centersY+1.5*dY],...
                     'labels',labels,...
                     'textProperties',[textProp,textProperties]);
    
  case {'sized-disks'}
    hpatch=bandCharts('centers',[centersX,centersY],...
                      'radii',radii,'colors',colors,...
                      'arc0',pi/2,'pointsPerCircle',pointsPerCircle,'patchProperties',patchProperties);
    htext=labelChart('centers',[centersX+1.5*dX,centersY+1.5*dY],...
                     'labels',labels,...
                     'textProperties',[textProp,textProperties]);
    
  case {'sized-slices'}
    hpatch=bandCharts('centers',[centersX,centersY],...
                      'arc0',pi/4+dR,'arcs',3*pi/2,...
                      'radii',radii,'colors',colors,...
                      'pointsPerCircle',pointsPerCircle,'patchProperties',patchProperties);
    htext=labelChart('centers',[centersX+1.5*dX,centersY+1.5*dY],...
                     'labels',labels,...
                     'textProperties',[textProp,textProperties]);
  otherwise
    error('bandChartLegend: unkown legendType "%s"\n',legendType);
end
