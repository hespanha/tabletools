function varargout=bandCharts(varargin);
% To get help, type bandCharts('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha


declareParameter(...
    'Help', {
        'Plots multiple "band charts" as a single patch object.'
        ' '
        'A "band chart" consists of a set of concentric pie charts'
        '(each called a ''band'') with the same center point and the'
        'same angles for each pie slices. However, different bands'
        'may have different radii and different colors.'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','centers',...
    'Description', {
        '#Pies x 2 matrix with the centers of each pie chart.'
        'These centers are shared by all the bands.'
                   });
declareParameter(...
    'VariableName','arc0',...
    'DefaultValue',0,...
    'Description', {
        '1x1 scalar with the angle (in radians) from which arcs are measured.'
        'A zero value for arc0 corresponds to arcs starting from the 3 o''clock'
        'horizontal line. For normal axis, angles are measure counter clockwise'
        'and for ''axis ij'' they are measured clockwise.'
                   });
declareParameter(...
    'VariableName','arcs',...
    'DefaultValue',2*pi,...
    'Description', {
        'Variable that specifies the arcs for the different pies/slices'
        '(shared by all the bands). Should be one of the following:'
        '  1) #pies x #slices matrix - Arcs (in radians) of each slice of each'
        '                              pie chart.'
        '  2) #slices x 1 vector     - Arcs of each slice. This presumes that'
        '                              all pies have similar slices.'
        '  3) 1 x 1 scalar - Single arc that is common to all the pies.'
        '                    This implies #slices=1.'
                   });
declareParameter(...
    'VariableName','radii',...
    'Description', {
        'Cell array of matrices specifying the radii of the different bands.'
        'Each matrix of the cell array must be one of the following:'
        '  1) #pies x #slices matrix - Radius of each slice of each pie chart.'
        '  2) #pies x 1 vector       - Radius of each pie. This presumes that all'
        '                              slices of the same pie have the same radius.'
        '  3) #slices x 1 vector     - Radius of each slice. This presumes that'
        '                              all pies have slices with similar radius.'
        '  4) 1 x 1 scalar           - Common radius for all slices of every pie.'
        'If #pies = #slides and a vector is given, 2) is assumed over 3).'
        ' '
        'Negative values are interpreted as radii relative to the previous band'
        '(inside it since the relative value is negative).'
        ' '
        'When a single band is desired, this function can either take a cell array'
        'with a single element or the matrix itself.'
                   });

declareParameter(...
    'VariableName','slicesSharpness',...
    'DefaultValue',0,...
    'Description', {
        'Cell array of matrices specifying the sharpness of the pie slices'
        'in each band. Each matrix of the cell array must be one of the following:'
        '  1) #pies x #slices matrix - Sharpness of each slice of each pie'
        '  2) #pies x 1 vector       - Sharpness of each pie. This presumes that all'
        '                              slices of the same pie have the same Sharpness.'
        '  3) #slices x 1 vector     - Sharpness of each slice. This presumes that'
        '                              all pies have slices with similar Sharpness.'
        '  4) 1 x 1 scalar           - Common sharpness value for every slice of every pie'
        'Positive sharpness values permit the construction of pie charts that'
        'look like flowers:'
        '  0  - regular pie slice'
        '  .1 - pie slice that looks like a fat flower petal'
        '  10 - pie slice that looks like a long and skinny flower petal'
        'Intermediate values are also allowed, producing different types of petals.'
                   });

declareParameter(...
    'VariableName','colors',...
    'Description', {
        'Cell array of matrices specifying the colors of the different bands.'
        'Each matrix of the cell array must be one of the following:'
        '  1) #pies x #slices x 3 matrix - RGB color of each slice of each pie chart.'
        '  2) #slices x 3 matrix         - RGB color of each slice, it presumes that'
        '                                  the slices of each pie have similar color'
        '  3) #pies x 3 matrix           - RGB color of each pie, it presumes that'
        '                                  all slices of the same pie have the same color'
        '  4) 1 x 3 vector               - Common RGB color for all slices of every pie.'
        'If #pies = #slides and a 2D matrix is given, 2) is assumed over 3).'
        ' '
        'When a single band is desired, this function can either take a cell array'
        'with a single element or the matrix itself.'
                   });

declareParameter(...
    'VariableName','pointsPerCircle',...
    'DefaultValue',36,...
    'Description', {
        'Scalar with the number of points used to draw the curved side of an whole.'
        'circle. A large number will result in a smoother curve, but slower drawing.'
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
    'VariableName','reorder',...
    'AdmissibleValues',{false,true},...
    'DefaultValue',true,...
    'Description', {
        'When ''true'' the circles are ordered so that the band charts with larger radii'
        'are drawn first and will therefore appear to be under bands drwan subsequently.'
        'By the ''radius of a band chart'', we mean the largest radii in  the first band'
        '(defined by radii{1}).'
        'The elements of the cell arrays radii and color are not reordered.'
        'Thus, within each band chart, the bands are drawn starting from radii{1}/color{1}'
        'and ending in radii{end}/color{end}.'
                   });

declareParameter(...
    'VariableName','patchProperties',...
    'DefaultValue',{},...
    'Description', {
        'Cell array of property/value pairs to be passed to the patch command.'
        'These could include:'
        '  ''LineStyle'', ''none'' or ''-'' or ''--'' or '':'' or ''-.'''
        '  ''LineWidth'', <value in pts>'
        'To see all properties use set(patch)'
                   });
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','h',...
    'Description', {
        'handle of the patch object'
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

nPies=size(centers,1);

if prod(size(arcs))==1
    arcs=repmat(arcs,nPies,1);
end
if ~iscell(radii)
    radii={radii};
end
if ~iscell(colors)
    colors={colors};
end

nSlices=size(arcs,2);

if length(radii)~=length(colors)
    error('bandChart: length(radii)=%d does not match length(colors)=%d\n',length(radii),length(colors));
end

nBands=length(radii);
fprintf(' bandChart: %d pies, %d slices, %d bands\n',nPies,nSlices,nBands);

if ~iscell(slicesSharpness)
    fprintf('    guessing: all pies have the same sharpness\n');
    slicesSharpness=repmat({slicesSharpness},nBands,1);
end

t0=clock;

% Completing arcs
if isequal(size(arcs),[nSlices,1]) || isequal(size(arcs),[1,nSlices])
    fprintf('    guessing: all pies have the same arcs\n');
    arcs=repmat(arcs(:)',nPies,1);
end
    
for i=1:nBands
    % Completing radii matrix
    if isequal(size(radii{i}),[nPies,1])
        fprintf('    guessing (band %d): all slices of each pie have same radius (size(radii)=%s',...
                i,mat2str(size(radii{i})));
        radii{i}=repmat(radii{i},1,nSlices);
        fprintf('->%s)\n',mat2str(size(radii{i})));
    elseif isequal(size(radii{i}),[nSlices,1]) || isequal(size(radii{i}),[1,nSlices])
        fprintf('    guessing (band %d): same slice of each pie has same radius (size(radii)=%s',...
                i,mat2str(size(radii{i})));
        radii{i}=repmat(radii{i}(:)',nPies,1);
        fprintf('->%s)\n',mat2str(size(radii{i})));
    elseif prod(size(radii{i}))==1
        fprintf('    guessing (band %d): all slices of all pies have same radius (size(radii)=%s',...
                i,mat2str(size(radii{i})));
        radii{i}=repmat(radii{i},nPies,nSlices);
        fprintf('->%s)\n',mat2str(size(radii{i})));
    end
    % Completing slicesSharpness matrix
    if isequal(size(slicesSharpness{i}),[nPies,1])
        fprintf('    guessing (band %d): all slices of each pie have same slicesSharpness (size(slicesSharpness)=%s',...
                i,mat2str(size(slicesSharpness{i})));
        slicesSharpness{i}=repmat(slicesSharpness{i},1,nSlices);
        fprintf('->%s)\n',mat2str(size(slicesSharpness{i})));
    elseif isequal(size(slicesSharpness{i}),[nSlices,1]) || isequal(size(slicesSharpness{i}),[1,nSlices])
        fprintf('    guessing (band %d): same slice of each pie has same slicesSharpness (size(slicesSharpness)=%s',...
                i,mat2str(size(slicesSharpness{i})));
        slicesSharpness{i}=repmat(slicesSharpness{i}(:)',nPies,1);
        fprintf('->%s)\n',mat2str(size(slicesSharpness{i})));
    elseif prod(size(slicesSharpness{i}))==1
        fprintf('    guessing (band %d): all slices of all pies have same slicesSharpness (size(slicesSharpness)=%s',...
                i,mat2str(size(slicesSharpness{i})));
        slicesSharpness{i}=repmat(slicesSharpness{i},nPies,nSlices);
        fprintf('->%s)\n',mat2str(size(slicesSharpness{i})));
    end
    % Completing color matrix
    if isequal(size(colors{i}),[nSlices,3])
        fprintf('    guessing (band %d): same slice of each pie has same color (size(colors)=%s',...
                i,mat2str(size(colors{i})));
        colors{i}=repmat(reshape(colors{i},1,nSlices,3),nPies,1,1);
        fprintf('->%s)\n',mat2str(size(colors{i})));
    elseif isequal(size(colors{i}),[nPies,3])
        fprintf('    guessing (band %d): all slices of same pie have same color (size(colors)=%s',...
                i,mat2str(size(colors{i})));
        colors{i}=repmat(reshape(colors{i},nPies,1,3),1,nSlices,1);
        fprintf('->%s)\n',mat2str(size(colors{i})));
    elseif isequal(size(colors{i}),[1,3])
        fprintf('    guessing (band %d): all slices of all pies have same color (size(colors)=%s',...
                i,mat2str(size(colors{i})));
        colors{i}=repmat(reshape(colors{i},1,1,3),nPies,nSlices,1);
        fprintf('->%s)\n',mat2str(size(colors{i})));
    end
end

%% Check sizes
if length(size(centers))~=2
    size(centers)
    error('bandChart: centers %d-dimensional array, expected %d-dimensional\n',length(size(centers)),2)
end
if ~isequal(size(centers),[nPies,2])
    error('bandChart: size of centers %dx%d, expected %dx%d\n',size(centers),[nPies,2])
end
if length(size(arcs))~=2
    size(arcs)
    error('bandChart: arcs %d-dimensional array, expected %d-dimensional\n',length(size(arcs)),2)
end
if ~isequal(size(arcs),[nPies,nSlices])
    error('bandChart: size of arcs %dx%d, expected %dx%d\n',size(arcs),[nPies,nSlices])
end
for i=1:length(radii)
    if length(size(radii{i}))~=2
        size(radii{i})
        error('bandChart: radii{%d} %d-dimensional array, expected %d-dimensional\n',i,length(size(radii{i})),2)
    end
    if ~isequal(size(radii{i}),[nPies,nSlices])
        error('bandChart: size of radii{%d} %dx%d, expected %dx%d\n',i,size(radii{i}),[nPies,nSlices])
    end
end
for i=1:length(colors)
    if length(size(colors{i}))~=3
        size(colors{i})
        error('bandChart: colors{%d} %d-dimensional array, expected %d-dimensional\n',i,length(size(colors{i})),3)
    end
    if ~isequal(size(colors{i}),[nPies,nSlices,3])
        error('bandChart: size of colors{%d} %dx%dx%d, expected %dx%dx%d\n',i,size(colors{i}),[nPies,nSlices,3])
    end
end

for i=2:nBands
    % negative radii rare relative to previous one
    k=(radii{i}<0);
    radii{i}(k)=radii{i-1}(k)+radii{i}(k);
end

if reorder
    % order centers by largest outer radius
    [~,k]=sort(max(radii{1},[],2),'descend');
    centers=centers(k,:);
    arcs=arcs(k,:);
    for i=1:nBands
        radii{i}=radii{i}(k,:);
        slicesSharpness{i}=slicesSharpness{i}(k,:);
        colors{i}=colors{i}(k,:,:);
    end
end

% duplicate centers (2nd is a copy of the first all white)
centers2=zeros(nBands*nPies,2);
arcs2=zeros(nBands*nPies,nSlices);
radii2=zeros(nBands*nPies,nSlices);
slicesSharpness2=zeros(nBands*nPies,nSlices);
colors2=zeros(nBands*nPies,nSlices,3);
for i=1:nBands
    centers2(i:nBands:end,:)=centers;
    arcs2(i:nBands:end,:)=arcs;
    radii2(i:nBands:end,:)=radii{i};
    slicesSharpness2(i:nBands:end,:)=slicesSharpness{i};
    colors2(i:nBands:end,:,:)=colors{i};
end

h=pieCharts(centers2,arcs2,arc0,radii2,slicesSharpness2,colors2,...
            halfSlices,pointsPerCircle,patchProperties{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);
