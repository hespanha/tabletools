function h=pieCharts(centers,arcs,arc0,radius,slicesSharpness,colors,...
                     halfSlices,pointsPerCircle,varargin)
% h=pieChars(centers,arcs,arc0,radius,slicesSharpness,colors,...
%                    halfSlices,pointsPerCircle,...)
%
% Plots multiple pie charts as a single patch object.
%
% Inputs:
%   centers (nPies x 2 matrix)       - coordinates of the centers
%                                      of each pie chart 
%   arcs    (nPies x nSlices matrix) - arcs (in radians) of each
%                                      slice of each pie chart 
%   arc0    (scalar)                 - angle (in radians, measured
%                                      counter clockwise) from
%                                      which arcs are measured. 
%                                      When arc0=0, the arcs start
%                                      from the 3 o'clock
%                                      horizontal line. 
%   radius  (nPies x nSlices matrix) - radius of each slice of each pie
%   slicesSharpness (nPies x nSlices matrix) - sharpness of each slice of each pie:
%                                           0  - regular pie slice
%                                           .1 - pie slices look like
%                                                fat flower petals
%                                           10 - pie slices look like 
%                                                long and skinny petals
%   colors  (nPies x nSlices x 3 matrix) - RGB color of each slice
%                                          of each pie 
%   pointsPerCircle  - number of points used to draw the curved side
%                      of an whole circle. (a large number will
%                      result in a smoother curve, but slower
%                      drawing) 
%   halfSlices - scalar specifying the type of slice
%               0 - regular slice
%               1 - half slice (clockwise)
%               -1 - half slice (anti-clockwise)
%   ... (optional)   - optional properties passed directly to the
%                      matlab function 'patch' 
%
% Output:
%   h - handle of the patch object
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

nPies=size(centers,1);

if prod(size(arcs))==1
    arcs=repmat(2*pi,nPies,1);
end

nSlices=size(arcs,2);

fprintf('  pieChart: %d pies, %d slices... ',nPies,nSlices);
t0=clock;

if prod(size(slicesSharpness))==1
    slicesSharpness=slicesSharpness*ones(nPies,nSlices);
end

%% Check sizes
if length(size(centers))~=2
    size(centers)
    error('pieChart: centers %d array, expected %d\n',length(size(centers)),2)
end
if ~isequal(size(centers),[nPies,2])
    error('pieChart: size of centers %dx%d, expected %dx%d\n',size(centers),[nPies,2])
end
if length(size(arcs))~=2
    size(arcs)
    error('pieChart: arcs %d array, expected %d\n',length(size(arcs)),2)
end
if ~isequal(size(arcs),[nPies,nSlices])
    error('pieChart: size of arcs %dx%d, expected %dx%d\n',size(arcs),[nPies,nSlices])
end
if length(size(radius))~=2
    size(radius)
    error('pieChart: radius %d array, expected %d\n',length(size(radius)),2)
end
if ~isequal(size(radius),[nPies,nSlices])
    error('pieChart: size of radradiusius %dx%d, expected %dx%d\n',size(radius),[nPies,nSlices])
end
if length(size(slicesSharpness))~=2
    size(slicesSharpness)
    error('pieChart: slicesSharpness %d array, expected %d\n',length(size(slicesSharpness)),2)
end
if ~isequal(size(slicesSharpness),[nPies,nSlices])
    error('pieChart: size of slicesSharpness %dx%d, expected %dx%d\n',size(slicesSharpness),[nPies,nSlices])
end
if length(size(colors))~=3
    size(colors)
    error('pieChart: colors %d array, expected %d\n',length(size(colors)),3)
end
if ~isequal(size(colors),[nPies,nSlices,3])
    error('pieChart: size of colors %dx%dx%d, expected %dx%dx%d\n',size(colors),[nPies,nSlices,3])
end
if  prod(size(halfSlices))~=1
    size(halfSlices)
    error('pieChart: scalar expected in halfSlices, matrix with %d entries found instead\n',...
          prod(size(halfSlices)));
end
if  prod(size(pointsPerCircle))~=1
    size(pointsPerCircle)
    error('pieChart: scalar expected in pointsPerCircle, matrix with %d entries found instead\n',...
          prod(size(pointsPerCircle)));
end

%% compute wedges arcs 
cumArcs=cumsum([repmat(arc0,nPies,1),arcs],2);  % nPies x (nSlices+1)
wedgeStarts=cumArcs(:,1:end-1);                 % nPies x nSlices
wedgeEnds=cumArcs(:,2:end);                     % nPies x nSlices
wedgeZ=-0.002*repmat((nPies:-1:1)',1,nSlices);  % nPies x nSlices

%% stack all wedges (one per column) 

nWedges=nPies*nSlices;

% (stacking all wedges of the same pie together)
wedgeWidth=reshape(arcs',1,nWedges);          % 1 x nWedges
wedgeStarts=reshape(wedgeStarts',1,nWedges);  % 1 x nWedges
wedgeEnds=reshape(wedgeEnds',1,nWedges);      % 1 x nWedges
wedgeRadius=reshape(radius',1,nWedges);       % 1 x nWedges
wedgeSliceSharpness=reshape(slicesSharpness',1,nWedges);         % 1 x nWedges
wedgeCentersX=reshape(repmat(centers(:,1),1,nSlices)',1,nWedges);% 1 x nWedges
wedgeCentersY=reshape(repmat(centers(:,2),1,nSlices)',1,nWedges);% 1 x nWedges
wedgeZ=reshape(wedgeZ',1,nWedges);            % 1 x nWedges
wedgeColors=zeros(1,nWedges,3);               % 1 x nWedges x 3
wedgeColors(:,:,1)=reshape(colors(:,:,1)',1,nWedges,1);
wedgeColors(:,:,2)=reshape(colors(:,:,2)',1,nWedges,1);
wedgeColors(:,:,3)=reshape(colors(:,:,3)',1,nWedges,1);

%% remove empty edges
k=(wedgeStarts==wedgeEnds | wedgeRadius==0);
if sum(k)>0
    fprintf('removing %d/%d empty edges... ',sum(k),length(k));
    wedgeStarts(:,k)=[];
    wedgeEnds(:,k)=[];
    wedgeWidth(:,k)=[];
    wedgeRadius(:,k)=[];
    wedgeSliceSharpness(:,k)=[];
    wedgeCentersX(:,k)=[];
    wedgeCentersY(:,k)=[];
    wedgeZ(:,k)=[];
    wedgeColors(:,k,:)=[];
    nWedges=size(wedgeStarts,2);
end


%[wedgeStarts;wedgeEnds;wedgeRadius;reshape(wedgeColors,nWedges,3)']'

%% Computes angles for each point in the wedge boundary
if 0
    %% fixed angle between points
    % angles for each wedge (starting at zero)
    angles=linspace(0,2*pi,pointsPerCircle);
    angles=repmat(angles',1,size(wedgeStarts,2)); % pointsPerCircle x nWedges
    % collaps all angles exceeding wedgeWidth
    angles=min(angles,repmat(wedgeWidth,pointsPerCircle,1));
else
    %% fixed # of points per edge
    angles=linspace(0,1,pointsPerCircle);
    angles=angles'*wedgeWidth;
end
% cumulative angles for each wedge (starting at wedge Start)
cumAngles=repmat(wedgeStarts,pointsPerCircle,1)+angles; % pointsPerCircle x nWedges

%% coordinates of the points w.r.t. center of pie chart
wedgeRadius=repmat(wedgeRadius,pointsPerCircle,1)...
    .*abs(sin(pi*angles./repmat(wedgeWidth,pointsPerCircle,1)))...
    .^repmat(wedgeSliceSharpness,pointsPerCircle,1);
if halfSlices~=0
    % find angles closest to mid-point
    [~,k]=min(abs(angles-repmat(wedgeWidth,pointsPerCircle,1)/2),[],1);
    % clear angles below/above
    if halfSlices>0
        kk=repmat((1:size(angles,1))',1,nWedges);
    else
        kk=repmat((size(angles,1):-1:1)',1,nWedges);
    end
    i=kk>repmat(k,size(kk,1),1);
    wedgeRadius(i)=0;
end    
% if halfSlices>0
%     k=find(angles<repmat(wedgeWidth,pointsPerCircle,1)/2-pi/180);
%     wedgeRadius(k)=0;
% elseif halfSlices<0
%     k=find(angles>repmat(wedgeWidth,pointsPerCircle,1)/2+pi/180);
%     wedgeRadius(k)=0;
% end

%% remove rows full of empty radius
k=all(abs(wedgeRadius)<1e-5*max(wedgeRadius(:)),2);
wedgeRadius(k,:)=[];
cumAngles(k,:)=[];

wedgesX=[zeros(1,nWedges);wedgeRadius.*cos(cumAngles);zeros(1,nWedges)];
wedgesY=[zeros(1,nWedges);wedgeRadius.*sin(cumAngles);zeros(1,nWedges)];

%% add center coordinates
wedgesX=wedgesX+repmat(wedgeCentersX,size(wedgesX,1),1);
wedgesY=wedgesY+repmat(wedgeCentersY,size(wedgesX,1),1);
wedgeZ=repmat(wedgeZ,size(wedgesX,1),1);

wedgeColors=repmat(wedgeColors,size(wedgesX,1),1,1);

% Scale by radius and add center

%% plot
if size(wedgesX,2)>60000
    fprintf('TOO MANY PATCHES (%d>60000), expect trouble with .eps file\n',...
            size(wedgesX,2));
end

if 0
    % rely on z for correct order -- problems with multiple pieCharts
    h=patch(wedgesX,wedgesY,wedgeZ,wedgeColors,...
            'faceColor','flat','marker','none','linestyle','none',varargin{:});
else
    % rely on order of plot
    h=patch(wedgesX,wedgesY,wedgeColors,...
            'faceColor','flat','marker','none','linestyle','none',varargin{:});
    %set(gca,'DrawMode','fast'); MATLAB 2013
    set(gca,'SortMethod','child'); % MATLAB 2014
end
%view(2);
fprintf('done (%d points, %d patches, %.3f)\n',...
        prod(size(wedgesX)),size(wedgesX,2),etime(clock,t0));