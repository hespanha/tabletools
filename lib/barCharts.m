function h=barCharts(centers,heights,colors,barWidth,barSeparation,edgeColor,lineStyle,lineWidth)
% h=barChars(centers,heights,colors,barWidth,barSeparation,edgeColor,lineStyle,lineWidth)
%
% Plots multiple bar charts as a single patch object.
%
% Inputs:
%   centers (nCharts x 2 matrix)     - coordinates of the centers of each bar chart
%   heights (nCharts x nBars matrix) - heights of each bar of each chart
%   colors  (nCharts x nBars x 3 matrix) - RGB color of each bar of each chart
%                                      If a 3x1 vector is provided, all bars are assumed
%                                      to have the same color.
%                                      If a nBarsx3 vector is provided, all charts are assumed
%                                      to have the same colors for their bars
%   edgeColor - edgeColor for the bar edges (see set(patch) for options)
%   lineStyle - lineStyle for the bar edges (see set(patch) for options)
%   lineWidth - lineWidth for the bar edges (see set(patch) for options)
% Output:
%   h - handle of the patch object
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

[nCharts,nBars]=size(heights);

fprintf('  barChart: %d charts, %d bars... ',nCharts,nBars);
t0=clock;

if isequal(size(colors),[nBars,3])
    colors=repmat(reshape(colors,1,nBars,3),nCharts,1,1);
end

if isequal(size(colors),[1,3])
    colors=repmat(reshape(colors,1,1,3),nCharts,nBars,1);
end

%% compute bar's coordinates
barXs=repmat([barWidth,barSeparation],1,nBars);
barXs=[0,cumsum(barXs)];
barStarts=repmat(barXs(1:2:end-1),nCharts,1);
barEnds=repmat(barXs(2:2:end),nCharts,1);
% center horizontaly
barStarts=barStarts-barEnds(end)/2;
barEnds=barEnds-barEnds(end)/2;

%% stack all bars
barStarts=reshape(barStarts,1,nCharts*nBars);
barEnds=reshape(barEnds,1,nCharts*nBars);
barHeights=reshape(heights,1,nCharts*nBars);
barCentersX=reshape(repmat(centers(:,1),1,nBars),1,nCharts*nBars);
barCentersY=reshape(repmat(centers(:,2),1,nBars),1,nCharts*nBars);
barColors=reshape(colors,1,nCharts*nBars,3);
% center vertically
mx=max(max(barHeights(:)),0);
mn=min(min(barHeights(:)),0);
barZero=-(mx+mn)/2;
barTops=barHeights+barZero;

%[barStarts;barEnds;barHeights;reshape(barColors,nCharts*nBars,3)']'

%% Compute points for each bar
barXs=[barStarts;barEnds;barEnds;barStarts;barStarts];
barYs=[barZero*ones(2,nCharts*nBars);barTops;barTops;barZero*ones(1,nCharts*nBars)];
barCentersX=repmat(barCentersX,5,1);
barCentersY=repmat(barCentersY,5,1);
barColors=repmat(barColors,5,1,1);

% Add center
barXs=barXs+barCentersX;
barYs=barYs+barCentersY;

%% plot
if size(barXs,2)>60000
    fprintf('TOO MANY PATCHES (%d>60000), expect trouble with .eps file\n',...
            size(barXs,2));
end
h=patch(barXs,barYs,barColors,...
        'lineStyle',lineStyle,'lineWidth',lineWidth,'edgeColor',edgeColor,...
        'marker','none');

axis equal

fprintf('done (%d points, %d pacthes, %.3f)\n',prod(size(barXs)),size(barXs,2),etime(clock,t0));