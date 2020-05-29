function varargout=tableLinkageColumnClustering(varargin);
% To get help, type tableLinkageClustering('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This script clusters takes a table variable (i.e., column) that'
        'contains a matrix and performs linkage clustering on the columns'
        'of such matrix.'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','tbl',...
    'Description', {
        'Table to be processed.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','tableColumn',...
    'Description', {
        'Name of the table column that contains the features used by the'
        'linkage classifier. This column should contain a matrix of numbers,'
        'whose columns are to be clustered based on the entries of the rows.'
        'In particular, tbl.(tableColumn) should result in a matrix,'
        'with one vector/features to be clustered per column.'
                   });
    declareParameter(...
    'VariableName','nClusters',...
    'DefaultValue',8,...
    'Description', {
        'Desired number of clusters. A vector can be provided to perform clustering'
        'with several different number of clusters.'
                   });
declareParameter(...
    'VariableName','clusteringDistance',...
    'DefaultValue','cosine',...
    'AdmissibleValues',{ 'euclidean','seuclidean','cityblock','minkowski','chebychev','mahalanobis','cosine','correlation','spearman','hamming','jaccard'},... % for linkage
    'Description', {
        'This variable specifices the metric used for clustering. See ''help linkage''.'
                   });
declareParameter(...
    'VariableName','seed',...
    'DefaultValue',0,...
    'Description', {
        'Seed used to initialize the random number generator.'
                   });


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','dendogramPrefix',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for the .eps output files, used for write access (output)'
        'Prefix for the name of the files containing images of the dendograms.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','tblcluster',...
    'Description', {
        'Result of the clustering, in the form of a tablecluster object.'
        'When nClusters is a vector, a cell array of tables is returned.'
                   });
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,params]=setParameters(nargout,varargin);
if stopNow
    return 
end

%verboseLevel=4;

streamSelect=RandStream.create('mt19937ar','seed',seed);
RandStream.setGlobalStream(streamSelect);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('tableLinkageClustering: starting\n');
t0=clock;


dataMatrix=full(tbl.(tableColumn))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Remove documents with NaN values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kNaN=find(all(~isnan(dataMatrix),2));

if size(dataMatrix,1)>length(kNaN)
    fprintf('  %d out of %d columns have NaN counts - will be ignored\n',...
            size(dataMatrix,1)-length(kNaN),size(dataMatrix,1));
    dataMatrix=dataMatrix(kNaN,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% kmeans clustering (with given distance)
fprintf('  computing tree (distance ''%s'', %d observations to be clustered, %d data variabels)... ',clusteringDistance,size(dataMatrix,1),size(dataMatrix,2));
t1=clock;
tree=linkage(dataMatrix,'average',clusteringDistance);
fprintf('done (%.3f sec)\n',etime(clock,t1));

fprintf('  clustering: ');
fprintf('%d... ',nClusters);
t1=clock;
clusteringIDs=cluster(tree,'MaxClust',nClusters);
fprintf('done (%.3f sec)\n',etime(clock,t1));


fprintf('  constructing clustering tables:\n');
t1=clock;
nClusters=zeros(size(clusteringIDs,2),1);
primaryKey=(1:size(dataMatrix,1))';
for i=1:size(clusteringIDs,2)
    nClusters(i)=length(unique(clusteringIDs(:,i)));
    clusterIDs=(1:nClusters(i))';
    if size(clusteringIDs,2)>1
        tblcluster{i}=tableclusters(primaryKey,clusteringIDs(:,i),[],...
                                    clusterIDs,[],[]);
    else
        tblcluster=tableclusters(primaryKey,clusteringIDs(:,i),[],...
                                 clusterIDs,[],[]);
    end
end
fprintf('      done (%.3f sec)\n',etime(clock,t1));


% regular dendrogram
fprintf('  dendogram...',nClusters);
t1=clock;
clearFigure('orientation','landscape','figureNumber',1,'figureName','dendogram');
dendrogram(tree,0);
fprintf('done  (%.3f sec)\n',etime(clock,t1));
fontsize=ceil(min(10,1600/size(clusteringIDs,1)));
fprintf('print dendrogram (fontsize = %g): %s... ',fontsize,[dendogramPrefix,'+dendrogram']);
t1=clock;
xticklabels=cellstr(get(gca,'xticklabel'));
format_ticks(gca,xticklabels,[],... % x/y labels
             [],[],... % x/y tick positions
             90,[],... % x/y label rotations
             .5,...  % label offset
             'FontSize',fontsize,...
             'FontName','Helvetica');
saveFigure('filename',[dendogramPrefix,'+dendrogram'],'saveFormat','all');
fprintf('done  (%.3f sec)\n',etime(clock,t1));

% polar dendrogram
fprintf('  polardendogram...',nClusters);
t1=clock;
clearFigure('orientation','landscape','figureNumber',2,'figureName','polardendogram');
fontsize=ceil(min(10,1600/size(clusteringIDs,1)));
fprintf('print polardendrogram (fontsize = %g): %s... ',fontsize,[dendogramPrefix,'+polardendrogram']);
set(gca,'FontSize',fontsize,...
        'FontName','Helvetica');
polardendrogram(tree,0);
saveFigure('filename',[dendogramPrefix,'+polardendrogram'],'saveFormat','all');
fprintf('done  (%.3f sec)\n',etime(clock,t1));

clearFigure('orientation','landscape','figureNumber',2,'figureName',sprintf('dendogram (%d clusters)',nClusters(end)));
old=get(0,'RecursionLimit');
set(0,'RecursionLimit',1000);
dendrogram(tree,30,'colorthreshold','default','label',num2str(clusteringIDs(:,end)))
set(0,'RecursionLimit',old);
drawnow

fprintf('done tableKmeansClustering (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);



