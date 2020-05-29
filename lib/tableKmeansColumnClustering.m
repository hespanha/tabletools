function varargout=tableKmeansColumnClustering(varargin);
% To get help, type tableKmeansClustering('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This script clusters takes a table variable (i.e., column) that'
        'contains a matrix and performs k-means clustering on the columns'
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
        'k-means classifier. This column should contain a matrix of numbers,'
        'whose columns are to be clustered based on the entries of the rows.'
        'In particular, tbl.(tableColumn) should result in a matrix,'
        'with one vector/features to be clustered per column.'
                   });
    declareParameter(...
    'VariableName','nClusters',...
    'DefaultValue',8,...
    'Description', {
        'Desired number of clusters.'
                   });
declareParameter(...
    'VariableName','clusteringDistance',...
    'DefaultValue','cosine',...
    'AdmissibleValues',{'sqEuclidean','cityblock','cosine','correlation','Hamming'},... % for kmeans
    'Description', {
        'This variable specifices the metric used for clustering. See ''help kmeans''.'
                   });
declareParameter(...
    'VariableName','replicates',...
    'DefaultValue',10,...
    'Description', {
        'This variable specifies the number of replicates, i.e., the'
        'number of times to repeat the clustering, each with a'
        'new set of initial centroids. See ''help kmeans.'''
                   });
declareParameter(...
    'VariableName','seed',...
    'DefaultValue',0,...
    'Description', {
        'Seed used to initialize the random number generator.'
                   });


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','tblcluster',...
    'Description', {
        'Result of the clustering, in the form of a tablecluster object.'
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

fprintf('tableKmeansClustering: starting\n');
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
fprintf('  kmeans: # of clusters = %d\n',nClusters);
if 0
    [clusteringIDs,centroids]=litekmeans(dataMatrix',nClusters);
    [clusteringIDs,clusterCentroids,dummy,distances]=kmeans(...
        dataMatrix,nClusters,...
        'distance',clusteringDistance,'replicates',1,...
        'EmptyAction','singleton',...
        'start',centroids',...
        'MaxIter',200,'display','final');
    %'MaxIter',200,'display','iter');
else
    %[clusteringIDs,centroids]=litekmeans(dataMatrix',nClusters);
    [clusteringIDs,clusterCentroids,dummy,distances]=kmeans(...
        dataMatrix,nClusters,...
        'distance',clusteringDistance,'replicates',replicates,...
        'EmptyAction','singleton',...
        'start','sample',...
        'MaxIter',200,'display','final');
    %'MaxIter',200,'display','iter');
end
if 0
    clusteringQerror=distances(sub2ind(size(distances),(1:length(clusteringIDs))',clusteringIDs));
    if any(clusteringQerror~=min(distances,[],2))
        error('kmeans classification does not match minimum distance\n');
    end
else
    clusteringQerror=min(distances,[],2);
end

clusterIDs=(1:nClusters)';
primaryKey=(1:size(dataMatrix,1))';
tblcluster=tableclusters(primaryKey,clusteringIDs,clusteringQerror,...
                         clusterIDs,clusterCentroids,[]);

fprintf('done tableKmeansClustering (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);



