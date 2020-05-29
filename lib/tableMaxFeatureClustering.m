function varargout=tableMaxFeatureClustering(varargin);
% To get help, type tableMaxFeatureClustering('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This script clusters the rows of a table in as many classes as'
        'the elements of the feature vector. Each row is assigned to a class'
        'based on which component of the feature vector is largest.'
        'In particular, a row is assigned to class c, if the c-th entry'
        'of the feature vector for that row is the largest.'
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
        'classifier. This column should contain rows of numbers.'
        'In particular, tbl.(tableColumn) should result in a matrix,'
        'with one feature vector per row and one feature per column.'
                   });
declareParameter(...
    'VariableName','featureNames',...
    'DefaultValue',[],...
    'Description', {
        'String array with the names of the features used by the classifier.'
        'This vector should have as many elements as the number of columns'
        'of the matrix tbl.(tableColumn).'
        'These names are used to name the clusters.'
        'In particular, the cluster corresponding to the rows whose'
        'i-th entry is the largest, will be named featureMames{i}.'
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('tableMaxFeatureClustering: starting\n');
t0=clock;

dataMatrix=full(tbl.(tableColumn));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Remove rows with NaN values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kNaN=find(all(~isnan(dataMatrix),2));

if size(dataMatrix,1)>length(kNaN)
    fprintf('  %d out of %d rows have NaN counts - will be ignored\n',...
            size(dataMatrix,1)-length(kNaN),size(dataMatrix,1));
    dataMatrix=dataMatrix(kNaN,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Do clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,clusteringIDs]=max(dataMatrix,[],2);

clusterIDs=(1:size(dataMatrix,2))';
tblcluster=tableclusters(tbl.primaryKey,clusteringIDs,[],...
                         clusterIDs,[],[],featureNames);

fprintf('done tableMaxFeatureClustering (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);



