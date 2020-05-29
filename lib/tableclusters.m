classdef tableclusters
% tablecluster class:
%
% Fields:
%  clusters - table with one row per cluster and the following columns:
%             id       - cluster id (integer or categorical)
%             name     - cluster name (string)
%             centroid - array of double characterizing the centrois
%                        (for kmeans is the centroid of cluster,
%                         from SOMs is the codebook value, etc.)
%             position - position of the cluster for graphical representations
%             maxDimension - when equal to an integer 'i', indicates that this  
%                            cluster has the maximum value for 'i'-th dimension
%                            of the centroid, otherwise is equal to NaN
%             maxDimensionWeight - when maxDimension is equal to an integer 'i',
%                            gives the value of the 'i'-the dimension of the centroid,
%                            (which is the largest among all clusters)
%  tbl - table with one row per object being clustered and the following columns:
%             primaryKey - primaryKey of the objects being clustered
%             id         - cluster ID for the object
%             error      - clustering error for that object
%
% tblcluster=tableclusters(primaryKeys,clusteringIDs,clusteringErrors,...
%                          clusterIDs,clusterCentroids,clusterPositions,clusterNames)
% Creates a tablecluster class based on clustering performed on the rows of the table 'tbl'
%  primaryKeys      - n x 1 primaryKeys for the table with clustering results
%  clusteringIDs    - n x 1 array with the result of the clustering
%  clusteringError  - n x 1 array with clustering errors (may be empty)
%  clusterIDs       - c x 1 array with the cluster IDs
%  clusterCentroids - c x ~ array with data pertaining to each cluster (may be empty)
%  clusterPositions - c x 2 or c x 3 array with positions for each cluster (may be empty)
%  clusterNames     - c x 1 array of strings with the names of the clusters (may be empty)
%
% tblcluster=addClusterNamesFromCentroids(tblcluster,centroidColumnNames,threshold)
%   
% Adds names for the clusters based on the centroid vectors.
%  tablecluster        - tablecluster to be processes
%  centroidColumnNames - m x 1 array of strings with the names of each column of the centriods vector
%  threshold           - m x 1 array with thresholds below which the name of the centroid column will be ignored
%
% The names of the clusters are constructed as follows:
%  . all clusters are initially given the empty string '' as their name
%  . for each one of the m columns of the centroid vector
%    . finds which cluster has the maximum value for that column of the centroid vector
%    . if that value is above the corresponding threshold, gives that cluster the corresponding name in centroidColumnNames

% Copyright (C) 2013-2014 Stacy & Joao Hespanha
    
    properties
        clusters    % table with one row per cluster
        tbl         % table with one row per item to be clustered
        
    end
    
    methods (Static)
        
        function tblcluster=loadobj(s)
            tblcluster=tableclusters();
            tblcluster.clusters=s.clusters;
            tblcluster.tbl=s.tbl;
        end
        
    end
    
    methods
        
        function s=saveobj(tblcluster)
            s=struct('clusters',tblcluster.clusters,'tbl',tblcluster.tbl);
        end
        
        function tblcluster=tableclusters(primaryKeys,clusteringIDs,clusteringErrors,...
                                          clusterIDs,clusterCentroids,clusterPositions,clusterNames)
            
            if nargin==0
                return
            end
            if nargin<7 || isempty(clusterNames)
                clusterNames=[];
            end
            if nargin<6
                clusterPositions=[];
            end
            if nargin<5
                clusterCentroids=[];
            end
            
            fprintf('  Creating clustering table... ');
            t1=clock;
            try
                if ~isempty(clusteringErrors)
                    tblcluster.tbl=table(clusteringIDs,clusteringErrors,...
                                         'VariableNames',{'id','error'});
                else
                    tblcluster.tbl=table(clusteringIDs,...
                                         'VariableNames',{'id'});
                end
                tblcluster.tbl.primaryKey=primaryKeys;
            catch me
                fprintf('\nERROR: size(clusteringIDs)=[%d,%d], size(clusteringErrors)=[%d,%d], size(primaryKeys)=[%d,%d]\n',size(clusteringIDs),size(clusteringErrors),size(primaryKeys));
                fprintf('ERROR creating tblcluster.tbl\n');    
                rethrow(me);
            end

            try
                tblcluster.clusters=table(clusterIDs,'VariableNames',{'id'});
                if ~isempty(clusterNames)
                    tblcluster.clusters.name=clusterNames;
                end
                if ~isempty(clusterCentroids)
                    tblcluster.clusters.centroid=clusterCentroids;
                end
                if ~isempty(clusterPositions)
                    tblcluster.clusters.position=clusterPositions;
                end
            catch me
                fprintf('\nERROR: size(clusterIDs)=[%d,%d], size(clusterNames)=[%d,%d], size(clusterCentroids)=[%d,%d], size(clusterPositions)=[%d,%d]\n',size(clusterIDs),size(clusterNames),size(clusterCentroids),size(clusterPositions));
                fprintf('ERROR creating tblcluster.tbl\n');
                rethrow(me);
            end

            if ~isempty(clusterCentroids)
                [weight,maxCluster]=max(tblcluster.clusters.centroid,[],1);
                maxDimension=NaN(size(tblcluster.clusters,1),1);
                maxDimensionWeight=NaN(size(tblcluster.clusters,1),1);
                for j=1:length(weight)
                    if isnan(maxDimension(maxCluster(j))) || weight(j)>maxDimensionWeight(maxCluster(j))
                        maxDimension(maxCluster(j))=j;
                        maxDimensionWeight(maxCluster(j))=weight(j);
                    end
                end
                tblcluster.clusters.maxDimension=maxDimension;
                tblcluster.clusters.maxDimensionWeight=maxDimensionWeight;
            end
            fprintf('done (#items=%d, #clusters=%d, %.3fsec)\n',size(tblcluster.tbl,1),size(tblcluster.clusters,1),etime(clock,t1));
        
        end % function tblcluster=tableclusters()
            
            function tblcluster=addClusterNamesFromCentroids(tblcluster,centroidColumnNames,threshold)
                if nargin<3
                    threshold=-inf;
                end
                
                k=~isnan(tblcluster.clusters.maxDimension) & (tblcluster.clusters.maxDimensionWeight>threshold);
                tblcluster.clusters.name=repmat({''},size(tblcluster.clusters,1),1);
                tblcluster.clusters.name(k)=centroidColumnNames(tblcluster.clusters.maxDimension(k));

            end % tblcluster=addClusterNames
            
    end
    
end