function varargout=tableTSNE(varargin)
% To get help, type tableTSNE('help')
%
% Copyright (C) 2014 Joao Hespanha & Stach Rebich Hespanha

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
        'This script uses the Diffusion Stochastic Neighbor Embedding'
        'algorithm (t-SNE) to perform dimensionality reduction on the'
        'values of a column of table.'
        'This script is essentially a wrapper that calls the fast'
        'Intel (IPP) implementation of t-SNE.'
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
    'VariableName','structField',...
    'DefaultValue','',...
    'Description', {
        'When this parameter is nonempty, it is assumed that ''tbl'' is actually'
        'a structure instead of a table. In this case, the table of interest is';
        'a field in the structure with name ''structField'', i.e., the table of';
        'insterest is:'
        '  tbl.(structField)';
        'Moreover, the output will again be a structure that is similar to the input';
        'structure, but with the output table in the same field ''structField''.';
        'This will be useful, e.g., to ran tableTSNE() on a tablecluster obtained';
        'from tableSOMClustering() or tableKmeansClustering().'
                   });

declareParameter(...
    'VariableName','tableHDColumn',...
    'Description', {
        'Name of the table''s column that contains the (high dimensional)'
        'data points. This column should contain rows of numbers.'
        'In particular, tbl.(tableColumn) should result in a matrix,'
        'with one data point vector per row.'
                   });

declareParameter(...
    'VariableName','tableLDColumn',...
    'DefaultValue','tSNEvector',...
    'Description', {
        'Name of the table''s column where the (low dimensional)'
        'data points will be stored.'
                   });

declareParameter(...
    'VariableName','tableCostColumn',...
    'DefaultValue','tSNEcost',...
    'Description', {
        'Name of the table''s column where the data point costs'
        'will be stored.'
                   });

declareParameter(...
    'VariableName','topicWeightThreshold',...
    'DefaultValue',0,...
    'Description', {
        'This variable specifies the topic weight below which topic proportions'
        'should be treated as 0.'
                   });

declareParameter(...
    'VariableName','landmarks',...
    'DefaultValue',[],...
    'Description', {
        'The percentage of points to use as landmarks, in the interval'
        'interval [0,1].'
        'One may also provide a vector with the indices of the rows'
        'to be used as landmarks.'
        'When empty [], the landmark selection is done automatically.'
        'Not used by Barnes-Hut t-SNE.'
                   });

declareParameter(...
    'VariableName','perplexity',...
    'DefaultValue',30,...
    'Description', {
        'Preplexity used in the Gaussian kernal. The perplexity is'
        'mainly of influence on small datasets where landmarks are'
        'not employed.'
        });

declareParameter(...
    'VariableName','BarnesHut',...
    'AdmissibleValues',{false,true},...
    'DefaultValue',true,...
    'Description', {
        'When ''true'' the Barnes-Hut t-SNE algorithm is used.'
        'This algorithm was designed to run on large (N > 5000)'
        'data sets. It  may give poor performance on very small'
        'data sets (it is better to use the standard t-SNE'
        'implementation on such data).'
        });
    
declareParameter(...
    'VariableName','theta',...
    'DefaultValue',.5,...
    'Description', {
        'This variable sets the trade-off parameter between speed'
        'and accuracy: theta = 0 corresponds to standard, slow t-SNE,'
        'while theta = 1 makes very crude approximations.'
        'Appropriate values for theta are between 0.1 and 0.7'
        'Only used by Barnes-Hut t-SNE.'
                   });

declareParameter(...
    'VariableName','nIterations',...
    'DefaultValue',1000,...
    'Description', {
        'Total number of iterations.'
        'Only used by Barnes-Hut t-SNE.'
                   });

declareParameter(...
    'VariableName','pValueIteration',...
    'DefaultValue',250,...
    'Description', {
        'Number of iteration at which the p-values starts to be computed.'
        'Generally equal to ''momentumSwitchIteration''.'
        'Only used by Barnes-Hut t-SNE.'
                   });

declareParameter(...
    'VariableName','initialMomentum',...
    'DefaultValue',.5,...
    'Description', {
        'Initial value of the momentum for the gradient update.'
        'Only used by Barnes-Hut t-SNE.'
                   });

declareParameter(...
    'VariableName','finalMomentum',...
    'DefaultValue',.8,...
    'Description', {
        'Final value of the momentum for the gradient update.'
        'Only used by Barnes-Hut t-SNE.'
                   });

declareParameter(...
    'VariableName','momentumSwitchIteration',...
    'DefaultValue',250,...
    'Description', {
        'Number of iteration at which the momentum for the gradient update'
        'switched from ''initialMomentum'' to ''finalMomentum''.'
        'Only used by Barnes-Hut t-SNE.'
                   });

declareParameter(...
    'VariableName','eta',...
    'DefaultValue',200,...
    'Description', {
        'Gain for the gradient update.'
        'Only used by Barnes-Hut t-SNE.'
                   });

declareParameter(...
    'VariableName','numberDimensions',...
    'DefaultValue',2,...
    'Description', {
        'Number of desired (reduced) dimensions.'
        });

declareParameter(...
    'VariableName','numberPCADimensions',...
    'DefaultValue',50,...
    'Description', {
        'Number of dimensions for an initial PCA-based dimensionality'
        'reduction, prior to calling t-SNE.'
                   });

declareParameter(...
    'VariableName','seed',...
    'DefaultValue',0,...
    'Description', {
        'Seed (non-negative integer) used to initialize the random'
        'number generator used to construct the initial SOMby tSNE.'
                   });


declareParameter(...
    'VariableName','tsneOutput',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an output file, used for write access (output)'
        'Output of fast_bhtsne, displaying the errors.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','tbl',...
    'Description', {
        'Table that was processed.'
                   });
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,parameters__]=setParameters(nargout,varargin);
if stopNow
    return 
end

if ~isempty(structField)
    structure=tbl;
    tbl=structure.(structField);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tbl(1000:end,:)=[];

fprintf('tableTSNE: starting...');
t0=clock;

if ismember(tableLDColumn,tbl.Properties.VariableNames)
    tbl.Properties.VariableNames
    warning('tableTSNE: tableLDColumn (''%s'') already exists\n',tableLDColumn);
end

if ismember(tableCostColumn,tbl.Properties.VariableNames)
    tbl.Properties.VariableNames
    error('tableTSNE: tableCostColumn (''%s'') already exists\n',tableCostColumn);
end

xHD=tbl.(tableHDColumn);

if topicWeightThreshold > 0
    k=(xHD<topicWeightThreshold);     % creates index of all elements of 'xHD' that are below topicWeightThreshold
    xHD(k)=0;                         % changes all elements of 'xHD' in 'k' to zero
end

if numberPCADimensions>min(size(xHD))
    numberPCADimensions=min(size(xHD));
end

if size(xHD,2)>size(xHD,1)
    fprintf('doing myPCA (%dx%d->%d)...',size(xHD),numberPCADimensions);
    t1=clock();
    xHD=myPCA(xHD,numberPCADimensions);
    fprintf('done myPCA %dx%d (%.3f sec)\n',size(xHD),etime(clock,t1));
else
    fprintf('not doing myPCA...');
end
if BarnesHut
    % for Barnes-Hut t-SNE
    [xLD,landmarks,costs]=fast_bhtsne(xHD,numberDimensions,numberPCADimensions,...
                                      perplexity,theta,seed,...
                                      nIterations,pValueIteration,...
                                      momentumSwitchIteration,initialMomentum,finalMomentum,...
                                      eta,tsneOutput);
else
    % original fast_tsne
    fprintf('tableTSNE: ignoring seed\n');
    [xLD,landmarks,costs]=fast_tsne(xHD,numberDimensions,...
                                    numberPCADimensions,landmarks,perplexity);
end

if length(landmarks)~=size(tbl,1)
    size(tbl)
    size(landmarks)
    fprintf('tableTSNE: table height (%d) ~= length of landmarks (%d)\n',...
            size(tbl,1),length(landmarks));
end
if size(xLD,1)~=size(tbl,1)
    size(tbl)
    size(xLD)
    fprintf('tableTSNE: table height (%d) ~= length of xLD (%d)\n',size(tbl,1),size(xLD,1));
end

tbl.(tableLDColumn)=nan(size(tbl,1),size(xLD,2));
tbl{landmarks,tableLDColumn}=xLD;
tbl.(tableCostColumn)=nan(size(tbl,1),1);
tbl{landmarks,tableCostColumn}=costs;

% Code goes here

fprintf('done tableTSNE (%.3f sec)\n',etime(clock,t0));

if ~isempty(structField)
    structure.(structField)=tbl;
    tbl=structure;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varargout=setOutputs(nargout,parameters__);

