function varargout=tableDiversity(varargin);
% To get help, type tableDiversity('help')


% This script is a Mablab implementation of a set of algorithms for calculating
% diversity metrics as described in:
%
% Rafols, I. & Meyer, M. (2010) Diversity and
% network coherence as indicators of interdisciplinarity: Case studies in
% bionanoscience. Scientometrics 82(2): 263-287.
% http://link.springer.com/article/10.1007%2Fs11192-009-0041-y
% or http://arxiv.org/abs/0901.1380
%
% and Stirling, A. (2007) A general framework for analysing diversity
% in science, technology and society. Journal of The Royal Society
% Interface, 4:707-719. DOI: 10.1098/rsif.2007.0213.
% 
% Script is based on an R implementation of these functions 'diversity_measures_1.R'
% version 0.96, 21/6/2011 by Diego Chavarro
% http://www.interdisciplinaryscience.net/pub_docs/idr-local-files/script/
%
% Copyright (C) 2014  Stacy Rebich Hespanha & Joao Hespanha
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

declareParameter(...
    'Help', {
        'This script processes a table class object by executing:'
        '1)  selects a column that represents a matrix with positive'
        '    or zero numbers upon which to perform the following operations'
        '2)  changes to 0 all elements of the input matrix that are below'
        '    a specified threshold. This threshold can be set manually by'
        '    specifying a topic weight or by selecting a percentile value'
        '    for all topic weight values in the matrix. Both thresholding'
        '    methods cannot, however, be used at the same time.'
        '3)  computes distance between columns of the matrix, performs'
        '    nonlinear scaling (raise to an exponent) of the resulting'
        '    values and saves pairwise distance values as new matrix'
        '4)  computes 1 minus the sum of all elements in each row of the'
        '    new matrix'
        '5)  normalizes each row so that its elements add up to one'
        '6)  computes Shannon entropy for each row and saves these values'
        '    to a new column named "shannon"'
        '7)  uses Shannon entropy values to compute evenness measure for'
        '    each row and saves these values to a new column named "evenness"'
        '8)  computes Simpson diversity for each row and and saves'
        '    these values to a new column named "simpson"'
        '9)  computes Stirling diversity for each row and saves these'
        '    values to a new column named "stirling"'
        '10) computes presence-absence based disparity measure for each row'
        '    and saves these values to a new column named "disparity"'
        '11) computes topic-weighted disparity measure for each row and saves these'
        '    values to a new column named "disparityWeighted"'
        'The operations are performed in ths order.'
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
        'diversity calculator. This column should contain a (sparse) matrix'
        'of non-negative numbers in the form of one observation per row and'
        'one feature per column of the matrix.'
                   });
declareParameter(...
    'VariableName','percentileThreshold',...
    'DefaultValue',0,...
    'Description', {
        'This variable specifies the percentile value below which topic proportions'
        'should be treated as 0.'
                   });
declareParameter(...
    'VariableName','topicWeightThreshold',...
    'DefaultValue',0,...
    'Description', {
        'This variable specifies the topic weight below which topic proportions'
        'should be treated as 0.'
                   });
declareParameter(...
    'VariableName','columnDistance',...
    'DefaultValue','cosine',...
    'AdmissibleValues',{'sqEuclidean','cityblock','cosine','correlation','Hamming'},...
    'Description', {
        'This variable specifices the metric used for computing the distance matrix.'
                   });
declareParameter(...
    'VariableName','distanceExponent',...
    'DefaultValue','1',...
    'Description', {
        ['This variable specifices the exponent (must be positive) to be used for'...
         'reshaping (flattening) the distribution of values in the distance matrix'...
         'for dichotomous version of diversity calculation.']
                   });
declareParameter(...
    'VariableName','topicWeightRemovedTableColumn',...
    'DefaultValue','topicWeightRemoved',...
    'Description', {
        'Name of the table column where the the sum of the topic proportions that have been removed will be stored'
                   });
declareParameter(...
    'VariableName','varietyTableColumn',...
    'DefaultValue','variety',...
    'Description', {
        'Name of the table column where the variety values will be stored'
                   });
declareParameter(...
    'VariableName','shannonTableColumn',...
    'DefaultValue','shannon',...
    'Description', {
        'Name of the table column where the Shannon entropy values will be stored'
                   });

declareParameter(...
    'VariableName','evennessTableColumn',...
    'DefaultValue','evenness',...
    'Description', {
        'Name of the table column where the evenness values will be stored'
                   });

declareParameter(...
    'VariableName','simpsonTableColumn',...
    'DefaultValue','simpson',...
    'Description', {
        'Name of the table column where the Simpson diversity values will be stored'
                   });

declareParameter(...
    'VariableName','stirlingTableColumn',...
    'DefaultValue','stirling',...
    'Description', {
        'Name of the table column where the Stirling diversity values will be stored'
                   });

declareParameter(...
    'VariableName','disparityTableColumn',...
    'DefaultValue','disparity',...
    'Description', {
        'Name of the table column where the disparity values will be stored'
                   });

declareParameter(...
    'VariableName','disparityWeightedTableColumn',...
    'DefaultValue','disparityWeighted',...
    'Description', {
        'Name of the table column where the disparity values weighted by topic weights will be stored'
                   });

declareParameter(...
    'VariableName','closestDistancesTableColumn',...
    'DefaultValue','closestDistances',...
    'Description', {
        'Name of the table column where the average of the N closest distances to the topic with'
        'maximum topic weight will be stored. The number N of closest distances considered is defined'
        'by the parameter ''numberClosestDistances''.'
                   });

declareParameter(...
    'VariableName','numberClosestDistances',...
    'Description', {
        'Number of closest distances used in the creation of the column ''closestDistancesTableColumn''.'
        'Can be an array, in which case the table column ''closestDistancesTableColumn'' will contain'
        'multiple values.'
                   });


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','proportionsHistogram',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an .eps output file, used for write access (output)'
        'Name of the figure containing a histogram of the values of all elements'
        'of the input table.'
                   });
declareParameter(...
    'VariableName','proportionsPercentiles',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an .ascii output file, used for write access (output)'
        'Name of the vector containing percentile values for all elements of'
        'the input table.'
                   });
declareParameter(...
    'VariableName','closestDistancesHistogram',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an .eps output file, used for write access (output)'
        'Name of the file containing images of the histograms with the the average'
        'of the N closest distances to each topic, for different values of N.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','tbl',...
    'Description', {
        'Table that was processed and has diversity measures appended as'
        'additional columns.'
                   });

declareOutput(...
    'VariableName','distanceTable',...
    'Description', {
        'Table representing a matrix containing distances between'
        'columns of tableColumn transformed by specified "distanceExponent".'
        'This table has one row for each column of tableColumn,' 
        'containing a vector of distances to all the columns.'
                   });
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,params]=setParameters(nargout,varargin); % retrives parameters passed by other script
if stopNow
    return 
end

%verboseLevel=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('tableDiversity: ');        % prints 'tableDiversity: ' on the screen
t0=clock;                           % sets clock to zero  

%tbl(1000:end,:)=[];                % used to specify part of table to process

matrix=tbl.(tableColumn);           % creates variable 'matrix' that contains contents of table column of interest

nozerosmatrix=matrix(:);            % converts matrix to column vector 'nozerosmatrix'
k=(nozerosmatrix == 0);             % creates index of all elements 'nozerosmatrix' that are equal to zero
nozerosmatrix(k)=[];                % removes all indexed elements of nozerosmatrix
clearFigure('orientation','landscape','figureNumber',1,'figureName','Proportions Histogram');
%figure(1);clf
histogram(nozerosmatrix,500);            % creates histogram of all remaining elements in 'nozerosmatrix' using 500 bins
%print(proportionsHistogram,'-depsc2')  % prints histogram figure to .eps file
saveFigure('filename',proportionsHistogram,'saveFormat','eps');

percentiles = prctile(nozerosmatrix,[10,20,30,40,50,60,70,80,90]);  % calculates percentiles for nozerosmatrix
save(sprintf('%s.txt',proportionsPercentiles),'percentiles','-ascii');   % saves percentile values to plain text file

if percentileThreshold > 0 && topicWeightThreshold == 0
    proportionThreshold = prctile(nozerosmatrix,percentileThreshold);   % saves value associated with specified
                                                                        % percentileThreshold
elseif percentileThreshold == 0 && topicWeightThreshold >= 0
    proportionThreshold = topicWeightThreshold;     % sets proportionThreshold equal to specified topicWeightThreshold
else
    error('invalid to set both percentileThreshold and topicWeightThreshold to positive values\n');
end

k=(matrix<proportionThreshold);     % creates index of all elements of 'matrix' that are below proportionThreshold
                                    % specified 'proportionThreshold'
matrix(k)=0;                        % changes all elements of 'matrix' in 'k' to zero

ncolumns=size(matrix,2);            % creates variable 'ncolumns' that specifies the number of columns in 'matrix'

fprintf('calculating distance matrix... ');               % prints 'calculating distance matrix... ' on the screen
dm=pdist(matrix',columnDistance);
fprintf('turning distance matrix into square form 1... '); 
distanceMatrix=squareform(dm.^distanceExponent); % creates variable 'distanceMatrix' that contains the
                                                                           % square symmetric form of the distance matrix calculated
                                                                           % using the distance metric specified by
                                                                           % columnDistance and distanceExponent

fprintf('calculating topic weight removed... ');      % prints 'calculating topic weight removed... ' on the screen
tbl.(topicWeightRemovedTableColumn)=1-sum(matrix,2);  % calculates the sum of the elements removed (by 'junk' topic
                                                      % removal or percentile-based changes of values to zero) for
                                                      % each row of the matrix by summing the values of all elements
                                                      % in that row and subtracting from one
                                                      % saves value to 'topicWeightRemovedTableColumn'

fprintf('calculating normalized (proportions) matrix (%d rows, %d columns, %.3f sec)\n',...
        size(matrix,1),size(matrix,2),etime(clock,t0));   % prints 'calculating normalized (proportions) table ...'
                                                          % and includes information about number of rows and columns
                                                          % in the input matrix and processing time for this step
normalizedTableColumn=matrix./repmat(sum(matrix,2),1,ncolumns);  % creates variable 'normalizedTableColumn' that
                                                                 % contains the normalized (proportions) matrix
                                                                 % generated by summing each row of the matrix and
                                                                 % dividing each element in a row by the sum for that
                                                                 % row 

fprintf('calculating variety... ');                       % prints 'calculating variety... ' on the screen
tbl.(varietyTableColumn)=sum(matrix>0,2);                 % calculates variety score for each row of the matrix by
                                                          % counting the number of nonzero elements in that row, and
                                                          % saves variety scores to 'varietyTableColumn'

fprintf('calculating Shannon... ');                       % prints 'calculating Shannon... ' on the screen
shannon=normalizedTableColumn;                            % creates variable named 'shannon' that contains the
                                                          % normalized (proportions) matrix
k=find(shannon);                                          % finds indices and values of non-zero elements of
                                                          % the normalized (proportions) matrix
shannon(k)=shannon(k).*log(shannon(k));                   % multiplies each non-zero element of the normalized
                                                          % (proportions) matrix by the natural log of its own value
tbl.(shannonTableColumn)=-sum(shannon,2);                 % sums values resulting from calculation of
                                                          % value*log(value) for each row of matrix and saves these
                                                          % sums (Shannon entropy) to 'shannonTableColumn'

fprintf('calculating evenness... ');                      % prints 'calculating evenness... ' on the screen 
tbl.(evennessTableColumn)=tbl.(shannonTableColumn)./log(tbl.(varietyTableColumn));   % calculates evenness score by
                                                                                     % dividing (entry by entry) the
                                                                                     % Shannon entropy for each row
                                                                                     % by the natural log of the
                                                                                     % variety score for the same row
                                                                                     % and saves evenness score in
                                                                                     % 'evennessTableColumn' 

fprintf('calculating Simpson... ');                       % prints 'calculating Simpson... ' on the screen
tbl.(simpsonTableColumn)=1 - sum(normalizedTableColumn.^2,2);   % calculates Simpson diversity by summing the squared
                                                                % values of all elements in each row of the
                                                                % normalized (proportions) matrix and subtracting
                                                                % the resulting sum from 1; saves the resulting
                                                                % Simpson diversity score in 'simpsonTableColumn'

fprintf('calculating Stirling... ');                      % prints 'calculating Stirling... ' on the screen
tbl.(stirlingTableColumn)=sum((normalizedTableColumn * distanceMatrix) .* normalizedTableColumn,2);   % calculates
                                                          % Stirling diversity by multiplying columns of the
                                                          % normalized (proportions) matrix by rows of the distance
                                                          % matrix and then multiplying (entry by entry) the
                                                          % resulting product by the normalized
                                                          % (proportions) matrix; each row of the resulting product
                                                          % is then summed to generate the Stirling diversity value,
                                                          % which is saved in 'stirlingTableColumn'

fprintf('calculating disparity... ');                     % prints 'calculating disparity... ' on the screen
dichotomousTableColumn=normalizedTableColumn;             % creates variable named 'dichotomousTableColumn that
                                                          % contains the normalized (proportions) matrix
k=find(dichotomousTableColumn);                           % finds indices and values of non-zero elements of the
                                                          % normalized (proportions) matrix
dichotomousTableColumn(k)=1;                              % creates dichotomous version of the normalized
                                                          % (proportions) matrix by changing values of all non-zero
                                                          % elements of the normalized (proportions) matrix to 1
tbl.(disparityTableColumn)=sum(((dichotomousTableColumn * distanceMatrix) .* dichotomousTableColumn)./ ...
                               repmat(tbl.(varietyTableColumn).*(tbl.(varietyTableColumn) - 1),1,ncolumns),2);  %
                                                          % calculates disparity score by 1) multiplying (entry by
                                                          % entry) the variety scores by the variety scores minus 1,
                                                          % 2) multiplying columns of the dichotomous version of the
                                                          % normalized (proportions) matrix by rows of the distance
                                                          % matrix and then multiplying (entry by entry) the
                                                          % resulting product by the dichotomous version of the
                                                          % normalized (proportions) matrix, 3) dividing (entry by
                                                          % entry) the product of 2) by the product of 1), and
                                                          % summing across rows; disparity scores are saved in
                                                          % 'disparityTableColumn'


%% Compute average of the closest distances to all topics
closestDistances=sort(distanceMatrix,2);
closestDistances(:,1)=[]; % first distance is always zero
closestMean=cumsum(closestDistances,2)./repmat(1:size(closestDistances,2),size(closestDistances,1),1);
% plot histogram of closest distances
clearFigure('orientation','landscape','figureNumber',2,'figureName','histograms of average closest distances');
%figure(2);clf
for c=1:min(9,size(closestMean,2))
    subplot(3,3,c)
    histogram(closestMean(:,c),20)
    if ismember(c,numberClosestDistances)
        title(sprintf('>>average of %d closest distances<<',c))
    else
        title(sprintf('average of %d closest distances',c))
    end
end
%print(closestDistancesHistogram,'-depsc2')  % prints histogram figure to .eps file
saveFigure('filename',[closestDistancesHistogram,'1to9'],'saveFormat','eps');

% plot histogram of closest distances
clearFigure('orientation','landscape','figureNumber',3,'figureName','histograms of average closest distances');
%figure(3);clf
d=1;
for c=10:10:min(150,size(closestMean,2))
    subplot(3,5,d);d=d+1;
    histogram(closestMean(:,c),20)
    if ismember(c,numberClosestDistances)
        title(sprintf('>>average of %d closest distances<<',c))
    else
        title(sprintf('average of %d closest distances',c))
    end
end
%print(closestDistancesHistogram,'-depsc2')  % prints histogram figure to .eps file
saveFigure('filename',[closestDistancesHistogram,'10to150'],'saveFormat','eps');

selectedClosestMean=closestMean(:,numberClosestDistances);
% find maximum topic
[~,maxTopics]=max(matrix,[],2);
% add table to column
tbl.(closestDistancesTableColumn)=selectedClosestMean(maxTopics);

%% create remaining outputs
distanceTable=table(distanceMatrix);    % creates 'distanceTable' variable that contains tabular form of the distance
                                        % matrix 

fprintf('done tableDiversity (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t0));   % prints
                                        %'done tableDiversity on the screen and reports the size of the resulting
                                        %output table and information about processing time required

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);

