function varargout=tableSOMClustering(varargin);
% To get help, type tableSOMClustering('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This script clusters the rows of a table into "neurons" of a'
        'Self Organizing Map (SOM).'
        'This script is essentially a wrapper to interface matlab'
        'with SOM_PAK. It calls the appropriate SOM_PAK executables'
        'passing as inputs files created by this script.'
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
        'Name of the table''s column that contains the features used by the'
        'SOM classifier. This column should contain rows of numbers.'
        'In particular, tbl.(tableColumn) should result in a matrix,'
        'with one feature vector per row.'
                   });
declareParameter(...
    'VariableName','featureNames',...
    'DefaultValue',[],...
    'Description', {
        'String array with the names of the features used by the SOM classifier.'
        'This vector should have as many elements as the number of columns'
        'of the matrix tbl.(tableColumn).'
        'These names are used to name the SOM cells.'
        'In particular, if a SOM cell has the centroid with the largest value'
        'for dimension i-th dimension, then this cluster is named'
        'featureNames{i}. SOM cells that do not have the highest centroid'
        'for any dimension are named ''''.'
                   });
declareParameter(...
    'VariableName','somSize',...
    'DefaultValue',[26,30],...
    'Description', {
        'Two-vector with the x and y-dimensions of the SOM.'
                   });
declareParameter(...
    'VariableName','neighborhoodType',...
    'DefaultValue','gaussian',...
    'AdmissibleValues',{'gaussian','bubble'},...
    'Description', {
        'The neighborhood function type given to randinit,'
        'to be used for training (vsom).'
        });
declareParameter(...
    'VariableName','alphaType',...
    'DefaultValue','linear',...
    'Description', {
        'String specifying the type of adjustment for the alpha parameter.'
        'Can be either ''linear'' or ''inverse_t''.'
        'The linear function is defined as '
        '    alpha(t) = alpha(0)(1.0 - t/rlen)'
        'and the inverse-time type function is defined as'
        '    alpha(t) = alpha(0)C/(C + t)'
        'with C = rlen/100.0.'
        'According to the som_pak manual,'
        ' "It is advisable to use the inverse-time type function'
        '  with large maps and long training runs, to allow more'
        '  balanced finetuning of the reference vectors.'
                   });
declareParameter(...
    'VariableName','trainingLengths',...
    'DefaultValue',[100000,100000,500000],...
    'Description', {
        'Vector with number of iterations for each training set.'
        'som_pak''s vsom will be called once for each entry of this array.'
        'Each time, with a training matrix that has been re-ordered.'
                   });
declareParameter(...
    'VariableName','somRadius',...
    'DefaultValue',3*[10,5,1],...
    'Description', {
        'Vector with SOM''s radius parameter for each training set.'
        'The length of this vector should match the length of the'
        '''trainingLengths'' vector.'
                   });
declareParameter(...
    'VariableName','somAlpha',...
    'DefaultValue',[.05,.04,.03],...
    'Description', {
        'Vector with SOM''s alpha parameter for each training set.'
        'The length of this vector should match the length of the'
        '''trainingLengths'' vector.'
                   });
declareParameter(...
    'VariableName','snapshootIntervals',...
    'DefaultValue',[],...
    'Description', {
        'Vector with the number of iterations between snapshots.'
        'When a *negative* scalar is given, snapshooIntervals is choosen to be'
        '   -trainingLengths/snapshootIntervals',
        'which means that each training length is divided into the given'
        'number of equal intervals.'
        'When empty, no snapshots are produced.'
        'The length of this vector should match the length of the'
        '''trainingLengths'' vector.'
                   });
declareParameter(...
    'VariableName','computeQuantizationErrors',...
    'DefaultValue',false,...
    'AdmissibleValues',{true,false},...
    'Description', {
        'When true, quantization errors are computed at the start and end of'
        'each training length and at the snapshots.'
                   });
declareParameter(...
    'VariableName','fixedRow',...
    'DefaultValue',[],...
    'Description', {
        'Index of a ''tableColumn'' row (integer or name of a categorial array)'
        'to be assigned a fixed location (x, y) on the SOM during training.'
        'Ignored when empty.';
                   });
declareParameter(...
    'VariableName','fixedRowPosition',...
    'DefaultValue',[],...
    'Description', {
        '2-array containing the desired x,y coordinates for the row'
        'specified by ''fixedRow''.'
                   });
declareParameter(...
    'VariableName','fixedRowWeight',...
    'DefaultValue',[],...
    'Description', {
        'Scalar containing the training rate for the row specified by ''fixedRow''.'
        'The ''fixedRow'' vector is multiplied by this parameter so that'
        'the reference vectors are up dated as if this input vector were'
        'repeated this many times during training.'
                   });

declareParameter(...
    'VariableName','seed',...
    'DefaultValue',0,...
    'Description', {
        'Seed used to initialize the random number generator'
        'used to construct the initial SOM.'
                   });

switch (computer)
  case {'MACI64'}
    default=[fileparts(which('trainSOM')),'/../som_pak-3.1-OSX'];
  case {'GLNX86','GLNXA64'}
    default=[fileparts(which('trainSOM')),'/../som_pak-3.1-linux'];
  otherwise
    default='';
end

declareParameter(...
    'VariableName','path2sompak',...
    'DefaultValue',default,...
    'Description', {
        'Path to the folder containing the SOM_PAK executables.'
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
%% Output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','trainingMatrix',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for the output files, used for write access (output)'
        'Files containing a version of tbl.(tableColumn) in a'
        'tab-separated format readble by SOM_PAK:'
        '1) The first row contains the number of columns (i.e., number of topics)'
        '2) Each subsequent row corresponds to one document, and contains'
        '   the weights of the different topics on different columns.'
        '   The last column of each row contains a unique document ID'
        '   number (essentially reflecting the document order in'
        '   the input matrix).'
        'These files are used by SOM_PAK''s ''vsom'' to train the SOM.'
        'Each training file will be named {trainingMatrix+1.txt}'
        'The different versions of the matrix in each file produced'
        'differ only by the order in which the documents appear.'
        'The construction of these multiple training files permits'
        'running ''vsom'' multiple times with different document orders.'
        'One training file is needed for each entry in the vector'
        'of ''trainingLengths.'''
                   });

declareParameter(...
    'VariableName','calibrateMatrix',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for the output files, used for write access (output)'
        'Matrix tbl.(tableColumn)in a tab-separated format readable'
        'by som_pack:'
        '1) The first row contains the number of columns (i.e., number of topics)'
        '2) Each subsequent row corresponds to one document, and contains'
        '   the weights of the different topics on different columns.'
        '   The last column of each row contains a unique document ID'
        '   number (essentially reflecting the document order in'
        '   the input matrix).'
        'This file is used by SOM_PAK''s ''visual'' to calibrate the SOM.'
        'This file is typically created by createSOM'
                   });

declareParameter(...
    'VariableName','quantizationError',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an output file, used for write access (output)'
        'Output of SOM_PAK''s ''qerror'', displaying the quantization error.'
        'When computeQuantizationErrors is false, this variable is ignored.'
                   });

declareParameter(...
    'VariableName','snap',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an output folder, used for write access (output)'
        'Snapshots of the code file produced by SOM_PAK''s ''vsom'''
        'during training are saved inside this folder with the name'
        'snap+X+YYYY.cod, where X stands for index of the training file'
        'and YYYY for the iteration number.'
        'These files contain the vectors associated with each neuron.'
        'of the SOM. Their formats are described in the SOM_PAK''s'
        'documentation.'
                   });

declareParameter(...
    'VariableName','neurons',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an output file, used for write access (output)'
        'Code file produced by SOM_PAK''s ''vsom'' after training.'
        'This file contains the vectors associated with each neuron'
        'of the SOM. Its format is described in the SOM_PAK''s'
        'documentation. This file can be converted to a shape file'
        'using the ArcGIS toolbox SOMAnalyst.'
                   });

declareParameter(...
    'VariableName','umat',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an output file, used for write access (output)'
        'Images for an umat visualization of the neurons.'
                   });

declareParameter(...
    'VariableName','docs',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an output file, used for write access (output)'
        'Code file produced by SOM_PAK''s ''visual'' after calibrating'
        'the documents. For each document (one per row), this file contains'
        'the x-y coordinates of the best-matching neuron and the associated'
        'quantization error. This file can be converted to a shape file'
        'using the ArcGID toolbox SOMAnalyst.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,params]=setParameters(nargout,varargin);
if stopNow
    return 
end

if length(snapshootIntervals)==1 && snapshootIntervals<0
    snapshootIntervals=-trainingLengths/snapshootIntervals;
end

if ~computeQuantizationErrors
    quantizationError='';
end

%verboseLevel=4;
execute=1; % execute commands

initializeAndTrain=1;  % initialize and train
qerrorDisplay=1;       % display quatization error
calibrate=1;           % calibrate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('tableSOMClustering: starting\n');
t0=clock;

dataMatrix=tbl.(tableColumn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Remove documents with NaN values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kNaN=find(all(~isnan(dataMatrix),2));

if size(dataMatrix,1)>length(kNaN)
    fprintf('  %d out of %d documents have NaN counts - will be ignored\n',...
            size(dataMatrix,1)-length(kNaN),size(dataMatrix,1));
    dataMatrix=dataMatrix(kNaN,:);
end

fixedRow=double(fixedRow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output matrices for SOM training with random row orders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if initializeAndTrain
    for i=1:length(trainingLengths)
        filename=sprintf('%s+%d',trainingMatrix,i);
        %t = getCurrentTask();
        %fprintf('worker %d: Saving SOM training data ''%s''\n',t.ID,filename);
        fprintf('  saving SOM training data ''%s''... ',filename);
        t1=clock;
        [dummy,k]=sort(rand(size(dataMatrix,1),1));
        labels=cellstr(num2str(k));
        if ~isempty(fixedRow)
            l=find(fixedRow==k);
            if ~isempty(fixedRowWeight)
                fprintf('  using weight %d for document %s with index %d, ',...
                        fixedRowWeight,fixedRow,l);
                labels{l}=sprintf('%s weight=%d',labels{l},fixedRowWeight);
            end
            fprintf('  using fixed location (%d,%d) for document %s with index %d, ',...
                    fixedRowPosition,fixedRow,l); 
            labels{l}=sprintf('%s fixed=(%d,%d)',labels{l},fixedRowPosition);
        end
        savematrix('outputPrefix',filename,...
                   'dataMatrix',dataMatrix(k,:),...
                   'labels',labels);
        fprintf('done (%.3f sec)\n',etime(clock,t1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output matrices for SOM calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('  saving SOM visualize docs data ''%s''... ',calibrateMatrix);
labels=cellstr(num2str((1:size(dataMatrix,1))'));
savematrix('outputPrefix',calibrateMatrix,...
           'dataMatrix',dataMatrix,...
           'labels',labels);
fprintf('done (%.3f sec)\n',etime(clock,t1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run sompak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(somSize)~=2
    somSize
    error('trainSOM: ''somSize'' input should be a 2-vector\n',length(somSize));
end

if length(trainingLengths)~=length(somRadius) ||...
        length(trainingLengths)~=length(somAlpha)
    trainingLengths
    somRadius
    somAlpha
    error('trainSOM: the inputs ''trainingLengths'', ''somRadius'', ''somAlpha'' should all have the same length\n');
end

t0=clock;

if initializeAndTrain
    if ~isempty(quantizationError)
        % remove quantizationError file
        cmd=sprintf('rm -f %s.txt',quantizationError);
        fprintf('%s\n',cmd);
        if execute
            system(cmd);            
        end
    end
    if ~ isempty(snapshootIntervals)
        % clear snap folder and create if it does not exist
        cmd=sprintf('rm -f %s/*.cod',snap);
        fprintf('%s\n',cmd);
        if execute
            system(cmd);            
        end
        cmd=sprintf('mkdir %s',snap);
        fprintf('%s\n',cmd);
        if execute
            system(cmd);            
        end
    end

    %% Initialize map
    fprintf('\nInitializing the SOM\n')
    filename=sprintf('%s+%d',trainingMatrix,1);
    cmd=sprintf(['randinit -din "%s.txt" -cout "%s.cod" ',...
                 '-xdim %.0f -ydim %.0f -topol hexa -neigh %s -rand %f'],...
                filename,neurons,somSize(1),somSize(2),neighborhoodType,seed);
    fprintf('%s\n',cmd);
    if execute
        t1=clock;
        system(sprintf('%s/%s',path2sompak,cmd));            
        fprintf('   done Initializing the SOM (%.2f secs)\n',etime(clock,t1));
    end

    fprintf('\nTraining the SOM\n')
    for i=1:length(trainingLengths)

        if ~isempty(quantizationError) && i==1
            %% Compute current quantization error
            fprintf('\nSOM initial quantization error %d/%d\n',i,length(trainingLengths));
            cmd=sprintf('qerror -din "%s.txt" -cin "%s.cod" | tee -a "%s.txt"',...
                        calibrateMatrix,neurons,quantizationError);
            fprintf('%s\n',cmd);
            if execute
                t1=clock;
                system(sprintf('%s/%s',path2sompak,cmd));            
                fprintf('   done SOM initial quantization error (%.2f secs)\n',etime(clock,t1));
            end
            if qerrorDisplay
                displayQerror(quantizationError,trainingLengths,somAlpha,somRadius)
                print('-depsc',sprintf('%s.eps',quantizationError))
            end
        end
        
        %% Train
        if trainingLengths(i)>0
            fprintf('\nSOM training %d/%d\n',i,length(trainingLengths));
            filename=sprintf('%s+%d',trainingMatrix,i);
            if isempty(snapshootIntervals)
                snapOption='';
            else
                snapOption=sprintf('-snapinterval %d -snapfile "%s/snap+%d+%%09d.cod" ',...
                             snapshootIntervals(i),snap,i);
            end
            cmd=sprintf(['vsom -din "%s.txt" -cin "%s.cod" '...
                         '-cout "%s.cod" %s ' ...
                         '-fixed %d -weight %d ' ...
                         '-rlen %.0f -alpha %f -radius %f -rand %f'],...
                        filename,neurons,neurons,snapOption,...
                        ~isempty(fixedRow),~isempty(fixedRowWeight),...
                        trainingLengths(i),somAlpha(i),somRadius(i),seed);
            fprintf('%s\n',cmd);
            if execute
                t1=clock;
                system(sprintf('%s/%s',path2sompak,cmd));            
                fprintf('   done training length %d (%.2f secs)\n',trainingLengths(i),etime(clock,t1));
            end
        end

        if ~isempty(quantizationError)
            %% Compute snapshot quantization error
            if ~isempty(snapshootIntervals)
                fprintf('\nSOM snapshots quantization errors %d/%d\n',i,length(trainingLengths))
                files=dir(sprintf('%s/snap+%d+*.cod',snap,i));
                for j=1:length(files)
                    iter=regexp(files(j).name,sprintf('snap\\+%d\\+(\\d*)\\.cod',i),'tokens');
                    if ~isempty(iter)
                        iter=str2num(iter{1}{1});
                        cmd=sprintf('qerror -din "%s.txt" -cin "%s/snap+%d+%09d.cod" >>"%s.txt"',...
                                    calibrateMatrix,snap,i,iter,quantizationError);
                        %fprintf('%s\n',cmd);
                        if execute
                            t1=clock;
                            system(sprintf('%s/%s',path2sompak,cmd));            
                            fprintf('   done SOM snapshots quantization error %d/%d (%.2f secs)\n',j,length(files),etime(clock,t1));
                        end
                        if qerrorDisplay
                            displayQerror(quantizationError,trainingLengths,somAlpha,somRadius)
                        end
                    end
                end
            end
        end
    
        if ~isempty(quantizationError)
            %% Compute final quantization error
            fprintf('\nSOM final quantization error %d/%d\n',i,length(trainingLengths));
            cmd=sprintf('qerror -din "%s.txt" -cin "%s.cod" | tee -a "%s.txt"',...
                        calibrateMatrix,neurons,quantizationError);
            fprintf('%s\n',cmd);
            if execute
                t1=clock;
                system(sprintf('%s/%s',path2sompak,cmd));            
                fprintf('   done SOM final quantization error (%.2f secs)\n',etime(clock,t1));
            end
            if qerrorDisplay
                displayQerror(quantizationError,trainingLengths,somAlpha,somRadius)
                print('-depsc',sprintf('%s.eps',quantizationError))
            end
        end
        
    end

end

if calibrate
    %% Calibrate
    fprintf('\nCalibrating the SOM\n')
    cmd=sprintf('visual -din "%s.txt" -cin "%s.cod" -dout "%s.bmu" -noskip 1',...
                calibrateMatrix,neurons,docs);
    fprintf('%s\n',cmd);
    if execute
        t1=clock;
        system(sprintf('%s/%s',path2sompak,cmd));            
        fprintf('   done Calibrating the SOM (%.2f secs)\n',etime(clock,t1));
    end
end

if ~isempty(umat)
    %% u-matrix visualization
    fprintf('u-matrix visualization of the SOM\n')
    cmd=sprintf('umat -cin "%s.cod" -paper A4 -fontsize 1.7 -o "%s.eps"',neurons,umat);
    fprintf('%s\n',cmd);
    if execute
        t1=clock;
        system(sprintf('%s/%s',path2sompak,cmd));            
        fprintf('   done u-matrix visualization of the SOM (%.2f secs)\n',etime(clock,t1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read SOM data about neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neuronsSOMfilename=sprintf('%s.cod',neurons);
fprintf('Reading SOM neurons (%s)... ',neuronsSOMfilename);
fid=fopen(neuronsSOMfilename);
nTopics=fscanf(fid,'%d');
neurons.shape=fscanf(fid,'%s',1);
neurons.nx=fscanf(fid,'%d',1);
neurons.ny=fscanf(fid,'%d',1);
neurons.influence=fscanf(fid,'%s',1);
neurons.topics=fscanf(fid,'%g',[nTopics,Inf])';
fclose(fid);
fprintf('done!\n');

neurons.nNeurons=neurons.nx*neurons.ny;

switch neurons.shape
  case 'hexa',
    neurons.xCoords=(0:neurons.nx-1)'*ones(1,neurons.ny);
    neurons.xCoords(:,2:2:end)=neurons.xCoords(:,2:2:end)+.5;
    neurons.yCoords=ones(neurons.nx,1)*(0:neurons.ny-1)*sqrt(3)/2;
  case 'rect',
    neurons.xCoords=(0:neurons.nx-1)'*ones(1,neurons.ny);
    neurons.yCoords=ones(neurons.nx,1)*(0:neurons.ny-1);
  otherwise
    error('tableSOMClustering: unknown neuron shape %s\n',neurons.shape);
end
neurons.xCoords=reshape(neurons.xCoords,neurons.nNeurons,1);
neurons.yCoords=reshape(neurons.yCoords,neurons.nNeurons,1);

if size(neurons.topics,1)~=neurons.nNeurons
    error('tableSOMClustering: %d neurons found, but %d expected\n',...
          size(neurons.topics,1),neurons.nNeurons);
end

clusterPositions=[neurons.xCoords,neurons.yCoords];
clusterCentroids=neurons.topics;
clusterIDs=(1:neurons.nNeurons)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read SOM data about documents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

docsSOMfilename=sprintf('%s.bmu',docs);
fprintf('Reading SOM documents (%s)... ',docsSOMfilename);
fid=fopen(docsSOMfilename);
x=fscanf(fid,'%d',1);
if x~=3
    error('tableSOMClustering: unexpected number of columns in ''%s'', 3 expected, %d found\n',...
          docsSOMfilename,x);
end
x=fscanf(fid,'%s',1);
if neurons.shape~=x
    error('tableSOMClustering: mismatch between neurons in ''%s'' and ''%s''\n',...
          neuronsSOMfilename,docsSOMfilename);
end
if neurons.nx~=fscanf(fid,'%d',1) || neurons.ny~=fscanf(fid,'%d',1) || ~strcmp(neurons.influence,fscanf(fid,'%s',1))
    error('tableSOMClustering: mismatch between neurons in ''%s'' and ''%s''\n',...
          neuronsSOMfilename,docsSOMfilename);
end
xyq=fscanf(fid,'%g',[3,Inf])';
fclose(fid);
fprintf('done!\n');

nDocuments=size(xyq,1);
        
clusteringIDs=xyq(:,1)+1+neurons.nx*xyq(:,2);
clusteringQerror=xyq(:,3);

tblcluster=tableclusters(tbl.primaryKey,clusteringIDs,clusteringQerror,...
                         clusterIDs,clusterCentroids,clusterPositions);
tblcluster=addClusterNamesFromCentroids(tblcluster,featureNames,-inf);

fprintf('done tableSOMClustering (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displayQerror(quantizationError,trainingLengths,somAlpha,somRadius)

%% Read and display quantization errors
filename=sprintf('%s.txt',quantizationError);
fid=fopen(filename);
qerror=[];
snapshots={};
while 1
    x=fgets(fid);
    if isequal(x,-1)
        break;
        %type(filename);
        %error('Unable to get quantization error from file ''%s''\n',filename);
    end
    x=textscan(x,'Quantization error of %[^ ] with map %[^ ] is %f per sample (%d samples)\n');
    if ~isempty(x{3})
        qerror(end+1,1)=x{3};
        s=regexp(x{2},'/snap\+(\d*)\+(\d*)\.cod','tokens');
        snapshots{end+1,1}=s{1};
    end
end
fclose(fid);
%qerror=x{3};
%snapshots=regexp(x{2},'/snap\+(\d*)\+(\d*)\.cod','tokens');

% find iterations & training file for finals
finals=cellfun(@(x) isempty(x),snapshots);
nfinals=[0;trainingLengths'];
iter(finals)=nfinals(1:sum(finals));
testFile(finals)=0:sum(finals)-1;

% find iterations & training file for snapshots
snaps=~finals;
if any(snaps)
    cellfun(@(x) str2num(x{1}{2}),snapshots(snaps));
    iter(snaps)=cellfun(@(x) str2num(x{1}{2}),snapshots(snaps));
    testFile(snaps)=cellfun(@(x) str2num(x{1}{1}),snapshots(snaps));
end

finals=find(finals);

% plot
figure(1);clf
h=[];
leg={};
iter0=0;
for i=1:length(trainingLengths)
    k=find(testFile==i);
    if ~isempty(k)
        iter(k)=iter(k)+iter0;
        k=min(k)-1:max(k); % start from last one
        h(end+1)=plot(iter(k),qerror(k),'-*b');
        grid on
        leg{end+1}=sprintf('train %d, \\alpha=%g, radius=%g',i,somAlpha(i),somRadius(i));
        hold on
        iter0=iter(k(end));
    end
end
if 0
    fprintf('iterations:\n');
    disp(iter)
    fprintf('final iterations:\n');
    disp(iter(finals))
end

h(end+1)=plot(iter(finals),qerror(finals),'*r');
leg{end+1}='finals';

% focus axis on smaller errors
axis([-1,max(iter)+1,0*min(qerror),1.1*min(2*min(qerror),max(qerror))]);
grid on
ylabel('quantization error');
xlabel('# of iterations');
title(sprintf('Final quantization error = %g (initial = %g, min = %g)',...
              qerror(finals(end)),qerror(finals(1)),min(qerror)));
legend(h,leg)

drawnow

