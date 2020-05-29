function varargout=labelChart(varargin);
% To get help, type labelChart('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'Plots multiple text labels as a single text object.'
            });

declareParameter(...
    'VariableName','centers',...
    'Description', {
        '#labels x 2 matrix with the centers of each label.'
                   });
declareParameter(...
    'VariableName','labels',...
    'Description', {
        '#Pies x 1 vector of cell strings with the labels.'
                   });
declareParameter(...
    'VariableName','wrapDelimiter',...
    'DefaultValue','',...
    'Description', {
        'When nonempty, this variable provides a regular expression that is used'
        'to find ''new-lines'' in the label. The label is then wrapped over multiple lines.'
                   });
declareParameter(...
    'VariableName','linesPerLabel',...
    'DefaultValue',3,...
    'Description', {
        'Maximum number of lines per label.'
                   });
declareParameter(...
    'VariableName','fontSize',...
    'DefaultValue',[],...
    'Description', {
        'When a scalar, specifies the font size for the labels.'
        'When #labels x 1, specifies the font size of each individual label.'
        'When empty, the default FontSize is used.'
                   });

declareParameter(...
    'VariableName','textProperties',...
    'DefaultValue',{},...
    'Description', {
        'Cell array of property/value pairs to be passed to the text command.'
        'These could include:'
        '  ''FontSize'', <value in pts> (overridden by ''fontSize'' parameter)'
        '  ''EdgeColor'', [red,green,blue] or ''k'', ''b'', etc.'
        '  ''LineWidth'', <value in pts>'
        'To see all properties use set(line)'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','h',...
    'Description', {
        'handle of the text object'
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

if ~isempty(wrapDelimiter)
    labels=regexprep(labels,wrapDelimiter,'\n');
end

labels=cellfun(@(str)sscanf(str,'%[^\n]%[\n]',2*linesPerLabel-1),...
                     labels,'UniformOutput',0);
if isscalar(fontSize)
    fontSize=fontSize*ones(length(labels),1);
end

% remove empty labels
k=cellfun('isempty',labels);
centers(k,:)=[];
labels(k)=[];

if isempty(fontSize)
    h=text(centers(:,1),centers(:,2),ones(length(labels),1),labels,...
           'interpreter','none',...
           'verticalalignment','middle','horizontalalignment','center',textProperties{:});
else
    fontSize(k)=[];
    sizes=unique(fontSize);
    for i=1:length(sizes)
        k=find(sizes(i)==fontSize);
        h=text(centers(k,1),centers(k,2),ones(length(k),1),labels(k),...
               'interpreter','none',...
               'verticalalignment','middle','horizontalalignment','center',textProperties{:},...
               'FontSize',sizes(i));
    end
end
view(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);
