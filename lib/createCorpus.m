function varargout=createCorpus(varargin);
% To get help, type createCorpus('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'Creates a corpus from documents stored in a table. Specifically:'
        '. Each document is constructed from a table row by concatenating'
        '  several table columns. '
        '. Each document is parsed into word'
        '. Collocation breaks are determined based on the'
        '  ''collocationBreak'' marker'
        '. All words are converted to lower case'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','tbl',...
    'Description', {
        'Table containing the documents.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','textFields',...
    'DefaultValue','',...
    'Description', {
        'String array with the text columns to be concatenated.'
        });

declareParameter(...
    'VariableName','wordDelimiter',...
    'DefaultValue','[ .,:;"()\[\]?!\t\n\r/\\]+',...
    'Description', {
        'Regular expression that finds the break between words.'
        });

declareParameter(...
    'VariableName','columnSeparator',...
    'DefaultValue',' SeNtEnCeBrEaK ',...
    'Description', {
        'String to be inserted between the columns as they are merged. Typically,'
        '. starts and ends with a wordDelimiter character to make sure'
        '  the column separator becomes an independent word,'
        '. should include the collocationBreak to make sure collocations'
        '  do not cross columns.'
        });

declareParameter(...
    'VariableName','collocationBreak',...
    'DefaultValue','SeNtEnCeBrEaK',...
    'Description', {
        'Word that separates words that cannot form a collocation.'
        'Must be a whole word and therefore cannot have any wordDelimiter'
        'character, otherwise it will never match any word.'
                   });

declareParameter(...
    'VariableName','stopList',...
    'DefaultValue','ENstoplist.txt',...
    'Description', {
        'Filename for a .txt file, used for read access (input)'
        'Stop words list.'
        'Note that, by default, MATLAB searches the file on the'
        'path, if it cannot find it in the current directory.'
                   });

declareParameter(...
    'VariableName','minWordCount',...
    'DefaultValue',10,...
    'Description', {
        'Remove words with word counts smaller than this'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','dictionary',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for an output file, used for write access (output)'
        'This files contains a summary of the word usage.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','crps',...
    'Description', {
        'Corpus to be created (object of class corpus).'
        'See ''help corpus'' for the class internal structure.'
                   });
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Retrieve parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,params]=setParameters(nargout,varargin);
if stopNow
    return 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(stopList);

if (fid<0)
    error('createCorpus: unable to open stop list file ''%s''\n',stopList);
end

stopList=textscan(fid,'%s','CommentStyle','%');
stopList=stopList{1};
fclose(fid);

fprintf('Creating corpus... ');
t0=clock;

%% merge columns
docs=tableMergeColumns(tbl,textFields,columnSeparator);

fprintf('  Parsing documents to words... ');
t1=clock;
%% parse to words
maxLength=max(cellfun('length',docs));
word=[];
document=zeros(0,1,'uint32');
index=zeros(0,1,'int32');
collocationAllowed=zeros(0,1,'int8');
tokenBlockSize=10000000;
nTokens=0;
for doc=1:length(docs)
    if mod(doc,10000)==0
        fprintf('%d/%d (%.0fsec) ',doc,length(docs),etime(clock,t1));
    end
    
    if 0
        str=textscan(docs{doc},'%s','Delimiter',wordDelimiter,...
                     'MultipleDelimsAsOne',1,'BufSize',maxLength,'CollectOutput',1);
        str=str{1};
    else
        [str,~]=regexp(docs{doc},wordDelimiter,'split','match');
    end
    if nTokens+length(str)>length(word)
        fprintf('growing variables ');
        word{end+tokenBlockSize,1}=collocationBreak;
        document(end+tokenBlockSize,1)=1;
        index(end+tokenBlockSize,1)=0;
        collocationAllowed(end+1:end+tokenBlockSize,1)=1; % by default allow collocations
    end
    word(nTokens+1:nTokens+length(str),1)=str;
    document(nTokens+1:nTokens+length(str),1)=doc;
    index(nTokens+1:nTokens+length(str),1)=1:length(str);
    collocationAllowed(nTokens+1,1)=0; % remove collocation flags at start of document
    nTokens=nTokens+length(str);
end
%% truncate corpus
word(nTokens+1:end)=[];
document(nTokens+1:end)=[];
index(nTokens+1:end)=[];
collocationAllowed(nTokens+1:end)=[];

document=categorical(tbl.primaryKey(document));

fprintf('done (%d tokens, %.3f sec)\n',length(word),etime(clock,t1));

%% remove words that cannot be converted to categorical
fprintf('  Removing non-categorical words... ');
k=isnan(double(categorical(lower(word))));
word(k)=[];
document(k)=[];
index(k)=[];
collocationAllowed(k)=[];
fprintf('done %d tokens removed (%d tokens, %.3f sec)\n',sum(k),length(word),etime(clock,t1));

%% Determine collocation breaks
fprintf('  Setting collocationAllowed flags (and removing collocationBreak)... ');
t1=clock;
% remove collocation flags right after collocationBreak marker
k=find(strcmp(word(1:end-1),collocationBreak));
collocationAllowed(k+1)=0;
% remove all collocation markers
k=find(strcmp(word,collocationBreak));
word(k)=[];
document(k)=[];
index(k)=[];
collocationAllowed(k)=[];
fprintf('done (%d tokens, %.3f sec)\n',length(word),etime(clock,t1));

%% convert words to lowercase & only then to categorical
word=categorical(lower(word));

%% check for problems
checkErrors(tbl,word,document,textFields);

%% Create corpus class
crps=corpus(word,document,index,collocationAllowed);

%% Remove works in stoplist
fprintf(' Removing words in the stoplist...\n');
t1=clock;
removeWords(crps,stopList,false); % do not break collocations at words removed
fprintf(' done removing words in stoplist (%d tokens, %d words, %d documents, %.3f sec)\n',...
        size(crps.tokens,1),length(categories(crps.tokens.word)),...
        length(categories(crps.tokens.document)),etime(clock,t1));

%% Remove low-count words
fprintf(' Removing words with count<%d...\n',minWordCount);
t1=clock;
k=find(countcats(crps.tokens.word)<minWordCount);
words=categories(crps.tokens.word);
removeWords(crps,words(k),true); % break collocations at words removed
fprintf(' done removing words with count<%d (%d tokens, %d words, %d documents, %.3f sec)\n',...
        minWordCount,size(crps.tokens,1),length(categories(crps.tokens.word)),...
        length(categories(crps.tokens.document)),etime(clock,t1));


%% Output dictionary
fprintf(' Writing dictionary...\n',minWordCount);
t1=clock;
dict=table(categories(crps.tokens.word),countcats(crps.tokens.word),...
           'VariableNames',{'word','count'});
tableMultiWrite('tbl',dict,'outputFormats','tab','outputTable',dictionary);
fprintf(' done (%.3f sec)\n',etime(clock,t1))

fprintf('Done creating corpus (%d tokens, %d words, %d documents, %.3f sec)\n',...
        size(crps.tokens,1),length(categories(crps.tokens.word)),...
        length(categories(crps.tokens.document)),etime(clock,t0));
whos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkErrors(tbl,word,document,textFields)

    k=any(double(word)<1);
    if any(k)
        fprintf('ERROR: <1 token %d, doc %s\n',k,char(document(k)));
        disp(tbl(k,textFields))
    end
    k=any(isnan(double(word)));
    if any(k)
        fprintf('ERROR: nan token %d, doc %s\n',k,char(document(k)));
        disp(tbl(k,textFields))
    end

    