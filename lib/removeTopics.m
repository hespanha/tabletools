function varargout=removeTopics(varargin)
% To get help, type ldacol('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This script removes certain topics from the results of'
        'LDA (Latent Dirichlet Allocation).'
        'Removed topics are assigned NAN in crpsLDA.tokens.topic'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','crpsLDA',...
    'Description', {
        'Corpus with the LDA results (object of class corpus).'
        'Typically created by corpusLDA.'
        'See ''help corpus'' for the class'' internal structure.'
                   });

declareParameter(...
    'VariableName','docsLDA',...
    'Description', {
        'Table of documents with the LDA results (object of class table).'
        'Typically created by corpusLDA.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','topics2remove',...
    'DefaultValue','',...
    'Description', {
        'Array with the indices of the topics that need to be removed.'
                   });

declareParameter(...
    'VariableName','topicConcentrationThresholds',...
    'DefaultValue','',...
    'Description', {
        'Remove topics with a topicConcentration below this threshold.'
        'The topicConcentration is a value in the [0,1] interval that quantifies how' 
        'much the documents with high values for the topic are clustered around'
        'the averageDocWeights for that topic:'
        '   1 - the documents with high values for the topic are'
        '       tightly clustered around averageDocWeights'
        '   0 - the documents are not much clustered'
        'This value is computed as follows:'
        '  topicConcentration = Sum_{d=1}^{#docs}  weight(document d,topic t) * '
        '                   dist(document d, averageDocWeights for topic t)'
        '        / [ Sum_{d=1}^{# docs} weight(document d,topic t) ]'
        'where the sum is over all documents; the weights '
        '     weight(document d,topic t) '
        'correspond to the probability of the topic t  in the document d;'
        'and'
        '     dist(document d, averageDocWeights for topic t)'
        'is the cosine distance between document d and averageDocWeights for topic t.'
        'See ''help corpus''.'
                   });

declareParameter(...
    'VariableName','minSumWeights',...
    'DefaultValue',0.2,...
    'Description', {
        'Minimum value for the sum of topic weights.'
        'Documents with a sum below this are excluded.'
                   });

declareParameter(...
    'VariableName','renormalizeWeights',...
    'AdmissibleValues',{true,false},...
    'DefaultValue',true,...
    'Description', {
        'When true, the weights matrix is renormalized so that its'
        'columns add back to 1.'
                   });

declareParameter(...
    'VariableName','topicNameLength',...
    'DefaultValue',40,...
    'Description', {
        'Number of (top-weighted) words used to construct the name of a topic.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output Filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','topics',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for the topics files, used for write access (output)'
        'The topic descriptions will be written to {topics}.tab and'
        'a plot with the topic concentrations will be written to {topics}.eps'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','crpsLDA',...
    'Description', {
        'New corpus to be created (object of class corpus) with'
        'the LDA results included in the tokens and topics table.'
        'See ''help corpus'' for the class'' internal structure.'
                   });

declareOutput(...
    'VariableName','docsLDA',...
    'Description', {
        'Table of documents to be created. It will be a table class,'
        'with one document per row and the topics weights as a columns.'
                   });

declareOutput(...
    'VariableName','topicsLDA',...
    'Description', {
        'Table of topics to be created. It will be a table class,'
        'with one topic per row, with the following columns:'
        '   . wordWeights - array of the probability of each word for the topic '
        '                   (one probability for each word in the dictionary)'
        '   . topicName   - string array with the top-weighted words in the topic'
        '   . averageDocWeights - typically topic weights (probability distribution over'
        '                   topics) for a document that is representative of the topic, '
        '                   computed using'
        '                     averageDocWeights for topic t = '
        '                        = Sum_{d=1}^{# docs} weight(document d,topic t) *  '
        '                                             (vector of weights for document d)'
        '                           / [ Sum_{d=1}^{# docs} weight(document d,topic t) ]'
        '                   where the sum is over all documents and the weights '
        '                        weight(document d,topic t) '
        '                   corresponds to the probability of the topic t  in the document d.'
        '   . topicConcentration - value in the [0,1] interval that quantifies how much the documents'
        '                   with high values for the topic are cluster around the averageDocWeights'
        '                   for that topic:'
        '                      1 - the documents with high values for the topic are '
        '                          tightly clustered around averageDocWeights'
        '                      0 - the documents are not much clustered'
        '                   This value is computed as follows:'
        '                     topicConcentration = Sum_{d=1}^{#docs}  weight(document d,topic t) * '
        '                                      dist(document d, averageDocWeights for topic t)'
        '                           / [ Sum_{d=1}^{# docs} weight(document d,topic t) ]'
        '                   where the sum is over all documents; the weights '
        '                        weight(document d,topic t) '
        '                   correspond to the probability of the topic t  in the document d;'
        '                   and '
        '                        dist(document d, averageDocWeights for topic t)'
        '                   is the cosine distance between document d and averageDocWeights'
        '                   for topic t.'
        ' '
        'NOTE: This table is a "repeat" of a structure that is internal to ''crpsLDA.'''
        '      Eventually, this table should be removed from crpsLDA, but for now it is'
        '      still kept there for backwards compatibility.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to revise %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% crpsLDA.topics.topicSpread in code below should be revised to 
% crpsLDA.topics.topicConcentration once outputs of corpusLDA.m
% have been updated with this new terminology

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

topics2remove1=topics2remove;
topics2remove2=find(crpsLDA.topics.topicSpread<topicConcentrationThresholds);
topics2remove=union(topics2remove1,topics2remove2);

fprintf('removeTopics: removing %d topics (%d given+%d below threshold)... ',...
        length(topics2remove),length(topics2remove1),length(topics2remove2));
t0=clock;

%% remove topics
docsLDA.topicWeights(:,topics2remove)=[];
crpsLDA.topics(topics2remove,:)=[];
topicsLDA=crpsLDA.topics;
%% fix tokens
nTopics=max(crpsLDA.tokens.topic);
topics2keep=setdiff(1:nTopics,topics2remove);
change=nan(nTopics,1);
change(topics2keep)=1:length(topics2keep);
crpsLDA.tokens.topic=change(crpsLDA.tokens.topic);
fprintf('done (%.3f)... ',etime(clock,t0));

%% remove documents
sums=sum(docsLDA.topicWeights,2);
k=find(sums<minSumWeights);
if ~isempty(k)
    fprintf('%d low-weight documents... ',length(k));
    docsLDA(k,:)=[];
    fprintf('done (%.3f)... ',etime(clock,t0));
end

%% renormalize weights
if renormalizeWeights
    fprintf('renormalizing... ');
    docsLDA.topicWeights=docsLDA.topicWeights./repmat(sum(docsLDA.topicWeights,2),1,size(docsLDA.topicWeights,2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save topics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nTopics=size(docsLDA.topicWeights,2);
fprintf('\nSaving topics to ''%s.tab''\n',topics); 
t1=clock;
fid=fopen(sprintf('%s.tab',topics),'w');
fprintf('\nMost likely words in the topics:\n');
fprintf(fid,'Most likely words in the topics:\n');
% compressed format
fprintf(fid,'topic\tweight in docsLDA\ttopic concentration\twords\n');
for k=1:nTopics
    fprintf('Topic %3d (w=%.3f, s=%.3f): %s\n',...
            k,full(sum(docsLDA.topicWeights(:,k))),crpsLDA.topics.topicSpread(k),crpsLDA.topics.topicName{k}(1:min(end,100)));
    
    fprintf(fid,'%d\t%f\t%f\t%s\n',k,...
            full(sum(docsLDA.topicWeights(:,k))),crpsLDA.topics.topicSpread(k),crpsLDA.topics.topicName{k});
end
% long format
words=categories(crpsLDA.tokens.word);
wordWeights=crpsLDA.topics.wordWeights';
for k=1:nTopics
    [w,word]=sort(wordWeights(:,k),1,'descend');
    fprintf(fid,'Topic =\t%3d\tweight =\t%f\n',k,full(sum(docsLDA.topicWeights(:,k))));
    
    for l=1:min(topicNameLength,length(word))
        fprintf(fid,'\t\t\t\t%f\t%s\n',full(w(l)),words{word(l)});
    end
end
clear words
fclose(fid);
fprintf('   done saving topics ''%s'' (%.3fsec)\n',topics,etime(clock,t1));


fprintf('done (%.3f sec)\n',etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);





