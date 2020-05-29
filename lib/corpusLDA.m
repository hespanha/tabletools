function varargout=corpusLDA(varargin)
% To get help, type ldacol('help')
%
% Copyright (C) 2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This script essentially performs LDA (Latent Dirichlet Allocation)'
        'with Collocations on a corpus class.'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','crps',...
    'Description', {
        'Corpus to be updated (object of class corpus).'
        'See ''help corpus'' for the class'' internal structure.'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','nTopics',...
    'DefaultValue',150,...
    'Description', {
        'Number of topics'
                   });

declareParameter(...
    'VariableName','collocations',...
    'AdmissibleValues',{false,true},...
    'DefaultValue',true,...
    'Description', {
        'When ''true'' collocations are allowed.'
                   });

declareParameter(...
    'VariableName','alphaTopics',...
    'DefaultValue',0,...
    'Description', {
        'Parameter of the Dirichelet prior on the per-document topic distribution.'
        'Small values favor sparse distributions with most probability concentrated'
        'on just a few topics.'
        'When 0, the value suggested in the topictoolbox of 50/(# topics) is used.'
                   });

declareParameter(...
    'VariableName','alphaWords',...
    'DefaultValue',0,...
    'Description', {
        'Parameter of the Dirichelet prior on the per-topic word distribution.'
        'Small values favor sparse distributions with most probability concentrated'
        'on just a few words.'
        'When 0, the value suggested in the topictoolbox of 200/(# words in vocabulary) is used.'
                   });

declareParameter(...
    'VariableName','Delta',...
    'DefaultValue',.1,...
    'Description', {
        'Parameter of the Dirichelet prior on the per-word next-word-in-collocation'
        'distribution. Small values favor sparse next-word distributions with most'
        'probability concentrated on just a few words.'
        'The value used in the topictoolbox examples is .1'
                   });

declareParameter(...
    'VariableName','Gamma0',...
    'DefaultValue',.1,...
    'Description', {
        'First parameter of the Beta prior on the per-word probability'
        'of next word being a collocation.'
        'For Gamma0=Gamma1<1, the distribution is U-shaped around and symmetric'
        'around 1/2, leading to an equal prior to form a collocation or not'
        'but favoring each word to consistently be in a collocation or not.'
        'The smaller Gamma0=Gamma1, the more one favors consistency.'
        'For Gamma0<Gamma1, the distribution is biased towards more colocations.'
        'The values used in the topictoolbox examples are Gamma0=Gamma1=.1'
                   });

declareParameter(...
    'VariableName','Gamma1',...
    'DefaultValue',.1,...
    'Description', {
        'Second parameter of the Beta prior on the per-word probability'
        'of next word being a collocation.'
        'For Gamma0=Gamma1<1, the distribution is U-shaped around and symmetric'
        'around 1/2, leading to an equal prior to form a collocation or not'
        'but favoring each word to consistently be in a collocation or not.'
        'The smaller Gamma0=Gamma1, the more one favors consistency.'
        'For Gamma0<Gamma1, the distribution is biased towards more colocations.'
        'The values used in the topictoolbox examples are Gamma0=Gamma1=.1'
                   });

declareParameter(...
    'VariableName','topicNameLength',...
    'DefaultValue',40,...
    'Description', {
        'Number of (top-weighted) words used to construct the name of a topic.'
                   });

declareParameter(...
    'VariableName','nIterations',...
    'DefaultValue',[200,100,100,100],...
    'Description', {
        'Number of iterations.'
        'When a vector is given, multiple intermediate outputs are generated:'
        'the first after nIterations(1) iterations, the second after '
        'nIterations(1)+nIterations(2) iterations, etc.'
        'This is useful to check the convergence progress.'
                   });

declareParameter(...
    'VariableName','nReplicates',...
    'DefaultValue',1,...
    'Description', {
        'Number of times the words of each document is replicated for Gibbs sampling.'
        'A value ofr ''nReplicates'' larger than 1 may be useful when some documents'
        'have a small number of words for which a single (Gibbs) sampling would not'
        'suffice to construct a good estimate of the topic probabilities for the document.'
        ''
        'The total size of the corpus (in number of tokens) is multiplied by ''nReplicates'''
        '(making each iteration slower), but the total number of  iterations needed should'
        'roughly decrease by the same amount.'
                   });

declareParameter(...
    'VariableName','nThreads',...
    'DefaultValue',1,...
    'Description', {
        'Number of threads to be used (should be no larger than number of CPUs/Cores).'
        'ATTENTION: when nThreads>1, the results are not repeatable (even with the same seed)'
        '           because of ''races'' among threads. This is not an error.'
                   });

declareParameter(...
    'VariableName','checkSums',...
    'DefaultValue',false,...
    'AdmissibleValues',{false,true},...
    'Description', {
        'When 1 performs a check that can be used to detect ''races'' in multi-processor'
        'systems. Makes the computation slower and at this point the check should not'
        'be necessary since lock''s are being used to prevent ''races''.'
                   });

declareParameter(...
    'VariableName','seed',...
    'DefaultValue',0,...
    'Description', {
        'Seed used to initialize the random number generator for Gibbs sampling.'
        'See comment on nThreads regarding repeatibility of results.'
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
        'a plot with the topic spreads will be written to {topics}.eps'
                   });

declareParameter(...
    'VariableName','convergence',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for the convergence plots, used for write access (output)'
        'A plot showing the LDA convergence will be written to this file.'
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
        '   . topicSpread - value in the [0,1] interval that quantifies how much the documents'
        '                   with high values for the topic are cluster around the averageDocWeights'
        '                   for that topic:'
        '                      1 - the documents with high values for the topic are '
        '                          tightly clustered around averageDocWeights'
        '                      0 - the documents are not much clustered'
        '                   This value is computed as follows:'
        '                     topicSpread = Sum_{d=1}^{#docs}  weight(document d,topic t) * '
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
%% Retrieve parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stopNow,params]=setParameters(nargout,varargin);
if stopNow
    return 
end

nWords=length(categories(crps.tokens.word));
nTokens=size(crps.tokens,1);
nDocs=length(categories(crps.tokens.document));

if alphaTopics==0
    alphaTopics=50/nTopics;
end
if alphaWords==0
    alphaWords=200/nWords;
end

streamSelect=RandStream.create('mt19937ar','seed',seed);
RandStream.setGlobalStream(streamSelect);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function body
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fake corpus for testing re-order
if 0
    nDocs=5;
    crps.tokens=table([1;1;1;3;3;2;2;2;4;4;4;5;5;5],...
                      [11;12;13;31;32;21;22;23;41;42;43;51;52;53],...
                      [1;2;3;1;2;1;2;3;1;2;3;1;2;3],...
                      'VariableNames',{'document','word','index'});
    crps.docNames={'1';'2';'3';'4';'5'}
    crps.tokens
end

if ~collocations
    crps.tokens.collocationAllowed=zeros(nTokens,1,'uint8'); 
end

%% replication of documents
if nReplicates>1
    t0=clock;
    fprintf('Replicating tokens %d times, before: %d tokens,... ',nReplicates,nTokens);
    crps.tokens=repmat(crps.tokens,nReplicates,1);
    % prevent duplicate indices
    k=crps.tokens.index;
    k=max(-k(2:end)+k(1:end-1),0);
    k(k>0)=k(k>0)+1;
    k=[0;cumsum(k)];
    crps.tokens.index=crps.tokens.index+k;
    nTokens=size(crps.tokens,1);
    fprintf('after: %d tokens, done (%.3fsec)\n',nTokens,etime(clock,t0));
end

%% random re-order of documents
k=uint32(randperm(nDocs)');
crps.tokens.reorder=k(crps.tokens.document);
crps.tokens=sortrows(crps.tokens,{'reorder','index'},'ascend');
crps.tokens.reorder=[];

%% Initialize latent variables

% ATTENTION: z is changed by reference so one cannot simply assign it equal to another variable
z=uint32(ceil(nTopics*rand(nTokens,1)));

% ATTENTION: x is changed by reference so one cannot simply assign it equal to another variable
if collocations
    % initialize with every collocation (seems to converge a bit faster)
    x=1*uint32(crps.tokens.collocationAllowed); % the 1* is needed!!!!
else
    % initialize with no collocations
    x=zeros(nTokens,1,'uint32');  
end

%% Create word-pairs counts:
k=find(crps.tokens.collocationAllowed);
w1=double(crps.tokens.word(k-1));
w2=double(crps.tokens.word(k));
wordPairs=sparse(w2,w1,ones(length(k),1),nWords,nWords);
ind=sub2ind(size(wordPairs),double(crps.tokens.word(2:end)),double(crps.tokens.word(1:end-1)));
ind=uint32(full(wordPairs(ind)));
ind(end+1)=0;
crps.tokens.wordPairs=ind;

if any(sum(crps.tokens.collocationAllowed)~=sum(sum(wordPairs)))
    error('corpusLDA: mistmatch between sum(crps.tokens.collocationAllowed)=%d and sum(wordPairs)=%d\n',sum(crps.tokens.collocationAllowed),sum(sum(wordPairs)))
end
clear ind wordPairs w2 w1 k

%% Create sum variable
sumTW=zeros(nTopics,nWords); % ATTENTION: changed by ref, do not assign equal to other variable
sumTD=zeros(nTopics,nDocs); % ATTENTION: changed by ref, do not assign equal to other variable
sumT=zeros(nTopics,1);  % ATTENTION: changed by ref, do not assign equal to other variable
sumW0=zeros(nWords,1); % ATTENTION: changed by ref, do not assign equal to other variable
sumW1=zeros(nWords,1); % ATTENTION: changed by ref, do not assign equal to other variable


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% loop over nIterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verboseLevel
    whos
end

t0=clock;

ldacol.logL=zeros(sum(nIterations),1);
ldacol.percentCollocations=zeros(sum(nIterations),1);

iter=0;
skipInitialSums=0;
maxCollocations=sum(crps.tokens.collocationAllowed);

for i=1:length(nIterations)
    for j=1:nIterations(i)
        iter=iter+1;
        fprintf('Iteration %5d/%d (%5.0fsec)... ',iter,sum(nIterations),etime(clock,t0));

        %%%%%%%%%%%%%%%%%%%%
        %% Gibbs Sampling
        %%%%%%%%%%%%%%%%%%%%
        t1=clock;
        u=rand(nTokens,1);
        v=rand(nTokens,1);
        [ldacol.logL(iter)]=myGibbsSamplingLDAcol_updatecounts(...
            sumTW,sumTD,sumT,sumW0,sumW1,...
            uint64(crps.tokens.word),uint64(crps.tokens.document),...
            uint32(crps.tokens.collocationAllowed),crps.tokens.wordPairs,...
            u,v,z,x,...
            alphaTopics,alphaWords,Delta,Gamma0,Gamma1,...
            skipInitialSums,nThreads);
        skipInitialSums=1;
        fprintf('mex (%.3fsec)... ',etime(clock,t1));

        %%%%%%%%%%%%%%%%%%%%
        %% Check Sums
        %%%%%%%%%%%%%%%%%%%%
        if checkSums
            t1=clock;
            %% Check W sums
            sW0=full(sparse(double(crps.tokens.word(1:end-1)),ones(nTokens-1,1),double(x(2:end)==0),nWords,1));
            sW1=full(sparse(double(crps.tokens.word(1:end-1)),ones(nTokens-1,1),double(x(2:end)==1),nWords,1));
            %[sumW0,sumW1,sumW0+sumW1,sW0,sW1,sW0+sW1]
            if any(sW0-sumW0)
                fprintf('\nERROR: sum(sumW0)=%10d, real sum=%10d\n',sum(sumW0),sum(sW0));
                %error('Error in sumW0\n');
            end
            if any(sW1-sumW1)
                fprintf('\nERROR: sum(sumW1)=%10d, real sum=%10d\n',sum(sumW1),sum(sW1));
                %error('Error in sumW0\n');
            end
            
            %% Check T sums
            sTW=full(sparse(double(z),double(crps.tokens.word),...
                            double(x==0),nTopics,nWords));
            sTD=full(sparse(double(z),double(crps.tokens.document),...
                            double(x==0),nTopics,nDocs));
            sT=full(sparse(double(z),ones(nTokens,1),...
                           double(x==0),nTopics,1));
            if any(any(sTW-sumTW))
                fprintf('\nERROR: sum(sumTW)=%10d, real sum=%10d\n',sum(sum(sumTW)),sum(sum(sTW)));
                %error('Error in sumTW\n');
            end
            if any(any(sTD-sumTD))
                fprintf('\nERROR: sum(sumTD)=%10d, real sum=%10d\n',sum(sum(sumTD)),sum(sum(sTD)));
                %error('Error in sumTD\n');
            end
            if any(sT-sumT)
                fprintf('\nERROR: sum(sumT) =%10d, real sum=%10d\n',sum(sum(sumT)),sum(sum(sT)));
                %error('Error in sumT\n');
            end
            fprintf('sum-checks (%.3fsec)... ',etime(clock,t1));
        end

        if isnan(ldacol.logL(iter))
            error('ldacolClassesnew: NaN log-likelihood\n');
        end

        ldacol.percentCollocations(iter)=sum(sumW1)/maxCollocations;
        fprintf('likelihood = %12.1f, # colls = %10d/%d (%6.3f%%)\n',...
                ldacol.logL(iter),sum(sumW1),maxCollocations,100*ldacol.percentCollocations(iter));

        %%%%%%%%%%%%%%%%%%%%
        %% Likelihood plots
        %%%%%%%%%%%%%%%%%%%%
        if mod(iter,20)==0 || j==nIterations(i)
            clf
            subplot(3,1,1)
            plot(1:iter,ldacol.logL(1:iter),'.-');
            axis tight
            grid on
            xlabel('iteration')
            ylabel('log likelihood')
            title(sprintf('current log likelihood = %.1f',ldacol.logL(iter)));
            minIter=200;
            if iter>minIter
                subplot(3,1,2)
                plot(minIter:iter,ldacol.logL(minIter:iter),'.-');
                axis([minIter,iter,min(ldacol.logL(minIter:iter)),max(ldacol.logL(minIter:iter))])
                axis tight
                grid on
                xlabel('iteration')
                ylabel('log likelihood')
            end
            subplot(3,1,3)
            plot(1:iter,100*ldacol.percentCollocations(1:iter),'.-');
            axis tight
            grid on
            xlabel('iteration')
            ylabel('% of collocations')
            title(sprintf('current collocations = %.1f%%',100*ldacol.percentCollocations(iter)));
            drawnow
        end
    
    end % for j=1:nIterations(i)

    crps.tokens.topic=1*z;
    crps.tokens.collocation=1*x;

    %% Include collocations in corpus
    crpsLDA=corpusWithCollocations(crps);
    
    %% Add topics to corpus
    docsLDA=computeTopics(crpsLDA,topicNameLength,nTopics);
    topicsLDA=crpsLDA.topics;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save iterations plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Save iterations plot
    fprintf('\nSaving iterations plot to ''%s'' ... ',convergence)
    t1=clock;
    print(sprintf('%s.eps',convergence),'-depsc2');
    fprintf('done (%.2f sec)\n',etime(clock,t1));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Save topics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Print topics
    fprintf('\nSaving topics to ''%s.tab''\n',topics); 
    t1=clock;
    fid=fopen(sprintf('%s.tab',topics),'w');
    fprintf('\nMost likely words in the topics:\n');
    fprintf(fid,'Iteration %d:\n',iter);
    fprintf(fid,'Most likely words in the topics:\n');
    % compressed format
    fprintf(fid,'topic\tweight in docsLDA\ttopic spread\twords\n');
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
    
    %% Plot topic spreads
    [d,k]=sort(crpsLDA.topics.topicSpread,'descend');
    clf
    subplot(2,1,1)
    plot(1:length(d),crpsLDA.topics.topicSpread,'.')
    grid on
    ylabel('topic spread')
    xlabel('topic number')
    title('Topic spread (0 - much spread, 1 - no spread)');
    subplot(2,1,2)
    plot(1:length(d),d,'.')
    text(1:2:length(d),d(1:2:end),cellstr(num2str((k(1:2:end)))),...
         'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
    text(2:2:length(d),d(2:2:end),cellstr(num2str((k(2:2:end)))),...
         'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8)
    grid on
    title('Sorted topic spread (0 - much spread, 1 - no spread)');
    ylabel('topic spread')
    xlabel('rank order of topic spread')
    legend('topic spread')
    
    fprintf('\nSaving topic spreads plot to ''%s.eps'' ... ',topics)
    t1=clock;
    print(sprintf('%s.eps',topics),'-depsc2');
    fprintf('done (%.2f sec)\n',etime(clock,t1));
end

%% Done
fprintf('Done ldacol (%.2f sec)\n',etime(clock,t0));

whos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);





