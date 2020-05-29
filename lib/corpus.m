classdef corpus < handle
% Corpus of documents represented by
% tokens - table where each row corresponds to a token (i.e., non-unique word)
%          in the corpus, with the following columns:
%           . word        - categorical string variable with the word (in lowercase)
%           . document    - categorical string variable with the name of the document
%           . index       - integer position of the token within the document (from 1)
%           . collocationAllowed - boolean variable indicating if the token can form
%                                  a collocation with the *previous* one
%          The following (optional) columns are typically created from an LDA analysis
%           . topic       - topic of the token
%           . collocation - boolean variable indicating if the token actually forms
%                           a collocation with the *previous* one
%           . wordPairs   - word-pair counts: wordPairs(i) counts how often word(i)
%                           is followed by the word(i+1).
% topics - table where each row corresponds to a topic, with the following columns:
%           . wordWeights - array of the probability of each word for the topic 
%                           (one probability for each word in the dictionary)
%           . topicName   - string array with the top-weighted words in the topic
%           . averageDocWeights - typically topic weights (probability distribution over
%                           topics) for a document that is representative of the topic, 
%                           computed using
%                             averageDocWeights for topic t = 
%                                = Sum_{d=1}^{# docs} weight(document d,topic t) *  
%                                                     (vector of weights for document d)
%                                   / [ Sum_{d=1}^{# docs} weight(document d,topic t) ]
%                           where the sum is over all documents and the weights 
%                                weight(document d,topic t) 
%                           corresponds to the probability of the topic t  in the document d.
%           . topicSpread - value in the [0,1] interval that quantifies how much the documents
%                           with high values for the topic are clustered around the averageDocWeights
%                           for that topic:
%                              1 - the documents with high values for the topic are 
%                                  tightly clustered around averageDocWeights
%                              0 - the documents are not much clustered
%                           This value is computed as follows:
%                             topicSpread = Sum_{d=1}^{#docs}  weight(document d,topic t) * 
%                                              dist(document d, averageDocWeights for topic t)
%                                   / [ Sum_{d=1}^{# docs} weight(document d,topic t) ]
%                           where the sum is over all documents; the weights 
%                                weight(document d,topic t) 
%                           correspond to the probability of the topic t  in the document d;
%                           and 
%                                dist(document d, averageDocWeights for topic t)
%                           is the cosine distance between document d and averageDocWeights for topic t.
% Copyright (C) 2013-2014  Stacy & Joao Hespanha
    
    properties
        tokens 
        topics
    end
    
    methods (Static)
        
        function crps=loadobj(s)
            crps=corpus();
            crps.tokens=s.tokens;
            crps.topics=s.topics;
        end
    
    end
    
    methods 
        
        function s=saveobj(crps)
            s=struct('tokens',crps.tokens,'topics',crps.topics);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% corpus - object creation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function crps=corpus(word,document,index,collocationAllowed)
        % crps=corpus(word,document,index,collocationAllowed)
        % Create object of class corpus given all the (non-optional) variables for the table.
            
            if nargin==0
                return
            end
            
            %% Create corpus table
            fprintf('  Creating corpus table...');
            t1=clock;
            crps.tokens=table(categorical(word),categorical(document),...
                             uint64(index),uint32(collocationAllowed),...
                             'VariableNames',{'word','document','index','collocationAllowed'});
            crps.topics=table();
            fprintf('done (%d tokens, %d words, %d documents, %.3f sec)\n',...
                    size(crps.tokens,1),length(categories(crps.tokens.word)),...
                    length(categories(crps.tokens.document)),etime(clock,t1));
            
            %% update word pairs
            cleanup(crps)

        end
        
        function cleanup(crps)
        % cleanup(crps)
        % Remove all unused values of the categorical variables
            
            fprintf('  Cleaning up categories (%d tokens, %d words, %d documents) ... ',...
                    size(crps.tokens,1),length(categories(crps.tokens.word)),...
                    length(categories(crps.tokens.document)));
            t1=clock;
            
            % rebuild categories
            crps.tokens.word=removecats(crps.tokens.word);
            crps.tokens.document=removecats(crps.tokens.document);
            
            fprintf('(%d tokens, %d words, %d documents, %.3f sec)\n',...
                    size(crps.tokens,1),length(categories(crps.tokens.word)),...
                    length(categories(crps.tokens.document)),etime(clock,t1));
        end
        
        function docs=computeTopics(crps,topicNameLength,nTopics)
        % docs=computeTopics(crps,docsTable)
        % 1) Computes topics weights for each document in corpus and returns the 
        %    weights in a table of documents
        % 2) Computes word weights for each topic and adds the weights to the
        %    table of topic
        % 3) Computes the top words for each topic and adds them to the
        %    table of topic
            
            %% debug
            % crps.tokens.topic=ceil(10*rand(50,1));
            % crps.tokens.document=ceil(3*rand(50,1));
            % [crps.tokens.document,crps.tokens.topic]
            
            %% Compute topic weights per document and add them to the document's table
            nDocs=length(categories(crps.tokens.document));
            docs=table();
            docs.topicWeights=sparse(double(crps.tokens.document),...
                                     double(crps.tokens.topic),...
                                     ones(size(crps.tokens,1),1),nDocs,nTopics);
            docs.nTokens=sum(docs.topicWeights,2);
            docs.topicWeights=docs.topicWeights./...
                repmat(docs.nTokens,1,size(docs.topicWeights,2));
            docs.primaryKey=categorical(categories(crps.tokens.document));
            
            % sort topics by total importance in corpus
            [~,k]=sort(sum(docs.topicWeights,1),2,'descend');
            docs.topicWeights=docs.topicWeights(:,k);
            ik=zeros(length(k),1);
            ik(k)=1:length(k);
            crps.tokens.topic=ik(crps.tokens.topic);
            
            %% Compute word weights per topic and create topics table
            % computed in the transposed for to speedup the loop
            words=categories(crps.tokens.word);
            wordWeights=sparse(double(crps.tokens.word),...
                               double(crps.tokens.topic),...
                               ones(size(crps.tokens,1),1),length(words),nTopics);
            wordWeights=wordWeights./...
                repmat(sum(wordWeights,1),size(wordWeights,1),1);

            % Topic centers of mass  -- weighted average of distribuitions, so also a distribution
            % CM(topic,:) = sum_doc weight(doc,:) weight(doc,topic) / sum_doc weight(doc,topic)
            crps.topics.averageDocWeights=full(diag(1./sum(docs.topicWeights,1))*(docs.topicWeights'*docs.topicWeights));
            % Cosine distances to topics CM
            d2CMs=docs.topicWeights*crps.topics.averageDocWeights';
            % Weighted cosine distances to topic CMs
            crps.topics.topicSpread=full(sum(d2CMs.*docs.topicWeights,1)./sum(docs.topicWeights,1))'; 

            %% Compute topic names
            topicNames=cell(size(wordWeights,2),1);
            for k=1:size(wordWeights,2);
                [w,word]=sort(wordWeights(:,k),1,'descend');
                topicNames{k}=sprintf('%s ',words{word(1:min(end,topicNameLength))});
            end            
            crps.topics.wordWeights=wordWeights'; 
            crps.topics.topicName=topicNames;
        end

        function newCrps=copy(crps)
        % newCrps=copy(crps)
        % Returns a copy of the (handle) object crps. 
        % Note that a simple assignment does not work for handle objects
            newCrps.tokens=crps.tokens;
            newCrps.topics=crps.topics;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% methods to update the corpus
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function removeDocuments(crps,docs2remove)
        % removeDocuments(crps,docs2remove)
        % Remove set of documents
            fprintf('  Removing %d documents...',length(docs2remove));
            t1=clock;
            k=ismember(crps.tokens.document,docs2remove);
            % remove tokens
            crps.tokens(k,:)=[];
            fprintf('done (%d tokens, %d words, %d documents, %.3f sec)\n',...
                    size(crps.tokens,1),length(categories(crps.tokens.word)),...
                    length(categories(crps.tokens.document)),etime(clock,t1));

            cleanup(crps);
        end
        
        function removeWords(crps,words2remove,breakCollocations)
        % removeWords(crps,words2remove,breakCollocations)
        % Remove set of words. When breakCollocations is true, collocations
        % are not allowed accross words removed, otherwise they are allowed to
        % skip the words removed.
           
            fprintf('  Removing %d words (breakCollocations=%)...',...
                    length(words2remove),breakCollocations);
            t1=clock;
            k=find(ismember(crps.tokens.word,words2remove));
            if breakCollocations
                % no colocation across removed words: every word after
                % a removed word cannot be allowed a collocation with
                % the previous word
                crps.tokens.collocationAllowed(k(k<size(crps.tokens,1))+1)=0;
            else
                % allow collocations to skip removed words: only "cut"
                % collocations if the removed word was not allowed a
                % collocation with previous one.
                kS=k(crps.tokens.collocationAllowed(k)==0); 
                while ~isempty(kS)
                    % last word in corpus does not matter
                    kS(kS==size(crps.tokens,1))=[]; 
                    kS=kS+1;
                    % if stopword could not connect to previous, then this is passed to next word
                    crps.tokens.collocationAllowed(kS)=0;    
                    kS=kS(ismember(kS,k));
                end
            end
            % remove tokens
            crps.tokens(k,:)=[];
            fprintf('done (%d tokens, %d words, %d documents, %.3f sec)\n',...
                    size(crps.tokens,1),length(categories(crps.tokens.word)),...
                    length(categories(crps.tokens.document)),etime(clock,t1));

            cleanup(crps);
        end
        
        function newCrps=corpusWithCollocations(crps)
        % newCrps=corpusWithCollocations(crps)
        % Created a new corpus for which words with collocations have been
        % compressed into "meta-words".
            
            t1=clock;
            fprintf('  Creating collocation corpus... ');
            x=crps.tokens.collocation;
            x(end+1)=0;
            wstart=find(x==0);
            nTokens=length(wstart)-1;
            if 0
                %% ATTENTION: This is very slow. Perhpas it can be improved.
                words=cell(nTokens,1);
                parfor i=1:nTokens
                    if mod(i,10000)==0
                        t = getCurrentTask();
                        fprintf('worker %d: %d/%d (%.0fsec)\n',t.ID,i,nTokens,etime(clock,t1));
                    end
                    ww=wstart(i):wstart(i+1)-1;
                    ww=cellstr(crps.tokens.word(ww));
                    ww=sprintf('%s_',ww{:});
                    ww(end)=[];
                    words{i}=ww;
                end
            else
                % x=round(rand(size(crps.tokens,1)+1,1));
                words=char(crps.tokens.word);
                words(:,end+1)=char(0);
                words(find(x(2:end)),end)='_';
                words=words';
                words=words(:)';
                words(words==' ')=[];
                words(end)=[];
                [words,~]=regexp(words,char(0),'split','match');
                words=words';
                clear x;
            end
             
            % save to table
            wstart(end)=[];
            newCrps=corpus(words,crps.tokens.document(wstart),crps.tokens.index(wstart),zeros(nTokens,1,'uint32'));
            newCrps.tokens.topic=crps.tokens.topic(wstart); %% what if changes in middle of collocation?
            
            %cleanup(newCrps)
            fprintf('%d collocation words (from %d simple words), %d collocation tokens (from %d simple tokens) done (%.3fsec)\n',...
                    length(categories(newCrps.tokens.word)),length(categories(crps.tokens.word)),...
                    size(newCrps.tokens,1),size(crps.tokens,1),etime(clock,t1));
            
        end     
        
    end
end