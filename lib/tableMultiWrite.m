function varargout=tableMultiWrite(varargin);
% To get help, type tableMultiWrite('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This scripts writes a table to a file.'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','outputTable',...
    'DefaultValue',getFromPedigree(),...
    'Description', {
        'Path and basename for the output files, used for write access (output)'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','tbl',...
    'Description', {
        'Table to be written.'
                   });

declareParameter(...
    'VariableName','variableNames',...
    'DefaultValue',[],...
    'Description', {
        'Names of the columns from the ''baseTable'' to be included in'
        'the joint table. When empty (i.e.,[]) all columns are included.'
                   });
    
declareParameter(...
    'VariableName','firstRow',...
    'DefaultValue',1,...
    'Description', {
        'Index of the first row to write.'
                   });

declareParameter(...
    'VariableName','lastRow',...
    'DefaultValue',inf,...
    'Description', {
        'Index of the last row to write. If the table does not have this'
        'many rows, this value is replaced by the number of rows of the table.'
                   });

declareParameter(...
    'VariableName','outputFormats',...
    'Description', {
        'Format of the output file:'
        '  tab - tab separated file with the variable names in the 1st row.'
        '  csv - comma separated file with the variable names in the 1st row.'
        '  mat6 - version 6 matlab file with the table'
        '         (file will contain a single variable, named ''tbl'')'
        '  mat7 - version 7 matlab file with the table'
        '         (file will contain a single variable, named ''tbl'')'
        '  mat7.3 - version 7.3 matlab file with the table'
        '         (file will contain a single variable, named ''tbl'')'
        'When a cell array of strings is provided, the file is saved in multiple formats,'
        'all with the same base name, but with different extensions.'
                   });
    
declareParameter(...
    'VariableName','quoteCells',...
    'DefaultValue',true,...
    'AdmissibleValues',{true,false},...
    'Description', {
        'When true, the content of every cell is enclosed in double quotes ("...").'
        'In this case, any double quote characters that appear as part of a ccll'
        'are replaced by two double quote characters.'
                   });
declareParameter(...
    'VariableName','endOfLineWord',...
    'DefaultValue',' LiNeBrEaK ',...
    'Description', {
        'Word used to replace new lines in the tab and csv formats.'
        'Only done when quoteCells=false.'
                   });

declareParameter(...
    'VariableName','tabWord',...
    'DefaultValue',' ',...
    'Description', {
        'Word used to replace tabs in the tab format.'
        'Only done when quoteCells=false.'
                   });

declareParameter(...
    'VariableName','quoteWord',...
    'DefaultValue',' ',...
    'Description', {
        'Word used to replace quotes in the tab format.'
        'Only done when quoteCells=false.'
                   });

declareParameter(...
    'VariableName','commaWord',...
    'DefaultValue',' ',...
    'Description', {
        'Word used to replace commas in the csv format.'
        'Only done when quoteCells=false.'
                   });

%% xml, tab, csv
declareParameter(...
    'VariableName','encoding',...
    'DefaultValue','UTF-8',...
    'Description', {
        'Encoding default scheme assumed for the input file.'
        'This parameter is only used for the text formats:'
        '  {''xml'', ''csv'', ''tab''} '
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','tbl',...
    'Description', {
        'Table that was written.'
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

if ischar(outputFormats)
    outputFormats={outputFormats};
end

%% remove extension
[outputPath,outputName]=fileparts(outputTable);
outputTable=fullfile(outputPath,outputName);

fprintf('tableMultiWrite: %d formats\n',length(outputFormats))
t0=clock;

fprintf('  selecting rows in tbl(%dx%d): ',size(tbl));
t1=clock;
tbl=tbl(firstRow:min(end,lastRow),:);
fprintf('  done tbl(%dx%d) (%.3f sec)\n',size(tbl),etime(clock,t1));

if ~isempty(variableNames)
    if ~ismember('primaryKey',variableNames)
        variableNames={'primaryKey',variableNames{:}};
    end

    fprintf('  selecting %d columns in tbl(%dx%d): ',length(variableNames),size(tbl));
    fprintf('%s ',variableNames{:});
    t1=clock;
    VNs=intersect(variableNames,tbl.Properties.VariableNames);
    if length(VNs)<length(variableNames)
        ignored=setdiff(variableNames,VNs);
        fprintf('\n    WARNING: %d fields do not exist: ',length(ignored));
        fprintf('%s ',ignored{:});
        fprintf('\n');
    end
    tbl=tbl(:,VNs);
    fprintf('  done tbl(%dx%d) (%.3f sec)\n',size(tbl),etime(clock,t1));
end

%% Save
for thisFormat=1:length(outputFormats)
    fprintf('  Writing %d/%d %s file %-60s... ',...
            thisFormat,length(outputFormats),...
            outputFormats{thisFormat},outputTable);
    t1=clock;
    switch (outputFormats{thisFormat})
      case {'tab'}
        % 'quote' tabs
        if quoteCells
            tbl1=tableRegexprep(tbl,tbl.Properties.VariableNames,...
                                {'\t','"','[\n\r]'},{tabWord,quoteWord,endOfLineWord});
        else
            tbl1=tbl;
        end
        % expand categorical to strings for speed
        tbl1=expandCategorical(tbl1);
        
        if exist('get_param','builtin') && exist('set_param','builtin')
            old=get_param(0,'CharacterEncoding');
            set_param(0,'CharacterEncoding',encoding);
        else
            fprintf('tableMultiWrite: ''get_param'' not available, unable to guarantee encoding %s\n',encoding);
        end
        %whos
        %summary(tbl1)
        writetable(tbl1,sprintf('%s.tab',outputTable),'FileType','text',...
                   'WriteVariableNames',true,...
                   'WriteRowNames',false,...
                   'QuoteStrings',quoteCells,...
                   'Delimiter','\t');
        if exist('get_param','builtin') && exist('set_param','builtin')
            set_param(0,'CharacterEncoding',old);
        end
        clear tbl1
        
      case {'csv'}
        % remove tabs
        if quoteCells
            tbl1=tableRegexprep(tbl,tbl.Properties.VariableNames,...
                                {',','"','[\n\r]'},{commaWord,quoteWord,endOfLineWord});
        else
            tbl1=tbl;
        end
        % expand categorical to strings for speed
        tbl1=expandCategorical(tbl1);
        
        if exist('get_param','builtin') && exist('set_param','builtin')
            old=get_param(0,'CharacterEncoding');
            set_param(0,'CharacterEncoding',encoding);
        else
            fprintf('tableMultiWrite: ''get_param'' not available, unable to guarantee encoding %s\n',encoding);
        end
        writetable(tbl1,sprintf('%s.csv',outputTable),'FileType','text',...
                   'WriteVariableNames',true,...
                   'WriteRowNames',false,...
                   'QuoteStrings',quoteCells,...
                   'Delimiter',',');
        if exist('get_param','builtin') && exist('set_param','builtin')
            set_param(0,'CharacterEncoding',old);
        end
        clear tbl1
        
      case {'mat6'}
        save(sprintf('%s.mat',outputTable),'tbl','-v6'); % backward compatible, can be read into R using matlab.R package
      case {'mat7'}
        save(sprintf('%s.mat',outputTable),'tbl','-v7'); % much faster
      case {'mat7.3'}
        save(sprintf('%s.mat',outputTable),'tbl','-v7.3'); % slower but needed for large files
      otherwise
        error('tableMultiread: unkown outputFormat ''%s''\n',outputFormats{thisFormat});
    end
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
end
    
fprintf('done tableMultiWrite (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);

end

function tbl=expandCategorical(tbl)

for i=1:size(tbl,2)
    if strcmp(class(tbl{:,i}),'categorical')
        fprintf('tableMultiWrite: expanding column %d from class=%s ',i,class(tbl{:,i}));
        % must actually remove column and replace by new one since
        % otherwise type does ot change from categorical
        x=cellstr(tbl{:,i});
        xn=tbl.Properties.VariableNames{i};
        tbl(:,i)=[];
        tbl(:,end+1)=x;
        tbl.Properties.VariableNames{end}=xn;
        fprintf('to class=%s\n',class(tbl{:,i}));
    end
end

end
