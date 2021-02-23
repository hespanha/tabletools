function varargout=tableMultiRead(varargin);
% To get help, type tableMultiRead('help')
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

declareParameter(...
    'Help', {
        'This script contructs a table class object from multiple files.'
            });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','inputTables',...
    'Description', {
        'Directory or wildcard of input files, used for read access (input)'
                   });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareParameter(...
    'VariableName','inputFormat',...
    'AdmissibleValues',{'xml','tab','csv','mat','xlsx'},...
    'Description', {
        'Format of the input files:'
        '  xml - xml file with each row encapsulated by a given tag (see rowTag),'
        '        inside which the different variables are encapsulated'
        '        within tags (see variableTags).'
        '  tab - tab separated file with the variable names in the 1st row.'
        '  csv - comma separated file with the variable names in the 1st row.'
        '  mat - matlab file with the table (must contain a single variable)'
        'If needed, the variable names are transformed to make then valid'
        'Matlab identifiers (a warning is generated if a transformation is needed).'
                   });
    
declareParameter(...
    'VariableName','maxRows',...
    'DefaultValue',inf,...
    'Description', {
        'Maximum number of rows to read. When equal to inf, all rows should be read.'
                   });
    
declareParameter(...
    'VariableName','tableDescription',...
    'DefaultValue',{},...
    'Description', {
        'String with a description of the table.'
                   });
    
declareParameter(...
    'VariableName','variableDescriptions',...
    'DefaultValue',{},...
    'Description', {
        'Nx2 cell array of strings, where variableDescriptions{:,2} containst the cell'
        'descriptions of the variables in variableDescriptions{:,1}.'
                   });

declareParameter(...
    'VariableName','primaryKey',...
    'DefaultValue','',...
    'Description', {
        'Name of the table variable to be used as the primary key'
        'for linking tables. A new copy of the variable will be added'
        'with the name ''primaryKey''.'
        'When empty, the primary key is automatically constructed as'
        'a string of the form ''%d:%s'' where %s is the name of the'
        'file (without a path) from which the row was'
        'read and %d the number of the record within the file.'
        ' '
        'The primary key will be checked for repetitions and an error'
        'is generated if it has repeated values.'
        ' '
        'An error will also be generated if the table already has'
        'a column called ''primaryKey''.'
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

%% tab, csv specific
declareParameter(...
    'VariableName','numberVariables',...
    'DefaultValue',{},...
    'Description', {
        'Array of strings with names of variables (i.e., column headers) that'
        ' should be converted to numeric values.'
                   });

declareParameter(...
    'VariableName','ignoreVariables',...
    'DefaultValue',{},...
    'Description', {
        'Array of strings with names of variables (i.e., column headers) that'
        'should be ignored.'
                   });


%% xml specific
declareParameter(...
    'VariableName','rowTag',...
    'DefaultValue','row',...
    'Description', {
        'String with the tag that encapsulates each table row.'
        'This parameter is only used for the ''xml'' format.'
                   });
declareParameter(...
    'VariableName','stringTags',...
    'DefaultValue','{str,date}',...
    'Description', {
        'String with the tag that encapsulates a string variable.'
        'Multiple tags can be provided as a string array.'
        'The name of the variable is taken from the attribute "name" for these tags.'
        'This parameter is only used for the ''xml'' format.'
                   });
declareParameter(...
    'VariableName','numberTags',...
    'DefaultValue','{long}',...
    'Description', {
        'String with the tag that encapsulates a numerical variable.'
        'Multiple tags can be provided as a string array.'
        'The name of the variable is taken from the attribute "name" for these tags.'
        'This parameter is only used for the ''xml'' format.'
                   });
declareParameter(...
    'VariableName','arrayTags',...
    'DefaultValue','arr',...
    'Description', {
        'String with the tag that encapsulates an array of strings.'
        'The name of the variable is taken from the attribute "name" for these tags.'
        'This parameter is only used for the ''xml'' format.'
                   });
declareParameter(...
    'VariableName','arraySeparator',...
    'DefaultValue',' | ',...
    'Description', {
        'String used to separate elements of an array of strings.'
        'This parameter is only used for the ''xml'' format.'
                   });

%%% Checks
declareParameter(...
    'VariableName','skipCheckClasses',...
    'DefaultValue',false,...
    'AdmissibleValues',{false,true},...
    'Description', {
        'When true, skips checking whether table columns are of a consistent types.'
        });

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

declareOutput(...
    'VariableName','tbl',...
    'Description', {
        'Table that was read.'
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

if ischar(numberVariables)
    numberVariables={numberVariables};
end
if ischar(stringTags)
    stringTags={stringTags};
end
if ischar(numberTags)
    numberTags={numberTags};
end
if ischar(arrayTags)
    arrayTags={arrayTags};
end

if isdir(inputTables)
    inputTables=fullfile(inputTables,'*.*');
end
files=dir(inputTables);
inputPath=fileparts(inputTables);

fprintf('tableMultiRead: %d files\n',length(files))
if isempty(files)
    inputTables
    error('no files to read');
end
t0=clock;

allVariableNames=cell(0,1);
tbl=cell(0,0);
rowNames=cell(0,1);
for thisFile=1:length(files)
    fprintf('    %d/%d %s file %-60s (%6.0fKB)... ',...
            thisFile,length(files),...
            inputFormat,files(thisFile).name,files(thisFile).bytes/1024);
    if strncmp(files(thisFile).name,'._',2)
        fprintf(' ignoring OSX hidden file\n');
        continue;
    end

    thisFilename=fullfile(inputPath,files(thisFile).name);
    [~,shortName]=fileparts(thisFilename);
    
    if size(tbl,1)>=maxRows
        tbl=tbl(1:maxRows,:);
        rowNames=rowNames(1:maxRows,1);
        break;
    end
    
    switch (inputFormat)
      case {'tab','csv','mat','xlsx'}
        
        t1=clock;
        switch (inputFormat)
          case {'mat'}
            mfile=matfile(thisFilename);
            w=whos(mfile);
            for i=1:length(w)
                if strcmp(w(i).class,'table')
                    if istable(tbl)
                        error('tableMultiRead: file %s has multiple table variables\n',thisFilename);
                    else
                        t=mfile.(w(i).name);
                    end
                end
            end
            clear mfile
          case {'tab'}
            t=readTextFile(thisFilename,'\t',encoding);
          case {'csv'}
            opts = detectImportOptions(thisFilename);
            opts.Encoding=encoding;
            t=readtable(thisFilename,opts);
          case {'xlsx'}
            opts = detectImportOptions(thisFilename);
            t=readtable(thisFilename,opts);
        end            
        
        for v=1:length(ignoreVariables)
            if ismember(ignoreVariables{v},t.Properties.VariableNames)
                fprintf('~%s,',ignoreVariables{v});
                t.(ignoreVariables{v})=[];
            end
        end

        variableNames=t.Properties.VariableNames;
        while 1 % this loop should only repeat once
            [exists,locb]=ismember(variableNames,allVariableNames);
            if all(exists)
                empt=cell(size(t,1),size(tbl,2));
                empt(:)={''};
                tbl(end+1:end+size(t,1),:)=empt;
                clear empt
                tbl(end-size(t,1)+1:end,locb)=table2cell(t);
                break;
            else
                %fprintf('enlarging columns');
                allVariableNames(end+1:end+sum(~exists),1)=variableNames(~exists);
                empt=cell(size(tbl,1),length(allVariableNames)-size(tbl,2));
                empt(:)={''};
                tbl(:,end+1:length(allVariableNames))=empt;
                clear empt
            end
        end
        thisNrows=size(t,1);
        fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
        
      case {'xml'}
        
        t1=clock;
        if exist('get_param','builtin') && exist('set_param','builtin')
            old=get_param(0,'CharacterEncoding');
            set_param(0,'CharacterEncoding',encoding);
        else
            fprintf('tableMultiRead: ''get_param'' not available, unable to guarantee encoding %s\n',encoding);
        end
        xml=xmlread(thisFilename);
        if exist('get_param','builtin') && exist('set_param','builtin')
            set_param(0,'CharacterEncoding',old);
        end
        fprintf('done (%.3f sec)\n',etime(clock,t1));
        
        allRows = xml.getElementsByTagName(rowTag);
        
        fprintf('    parsing xml structure... ');
        t1=clock;
        firstRow=size(tbl,1)+1;
        % "reserve" memory for data
        if size(tbl,2)>0
            %fprintf('enlarging rows');
            empt=cell(allRows.getLength,size(tbl,2));
            empt(:)={''};
            tbl(firstRow:firstRow+allRows.getLength-1,1:end)=empt;
            clear empt
        else
            tbl=cell(allRows.getLength,0);
        end
        %tbl(end,:)
        %% loop over rows (defined by ''rowTag'')
        for row=0:allRows.getLength-1
            if firstRow+row>maxRows
                break;
            end
            
            if mod(row,250)==0
                fprintf('Row %d/%d ',firstRow+row,size(tbl,1));
            end
            
            thisRow=allRows.item(row);
            
            thisRow=thisRow.getChildNodes;
            
            variableNames=cell(1,thisRow.getLength);
            values=cell(1,thisRow.getLength);
            n=0;
            for tag=0:thisRow.getLength-1
                %% loop over items with the row
                thisTag=thisRow.item(tag);
                
                if isempty(get(thisTag,'Attributes'))
                    continue
                end
                
                %% found a variable Tag
                thisAttribute=get(thisTag,'Attributes');
                if thisAttribute.getLength>1
                    continue
                end
                %% get attribute name (to be used as the variable name)
                field=char(thisAttribute.item(0));
                if ~strncmp(field,'name="',6) || field(end)~='"'
                    continue
                end

                n=n+1;
                variableNames{n}=field(7:end-1);
                switch (char(thisTag.getTagName))
                  case stringTags
                    values{n}=char(thisTag.getTextContent);
                    %values{n}=utf2html(values{n});
                    %fprintf('%s(%s) = %s (%s)\n',variableNames{n},char(thisTag.getTagName),values{n},class(values{n}));
                  case numberTags
                    values{n}=str2num(thisTag.getTextContent);
                    %fprintf('%s(%s) = %d (%s)\n',variableNames{n},char(thisTag.getTagName),values{n},class(values{n}));
                  case arrayTags
                    values{n}=char(thisTag.item(0).getTextContent);
                    for i=1:thisTag.getLength-1
                        values{n}=[values{n},arraySeparator,char(thisTag.item(i).getTextContent)];
                    end
                    %values{n}=utf2html(values{n});
                    %fprintf('%s(%s) = %s (%s)\n',variableNames{n},char(thisTag.getTagName),values{n},class(values{n}));
                  otherwise
                    n=n-1;
                    continue
                end
            end
            variableNames=variableNames(1:n);
            values=values(1:n);

            % % find 2-byte encodings
            % k=find(cellfun(@(v)ischar(v) && any(double(v)>255),values));
            % if ~isempty(k)
            %     for l=1:length(k)
            %         j=find(double(values{k(l)})>255);
            %         fprintf('2-byte encoding in table(%d,''%s'') "',firstRow+row,variableNames{k(l)});
            %         fprintf('%s',utf2html(values{k(l)}(j)));
            %         fprintf('"\n');
            %     end
            % end
            
            while 1  % this loop should only repeat once
                [exists,locb]=ismember(variableNames,allVariableNames);
                if all(exists)
                    tbl(firstRow+row,locb)=values;
                    break;
                else
                    %fprintf('enlarging columns');
                    allVariableNames(end+1:end+sum(~exists),1)=variableNames(~exists);
                    empt=cell(size(tbl,1),length(allVariableNames)-size(tbl,2));
                    empt(:)={''};
                    tbl(:,end+1:length(allVariableNames))=empt;
                end
            end
        end
        tbl(end,:)
        thisNrows=allRows.getLength;
        fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
      otherwise
        error('tableMultiRead: unkown inputFormat ''%s''\n',inputFormat);
    end % switch
    rowNames(end+1:size(tbl,1),1)=strcat(cellstr(num2str((1:thisNrows)','%-d')),cellstr(repmat([':',shortName],thisNrows,1)));
end % for thisFile=1:length(files)
clear t

%% Construct table
fprintf('  converting to table...  ');
t1=clock;
if ~isempty(primaryKey)
    rowNames=tbl(:,strcmp(allVariableNames,primaryKey));
end
%tbl=cell2table(tbl,'VariableNames',allVariableNames,'RowNames',rowNames);
tbl=cell2table(tbl,'VariableNames',allVariableNames);
fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));


if ismember('primaryKey',allVariableNames)
    error('  table already has a column named ''primaryKey''\n');
else
    fprintf('  adding primary key... ');
    t1=clock;
    if length(rowNames)~=length(unique(rowNames))
        error('\ntableMultiRead: primary key has repeated entries, length(primaryKey)=%d, length(unique(primaryKey))=%d\n',length(rowNames),length(unique(rowNames)));
    end
    fprintf('  has no repetition... ');
    tbl.primaryKey=categorical(rowNames);
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
    allVariableNames=[{'primaryKey'};allVariableNames];
end

%% Assign cell descriptions
% start with descriptions = variable names
tbl.Properties.VariableDescriptions=allVariableNames;
% include descriptions provided by the user
if size(variableDescriptions,1)>0
    k=ismember(variableDescriptions(:,1),tbl.Properties.VariableNames);
    if any(~k)
        fprintf('tableMultiRead: ignoring descriptions of unexistent fields ');
        fprintf('''%s'' ',variableDescriptions{~k,1});
        fprintf('\n');
    end
    if any(k)
        tbl.Properties.VariableDescriptions(variableDescriptions(k,1))=variableDescriptions(k,2);
    end
end

if length(numberVariables)>0
    fprintf('  convertions to double... ');
    for i=1:length(numberVariables)
        if ismember(numberVariables{i},tbl.Properties.VariableNames)
            fprintf('%d/%d(%s) ',i,length(numberVariables),numberVariables{i});
            if iscell(tbl.(numberVariables{i}))
                k=cellfun(@ischar,tbl.(numberVariables{i}));
                v=nan(size(tbl,1),1);
                v(k)=str2double(cellstr(tbl.(numberVariables{i})(k)));
                v(~k)=[tbl.(numberVariables{i}){~k}];
                tbl.(numberVariables{i})=v;
            else
                if isnumeric(tbl.(numberVariables{i}))
                    fprintf('(num2double) ');
                    tbl.(numberVariables{i})=double(tbl.(numberVariables{i}));
                else
                    fprintf('(str2double) ');
                    x=cellfun(@class,tbl.(numberVariables{i}),'UniformOutput',false);
                    unique(x);
                    
                    tbl.(numberVariables{i})=str2double(tbl.(numberVariables{i}));
                end
            end
        else
            tbl.Properties.VariableNames
            fprintf('tableMultiRead: numberVariable ''%s'' not a table variable\n',numberVariables{i});
        end
    end
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t1));
end
if ~skipCheckClasses
    %% Check types
    fprintf('ANY TYPE MISMATCHES SHOULD NOT OCCUR, THEY INDICATE THAT THE CODE NEEDS CHECKING\n');
    tbl=tableCheckClasses(tbl); % not really needed here
end

%% Assign table description
if ~isempty(tableDescription)
    tbl.Properties.Description=tableDescription;
end

fprintf('done tableMultiRead (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t0));

%whos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout=setOutputs(nargout,params);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function tbl=readTextFile(filename,delimiter,encoding)

verboseLevel=0;

if verboseLevel>0
    fprintf('\n      readTextFile:... ');
end
t0=clock;
fid=fopen(filename,'rt','n',encoding);

if (fid<0) 
    error('file not found: "%s"',filename);
end

variableNames=fgets(fid);
[variableNames,~]=regexp(variableNames,delimiter,'split','match');
nVars=length(variableNames);
if verboseLevel>1
    fprintf('VariableNames\n');
    variableNames
end
if verboseLevel>0
    fprintf('%d headers... ',nVars);
end

% read the whole file
allFile=fread(fid,inf,'char=>char')';
fclose(fid);

if verboseLevel>0
    fprintf('%d bytes... ',length(allFile));
end

% break into lines
[allFile,~]=regexp(allFile,'[\n\r]*','split','match');
if verboseLevel>0
    fprintf('%d lines... ',length(allFile));
end
allFile=allFile';
if isempty(allFile{end})
    allFile(end)=[];
end

% break lines into variables
tbl=cell(length(allFile),nVars);
err=0;
for i=1:length(allFile)
    if mod(i,50000)==0
        fprintf('%d/%d... ',i,length(allFile));
    end
    [vars,~]=regexp(allFile{i},delimiter,'split','match');
    if length(vars)>nVars && all(cellfun('isempty',vars(nVars+1:end)))
       % extra empty cells
       if verboseLevel>0
         fprintf('\n\ttableMultiRead: ignoring %d empty cells at the end of line\n',length(vars)-nVars);
         if verboseLevel>1
            disp(vars(nVars+1:end));
         end
       end
       vars(nVars+1:end)=[];
    end	
    if length(vars)~=nVars
        fprintf('\nLine %d:',i);
        allFile{i}
        fprintf('vars{1:%d}=',length(vars));
        disp(vars)
        fprintf('\n\ttableMultiRead: line %d has %d variables (%d expected)\n',i+1,length(vars),nVars);
        for j=1:length(vars)
            if j<=nVars
                fprintf('\t\tvar{%d} (%s) = ''%s''\n',j,variableNames{j},vars{j});
            else
                fprintf('\t\tvar{%d} = ''%s''\n',j,vars{j});
            end
        end
        err=err+1;
    else
        tbl(i,:)=vars;
    end
end
if err>0
    error('\ntableMultiRead: wrong number of fields in %d lines\n',err);
end
if any(cellfun('isempty',variableNames))
    fprintf('\nVariableNames:\n');
    disp(variableNames')
    fprintf('tableMultiRead: some of the variable names are empty\n');
end
for i=1:length(variableNames)
    if verLessThan('matlab','8.3')
        s=regexprep(variableNames{i},'^/','');  % remove leading /
        s=regexprep(s,'[^A-Za-z0-9_]','_');     % remove characters not allowed
        s=regexprep(s,'^([0-9])','_$1');        % prevent leading numeric
        s=regexprep(s,'^_','x_');               % prevent leading _
    else
        s=matlab.lang.makeValidName(variableNames{i});
    end
    if ~strcmp(s,variableNames{i})
        fprintf('\ntableMultiRead: changed variableName from "%s" to "%s"  ',variableNames{i},s);
        variableNames{i}=s;
    end
end
tbl=cell2table(tbl,'VariableNames',variableNames);
if verboseLevel>0
    fprintf('done (%.3fsec)\n',etime(clock,t0));
end
end