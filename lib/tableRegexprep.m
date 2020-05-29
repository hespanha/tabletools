function tbl=tableRegexprep(tbl,fieldNames,expression,repl);
% tbl=tableRegexprep(tbl,fieldNames,expression,replace)
% 
% Returns a new table for which the fields in fieldNames have been
% processed by performing a regular expression replacement on their
% values.
%
% See regexprep, regarding the meaning of expression and replace.
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

if size(tbl,1)==0
    return;
end

if ~iscell(fieldNames)
    fieldNames={fieldNames};
end

%fprintf('tableRegexprep... ');

% Replace
%t0=clock;

if 1
    % Restrict attention to string fields
    [classes,names]=tableVariableClasses(tbl);
    names=names(strcmp(classes,'char'));
    fieldNames=intersect(names,fieldNames);
    if 1
        for i=1:length(fieldNames)
            try
                tbl.(fieldNames{i})=regexprep(tbl.(fieldNames{i}),expression,repl);
            catch
                %tbl.(fieldNames{i})
                %fprintf('(field %s failed) ',fieldNames{i})
            end
        end
    else
        if ~isempty(fieldNames)
            tbl{:,fieldNames}=regexprep(tbl{:,fieldNames},expression,repl);
        end
    end
else
    for i=1:length(fieldNames)
        try
            %tbl{:,fieldNames{i}}=regexprep(tbl{:,fieldNames{i}},expression,repl);
            tbl.(fieldNames{i})=regexprep(tbl.(fieldNames{i}),expression,repl);
            fprintf('%s ',fieldNames{i})
        catch
            fprintf('(field %s failed) ',fieldNames{i})
        end
    end
end

%fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t0));
