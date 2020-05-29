function cellstr=tableMergeColumns(tbl,fieldNames,separator)
% cellstr=tableMergeFields(tbl,fieldNames,separator)
%
% Returns a string array with all the table fields in the cell array
% fieldNames merges.  All the fields in fieldNames should be strings.
% 
% The (optional) parameter separator is a string that should be
% inserted between the fields to be merged.
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha


if nargin<3
    separator=[];
end

if ~iscell(fieldNames)
    fieldNames={fieldNames};
end

fprintf('tableMergeFields... ');
t1=clock;

x=tbl{:,fieldNames};

cellstr=cell(size(x,1),1);
if size(x,2)==1
    for row=1:size(x,1)
        cellstr{row}=x{row,1};
    end
else
    for row=1:size(x,1)
        y(1:2:2*size(x,2)-1)=x(row,:);
        y(2:2:2*size(x,2)-1)={separator};
        cellstr{row}=[y{:}];
        
    end
end

fprintf('done (%.3f sec)\n',etime(clock,t1));

