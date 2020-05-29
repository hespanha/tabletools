function k=tableSelect(tbl,query)
% k=tableSelect(tbl,query)
%
% Returns an array of table rows that match a given 'query'. 
% A query is a function whose input is the table and that returns the 
%
% E.g., 
% To find the rows whose length of the field name is smaller than 20, use:
%    k=tableSelect(tbl,@(tbl) find(length(tbl.name)<20) );
%
% Copyright (C) 2013-2014 Stacy & Joao Hespanha

%fprintf('tableSelect: %s ... ',char(query));
%t0=clock;
if isempty(query)
    k=(1:size(tbl,1))';
else
    k=find(query(tbl));
end
%fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t0));


