function tbl=tableCheckClasses(tbl,noexception,nowarnings)
% tbl=tableVariableClasses(tbl)
%
% Throws an exception if any row does not match the types of the 1st
% row, unless all rows are either numeric or character strings, in
% which case it tries to convert all types to char.
%
% when the optional parameter noexception is true, tris to fix but
% does not through exception if it cannot.
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

if nargin<2
    noexception=false;
end

if nargin<2
    nowarnings=false;
end

names=tbl.Properties.VariableNames;
if size(tbl,1)==0
    classes=names;
    classes(:)={''};
    return
end

% type of first row
classes=cellfun(@(x)class(x),table2cell(tbl(1,:)),'UniformOutput',0)';

if ~nowarnings
    fprintf('  tableCheckClasses: ');
    t0=clock;
end
for i=1:length(classes)
    if ~nowarnings
        fprintf('%s(%s) %d/%d ',names{i},classes{i},i,length(classes));
    end
    if iscell(tbl.(names{i}))
        c=cellfun(@(x)class(x),tbl.(names{i}),'UniformOutput',0);
        k=find(~strcmp(c,classes{i}));
        if any(k)
            % easy fix?
            if all(ismember(c,{'char';'double';'datetime'}))
                if ~nowarnings
                    fprintf('\n    ATTENTION: type mismatch in variable %s(%s), converting all to char\n  ',names{i},classes{i});
                end
                kDate=find(cellfun(@(x)isdatetime(x),tbl.(names{i}),'UniformOutput',1));
                kNum=find(cellfun(@(x)isnumeric(x),tbl.(names{i}),'UniformOutput',1));
                kNan=kNum(isnan(cell2mat(tbl{kNum,names{i}})));
                kNum=setdiff(kNum,kNan);
                % replace datetime by string
                if ~isempty(kDate) && ~nowarnings
                    fprintf('    using datestr() in %d lines (%d..%d)\n',...
                            length(kDate),min(kDate),max(kDate));
                end
                tbl{kDate,names{i}}=cellfun(@(x)datestr(x),tbl{kDate,names{i}},'UniformOutput',0);
                % replace numbers by string
                if ~isempty(kNum) && ~nowarnings
                    fprintf('    using num2str() in %d lines (%d..%d)\n',...
                            length(kNum),min(kNum),max(kNum));
                end
                tbl{kNum,names{i}}=cellfun(@(x)num2str(x),tbl{kNum,names{i}},'UniformOutput',0);
                % replace NaN by empty string   
                if ~isempty(kNan) && ~nowarnings
                    fprintf('    replacing nan by '''' in %d lines (%d..%d)\n',length(kNan),min(kNan),max(kNan));
                end
                tbl{kNan,names{i}}=repmat({''},length(kNan),1);
                classes{i}='char';
                continue;
            end
            if ~noexception
                kk=mat2cell(k,ones(size(k)),1);
                [kk(1:min(10,end)),c(k(1:min(10,end)))]
                error('\ntype mismatch in variable %s(%s)\n',names{i},classes{i});
            end
        end
    elseif ~strcmp(class(tbl.(names{i})),classes{i})
        if ~noexception
            class(tbl.(names{i}))
            error('\ntype mismatch in variable %s(%s)\n',names{i},classes{i});
        end
    end
end
if ~nowarnings
    fprintf('done (%d rows, %d columns, %.3f sec)\n',size(tbl,1),size(tbl,2),etime(clock,t0));
end

end

