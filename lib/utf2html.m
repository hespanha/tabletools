function str=utf2html(str);
%
% Copyright (C) 2013-2014 Stacy & Joao Hespanha

k=find(str>255);
if isempty(k)
    return
end
htmlcodes=strcat('&+',dec2hex(double(str(k(:)))),';');
ut=cellstr(htmlcodes)';
%fprintf('length(str)=%d, found at positions ',length(str));disp(k)
str(k)=0;
str=strsplit(str,char(0),'CollapseDelimiters',false);
%fprintf('str split into %d pieces with lengths',length(str));disp(cellfun('length',str))
str=strjoin(str,ut);