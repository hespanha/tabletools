function [classes,names]=tableVariableClasses(tbl)
% [classes,names]=tableVariableClasses(tbl)
%
% Returns the classes and names of the different variables in the table.
%
% Copyright (C) 2013-2014  Stacy & Joao Hespanha

%% problems since class is not consistent across rows
%% (e.g., one row may have empty as double and another a string)

classes=cellfun(@(x)class(x),table2cell(tbl(1,:)),'UniformOutput',0)';
names=tbl.Properties.VariableNames;
