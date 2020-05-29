function fix2xml(filename1,filename2);

fprintf('fix2xml: reading file... ');
t1=clock;
fid=fopen(filename1,'r');
file=fread(fid,inf,'char=>char')';
fclose(fid);
fprintf('done (%.3f sec)\n',etime(clock,t1));

fprintf('fix2xml: making replacements... ');
t1=clock;
file=regexprep(file,'<br[^>]*>','<br />');
file=regexprep(file,'<hr[^>]*>','<hr />');
fprintf('done (%.3f sec)\n',etime(clock,t1));

fprintf('fix2xml: writing file... ');
t1=clock;
fid=fopen(filename2,'w');
fwrite(fid,file);
fclose(fid);
fprintf('done (%.3f sec)\n',etime(clock,t1));

