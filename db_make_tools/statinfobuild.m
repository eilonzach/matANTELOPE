% Get station information from database and make a 'statinfo.txt' file for
% things like trxinbuild
dbnam='pngbodydb';
outfile='statinfo.txt';
odir='.';

%check for existence of file, and ask about deleting it
yn='n';
if exist(sprintf('%s/%s',odir,outfile),'file')==2
    yn=input('ofile already exists... overwrite? (y/n) ','s');
    if strcmp(yn,'y')==1
        delete(sprintf('%s/%s',odir,outfile))
    else
        return
    end
end

db=dbopen(dbnam, 'r');
dbsite=dblookup(db,'','site','','');
stas=dbgetv(dbsite,'sta');
dates=dbgetv(dbsite,'ondate');
slats=dbgetv(dbsite,'lat');
slons=dbgetv(dbsite,'lon');
selevs=dbgetv(dbsite,'elev');
cd(odir);

for is=1:length(stas)
fout = fopen(outfile,'a'); % or 'a' for append
% line of output
fprintf(fout,'%s\t 0.0\t %7.3f\t %8.3f\t %4.0f\t %4s\n',num2str(dates(is)),slats(is),slons(is),1000*selevs(is),char(stas(is)));
fclose(fout);
end