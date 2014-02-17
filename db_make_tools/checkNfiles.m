% Script to compare wfdisc rows to data on my disc and on poas
clear all
dbnam = '/Volumes/Portadb/PORTAPNGBODY/pngbodydb';
homedata = '/Volumes/Portadb/PORTAPNGBODY/pngbodydb_tables/data/';
poasdata = '/Volumes/zeilon/CdpapuaBodyWaveDb/pngbodydb_tables/data/';

db = dbopen(dbnam,'r');
dbwf = dblookup_table(db,'wfdisc');
dbsite = dblookup_table(db,'site');
stas = dbgetv(dbsite,'sta');
nstas = length(stas);
for is = 1:nstas
    dbs = dbsubset(dbwf,sprintf('sta =="%s"',char(stas(is))));
    nrecs = dbnrecs(dbs);
    [~,nhome] = system(sprintf('ls -l %s%s | wc -l',homedata,char(stas(is))));
    [~,naway] = system(sprintf('ls -l %s%s | wc -l',poasdata,char(stas(is))));
    fprintf('sta %s: %.f recs - %.f here and %.f on poas\n',char(stas(is)),nrecs,eval(nhome),eval(naway));
end
dbclose(db)