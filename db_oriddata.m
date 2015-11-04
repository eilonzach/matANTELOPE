function [ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam )
% [ norids,orids,elats,elons,edeps,evtimes,mags ]  = db_oriddata( dbdir,dbnam )
%   Grab origin data from database

if nargin<2 % if only one arg given, assume it's dbnam, and we're in the right dir 
    dbdir = './';
    dbnam = dbdir;
end
if ~strcmp(dbdir(end),'/')
    dbdir = [dbdir,'/']; % append final slash if none
end

db = dbopen([dbdir,dbnam],'r');
dbor = dblookup_table(db,'origin');
[orids,elats,elons,edeps,evtimes,mags] = dbgetv(dbor,'orid','lat','lon','depth','time','ms');
norids = dbnrecs(dbor);
dbclose(db)

end

