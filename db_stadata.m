function [ nstas,stas,slats,slons,selevs,ondate,offdate,staname ] = db_stadata( dbdir,dbnam )
% [ nstas,stas,slats,slons,selevs,ondate,offdate,staname ] = db_stadata( dbdir,dbnam )
%   Grab station data from database

if nargin<2 % if only one arg given, assume it's dbnam, and we're in the right dir 
    dbdir = './';
    dbnam = dbdir;
end
if ~strcmp(dbdir(end),'/')
    dbdir = [dbdir,'/']; % append final slash if none
end

db = dbopen([dbdir,dbnam],'r');
dbsi = dblookup_table(db,'site');
[stas,slats,slons,selevs,ondate,offdate,staname] = dbgetv(dbsi,'sta','lat','lon','elev','ondate','offdate','staname');
nstas = dbnrecs(dbsi);
dbclose(db)

end

