% Script to make an Earthquake Catalogue file that SplitLab can use - i.e.
% in the format of the harvardCMT.mat file - from the origin (or other)
% table of an antelope db

dbnam = 'local4split';
table = 'origin';
ofile = 'localEQcat.mat';
region = 'SE PAPUA NEW GUINEA';
mag_default = 5;

%% ################## NO NEED TO EDIT AFTER THIS LINE ################## %%

db = dbopen(dbnam,'r');
dbt = dblookup_table(db,table);
norids = dbnrecs(dbt);
evtimes = dbgetv(dbt,'time');
lats = dbgetv(dbt,'lat');
lons = dbgetv(dbt,'lon');
depths = dbgetv(dbt,'depth');
try
    Mbs = dbgetv(dbt,'Mb');
    MSs = dbgetv(dbt,'MS');
catch me 
    Mbs = ones(norids,1)*mag_default;
    MSs = ones(norids,1)*mag_default;
end

cmt = struct([]);

    cmt(1).ID    = char(epoch2str(evtimes,'C%Y%m%d%H%MA'));
    cmt(1).year  = str2double(epoch2str(evtimes,'%Y'));
	cmt(1).month = str2double(epoch2str(evtimes,'%m'));
    cmt(1).day   = str2double(epoch2str(evtimes,'%d'));
    cmt(1).jjj   = str2double(epoch2str(evtimes,'%j'));
    cmt(1).hour  = str2double(epoch2str(evtimes,'%H'));
    cmt(1).minute = str2double(epoch2str(evtimes,'%M'));
    cmt(1).sec   = str2double(epoch2str(evtimes,'%S.%s'));
    cmt(1).lat   = lats;
    cmt(1).long  = lons;
    cmt(1).depth = depths;
    cmt(1).Mb = Mbs;
    cmt(1).MS = MSs;
    cmt(1).M0 = ones(norids,1);
    cmt(1).Mw = mag_default*ones(norids,1);
        for j = 1:norids
    cmt(1).region(j,:) = region;
        end
    cmt(1).strike = zeros(norids,1);
    cmt(1).dip = zeros(norids,1);
    cmt(1).rake = zeros(norids,1);

save(ofile,'cmt');
dbclose(db)
