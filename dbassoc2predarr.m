% This script uses assoc and arrival tables (as well as site and orid)
% to make the predarr table with predicted arrival times using taup, and
% specifyint event and station locations. 

%% Starts
dbnam = 'pngbodydb';
dbdir = '/Users/Zach/Documents/MATLAB/PNG_tomog/PNGbodydb/';% remember final "/"
dbtabdir = '/Users/Zach/Documents/MATLAB/PNG_tomog/PNGbodydb/'; % remember final "/"

% see if predarr table already exists, if so move it to predarr_old
if exist(strcat(dbtabdir,dbnam,'.predarr'),'file') == 2
    yn = 'n'; %by default
    yn = input('Predarr table seems to exist. Do you want to overwrite? (y/n) ','s');
    if strcmp(yn,'y')==1,   movefile(strcat(dbtabdir,dbnam,'.predarr'),strcat(dbtabdir,dbnam,'.predarr_old')); 
    else return; end
end

db = dbopen(strcat(dbdir,dbnam),'r+');
dbsi = dblookup_table(db,'site');
dbor = dblookup_table(db,'origin');
dbarr = dblookup_table(db,'arrival');
dbass = dblookup_table(db,'assoc');

nass = dbnrecs(dbass);
arids = dbgetv(dbass,'arid');

dbpredarr = dblookup_table(db,'predarr');
dbrecord = dbaddnull(dbpredarr);
% now loop through arrivals putting in predictio0ns
for ias = 1:nass
    if any(ias==0:round(nass/20):nass), fprintf('%.0f%% done\n',100*ias/nass); end
    dbsass = dbsubset(dbass,sprintf('arid == %u',arids(ias)));
    dbjars = dbjoin(dbsass,dbarr);
    dbjarsor = dbjoin(dbjars,dbor);
    dbjarsorsi = dbjoin(dbjarsor,dbsi);
    dbj = dbjarsorsi;
    if dbnrecs(dbj)~=1; return; end
    [orid,elat,elon,edep,evtime,slat,slon,phase,seaz,esaz,delta] =...
     dbgetv(dbj,'orid','origin.lat','origin.lon','origin.depth','origin.time',...
                'site.lat','site.lon','iphase','seaz','esaz','delta');
%     ttaup = tauptime('p',phase,'z',edep,'s',[slat slon],'e',[elat elon]);
    ttaup = taupant('p',phase,'z',edep,'s',[slat slon],'e',[elat elon]);
    if isempty(ttaup)
        if strcmp(phase,'P')
            try
%             ttaup = tauptime('p','Pdiff','z',edep,'s',[slat slon],'e',[elat elon]);
            ttaup = taupant('p','Pdiff','z',edep,'s',[slat slon],'e',[elat elon]);
            catch, fprintf('No predar for arid %u\n',arids(ias)); continue
            end
        elseif strcmp(phase,'S')
            try
%             ttaup = tauptime('p','Sdiff','z',edep,'s',[slat slon],'e',[elat elon]);
            ttaup = taupant('p','Sdiff','z',edep,'s',[slat slon],'e',[elat elon]);
            catch, fprintf('No predar for arid %u\n',arids(ias)); continue
            end
%         elseif strcmp(char(phases(ip)),'PKP')
%             try
%             ttaup = tauptime('p','PKiKP','z',edep,'s',[slat slon],'e',[elat elon]);
%             catch, fprintf('No predar for arid %u\n',arids(ias)); continue
%             end
        else
            fprintf('No predar for arid %u\n',arids(ias)); continue
        end
    end
    % put new row into assoc table.
    db.record=dbaddv(dbpredarr,'arid',arids(ias),'orid',orid,'time',evtime+ttaup(1).time,...
    'slow',ttaup(1).rayparameter,'seaz',seaz,'esaz',esaz);
end % loop on arrivals
dbcrunch(dbpredarr);

dbclose(db)

% movefile(strcat(dbdir,dbnam,'.predarr'),dbtabdir);