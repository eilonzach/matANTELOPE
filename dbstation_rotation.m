%OBS_DO_ROTATION
%  Do all the rotations for all the OBS data to get into proper Z-N-E form
%  Input:   needs portapng db  with .wfdisc .site .origin  table
%           needs origin table to only have events which satisfy correct
%           mag,distance conditions
%  Output:  portapng_rot.wfdisc & data in data/rotSTA directory 
%
%	Z. Eilon 23 May 2013

%   Station names and corrections 
stas= {'B'};   
%OBS:   {'B','D', 'E', 'F', 'G', 'H', 'J'};
%Land:  {'PEMM', 'KEIA', 'JONE', 'GOGO'};
%% NB correction angles are the angles from true north that the station norths are pointing.
scors= [116.8];
%OBS:   [ 116.8  0.3  183.1  301  47.2  324.5  143.9]; 
%Land:  [15.6  12.1  21.4  17.2]; 

dbnam='pngbodydb';
dbout='pngbodydb_rot';
tdir = '/Volumes/Portadb/PORTAPNGBODY/pngbodydb_tables/';

ichans=char({'BH0','BH1','BHZ'}); % input channels in order e ,n ,z
ochans=char({'BHE','BHN','BHZ'}); % output channels in order e ,n ,z
fchans=char({'e','n','z'});       % channel suffices for sac files

ynall = 'yall'; %permission: yall/n
noplot = 0;

% Processing starts here
wd = pwd;
cd(tdir); % need to be in the dir w/ wfdisc table, assuming all dir fields are w. reference to this
try(dbclose(db)); end % make sure.
%% loop over stations
for is=1:length(stas)
sta=char(stas(is));
odir=sprintf('data/rot%s/',sta); % don't forget final "/"!
xdir=sprintf('data/pre_rot%s/',sta); % don't forget final "/"!
if exist(odir,'dir')==0, mkdir(odir); end
if exist(xdir,'dir')==0, mkdir(xdir); end

scor=scors(is);
sscor=sind(scor);
cscor=cosd(scor);
    
db=dbopen(dbnam,'r+');
dbwf = dblookup_table(db,'wfdisc');
dbs1 = dbsubset(dbwf,sprintf('sta =="%s"',sta));
nrecs = dbnrecs(dbs1);
wfids = dbgetv(dbs1,'wfid');
%% Sitechan edits + new rows if needed
dbsch = dblookup_table(db,'sitechan');
chanids = zeros(3,1);
%EAST
if ~dbnrecs(dbsubset(dbsch,sprintf('sta == "%s" && chan == "%s"',sta,ochans(1,:)))) % if ochan doesn't exist yet, make new row
chanids(1) = max(dbgetv(dbsch,'chanid'))+2;
dbsch.record = dbfind(dbsch,sprintf('sta == "%s" && chan == "%s"',sta,ichans(1,:)));
record = dbget(dbsch);
dbsch.record = dbaddnull(dbsch);
dbput(dbsch,record);
dbputv(dbsch,'hang',90.0,'vang',90.0,'chan',ochans(1,:),'chanid',chanids(1))
else % chan exists, make sure hang etc. ok
dbsch.record = dbfind(dbsch,sprintf('sta == "%s" && chan == "%s"',sta,ichans(1,:)));
chanids(1) = dbgetv(dbsch,'chanid');
dbputv(dbsch,'hang',90.0,'vang',90.0)
end
%NORTH
if ~dbnrecs(dbsubset(dbsch,sprintf('sta == "%s" && chan == "%s"',sta,ochans(2,:)))) % if ochan doesn't exist yet, make new row
chanids(2) = max(dbgetv(dbsch,'chanid'))+1;
dbsch.record = dbfind(dbsch,sprintf('sta == "%s" && chan == "%s"',sta,ichans(2,:)));
record = dbget(dbsch);
dbsch.record = dbaddnull(dbsch);
dbput(dbsch,record);
dbputv(dbsch,'hang',0.0,'vang',90.0,'chan',ochans(2,:),'chanid',chanids(2))
else % chan exists, make sure hang etc. ok
dbsch.record = dbfind(dbsch,sprintf('sta == "%s" && chan == "%s"',sta,ichans(2,:)));
chanids(2) = dbgetv(dbsch,'chanid');
dbputv(dbsch,'hang',0.0,'vang',90.0)
end
%VERTICAL
dbsch.record = dbfind(dbsch,sprintf('sta == "%s" && chan == "%s"',sta,ichans(3,:)));
chanids(3) = dbgetv(dbsch,'chanid');

dbclose(db)
usedwfids = [];

%% loop through records
% find all 3 components for each time, and cycle up in increasing wfid. 
for irec = 1:nrecs
wfid = wfids(irec);
if any(usedwfids == wfid)==1, continue; end % move on if already done this wfid
fprintf('Rotating sta %s, record %u/%u\n',sta,length(usedwfids)+1,nrecs);
%% find all other records starting coincident with this one
db=dbopen(dbnam,'r+');
dbwf = dblookup_table(db,'wfdisc');
dbs1 = dbsubset(dbwf,sprintf('sta =="%s"',sta));
dbs1.record = dbfind(dbs1,sprintf('wfid == %u',wfid));
rtime = dbgetv(dbs1,'time'); % get time of first wf
dbs1.record=-501; %reset to whole table

dbs2 = dbsubset(dbs1,sprintf('time - %f < 0 && time - %f > 0 ',rtime+0.5,rtime-0.5)); %subset to same record start time
if dbnrecs(dbs2)>3; fprintf('ERROR... DUPLICATE RECORDS?'); end
t0s = dbgetv(dbs2,'time');
t1s = dbgetv(dbs2,'endtime');
t0 = max(t0s); %starttime that will include all records
t1 = min(t1s); %endtime that will include all records

%% Check all chans are there and take relevant action
dbs2e = dbsubset(dbs2,sprintf('chan == "%s"',ichans(1,:))); dbs2e.record=0; %use only 1st
dbs2n = dbsubset(dbs2,sprintf('chan == "%s"',ichans(2,:))); dbs2n.record=0; %use only 1st
dbs2z = dbsubset(dbs2,sprintf('chan == "%s"',ichans(3,:))); dbs2z.record=0; %use only 1st
usedwfids = cat(1,usedwfids,dbgetv(dbs2,'wfid'));

yn = ynall;
if dbnrecs(dbs2e)~=1 && dbnrecs(dbs2n)==1 % no E yes N
    fprintf('No E data\n')
    % delete old file query?
        if strcmp(yn,'yall')==0
            yn=input('Do you want to delete/move file? (y/n/yall) ','s');
        end
        if strcmp(yn,'y')==1 || strcmp(yn,'yall')==1 
            movefile(strcat(dbgetv(dbs2n,'dir'),'/',dbgetv(dbs2n,'dfile')),...
                     strcat(xdir,dbgetv(dbs2n,'dfile')))
%             dbwf.record = dbfind(dbwf,sprintf('sta == "%s" && chan == "%s"',sta,ichans(2,:)));
%             dbdelete(dbwf);  dbwf = dblookup_table(db,'wfdisc'); %delete N from wfdisc        
        end

elseif dbnrecs(dbs2n)~=1 && dbnrecs(dbs2e)==1 % no N yes E
    fprintf('No N data\n')
    % delete old file query?
        if strcmp(yn,'yall')==0
            yn=input('Do you want to delete/move file? (y/n/yall) ','s');
        end
        if strcmp(yn,'y')==1 || strcmp(yn,'yall')==1 
            movefile(strcat(dbgetv(dbs2e,'dir'),'/',dbgetv(dbs2e,'dfile')),...
                     strcat(xdir,dbgetv(dbs2e,'dfile')))
%             dbwf.record = dbfind(dbwf,sprintf('sta == "%s" && chan == "%s"',sta,ichans(1,:)));
%             dbdelete(dbwf);  dbwf = dblookup_table(db,'wfdisc'); %delete E from wfdisc
        end

elseif dbnrecs(dbs2e)==1 && dbnrecs(dbs2n)==1 % have both horiz - get data!

%% Old East
[tte, dat_e, x, nsamps, samprate, wfide] = dbgetwfz(db, sta, t0, t1, 'epoch', ichans(1,:));
[calib(1,:) instype(1,:) segtype(1,:)] = dbgetv(dbs2e,'calib','instype','segtype');
movefile(strcat(dbgetv(dbs2e,'dir'),'/',dbgetv(dbs2e,'dfile')),...
         strcat(xdir,dbgetv(dbs2e,'dfile')))
% dbwf.record = dbfind(dbwf,sprintf('sta == "%s" && chan == "%s"',sta,ichans(1,:)));
% dbdelete(dbwf); dbwf = dblookup_table(db,'wfdisc'); %delete E from wfdisc
%% Old North
[ttn, dat_n, x, nsamps, samprate, wfidn] = dbgetwfz(db, sta, t0, t1, 'epoch', ichans(2,:));
[calib(2,:) instype(2,:) segtype(2,:)] = dbgetv(dbs2n,'calib','instype','segtype');
movefile(strcat(dbgetv(dbs2n,'dir'),'/',dbgetv(dbs2n,'dfile')),...
         strcat(xdir,dbgetv(dbs2n,'dfile')))
% dbwf.record = dbfind(dbwf,sprintf('sta == "%s" && chan == "%s"',sta,ichans(2,:)));
% dbdelete(dbwf); dbwf = dblookup_table(db,'wfdisc'); %delete N from wfdisc
%% Old Vertical - if it exists
if dbnrecs(dbs2z)
[ttz, dat_z, x, nsamps, samprate, wfidz] = dbgetwfz(db, sta, t0, t1, 'epoch', ichans(3,:));
[calib(3,:) instype(3,:) segtype(3,:)] = dbgetv(dbs2z,'calib','instype','segtype');
movefile(strcat(dbgetv(dbs2z,'dir'),'/',dbgetv(dbs2z,'dfile')),...
         strcat(xdir,dbgetv(dbs2z,'dfile')))
% dbwf.record = dbfind(dbwf,sprintf('sta == "%s" && chan == "%s"',sta,ichans(3,:)));
% dbdelete(dbwf);  dbwf = dblookup_table(db,'wfdisc'); %delete Z from wfdisc
else
    dat_z = zeros(size(dat_n)); ttz = ttn; wfidz = [];
    calib(3,:)=calib(2,:); instypye(3,:)=instype(2,:); segtype(3,:)=segtype(2,:);
end
dbclose(db)
% order the wfids for putting back in later
owfids = [wfide;wfidn;wfidz];

% ROTATE 
        tt = linspace(t0,t1-1./samprate,nsamps)';
%       linear interpolate to same times
        dat_n = interp1(ttn,dat_n,tt);
        dat_e = interp1(tte,dat_e,tt);
        dat_z = interp1(ttz,dat_z,tt); 
%       rotate
        data_N =  dat_n.*cscor - dat_e.*sscor;
        data_E =  dat_n.*sscor + dat_e.*cscor;
        dat=[data_E,data_N,dat_z];
        
if ~noplot
% view pre-rotated
    figure(13);
    subplot(3,2,1); plot(tte,dat_e,'b'); title(sprintf('Component %s pre',ichans(1,:)));
    subplot(3,2,3); plot(ttn,dat_n,'b'); title(sprintf('Component %s pre',ichans(2,:)));
    subplot(3,2,5); plot(ttz,dat_z,'b'); title(sprintf('Component %s pre',ichans(3,:)));
% view post-rotated
    subplot(3,2,2); plot(tt,dat(:,1),'r'); title(sprintf('Component %s post',ochans(1,:)));
    subplot(3,2,4); plot(tt,dat(:,2),'r'); title(sprintf('Component %s post',ochans(2,:)));
    subplot(3,2,6); plot(tt,dat(:,3),'r'); title(sprintf('Component %s post',ochans(3,:)));
end

%% Put new chans into wfdisc
yn = ynall;

% Save the trace data in a new database, with the underlying file in sac format:
dbo=dbopen(dbout,'r+');
dbow=dblookup_table(dbo,'wfdisc');
for ic=1:3 % E is 1, N is 2, Z is 3
% Set up trace and its parameters
tr=trnew;
tr=dblookup_table(tr,'trace');
wfchan=ochans(ic,:);
dfile=sprintf('%s.%s.sac.%s',epoch2str(t0,'%Y.%j.%H.%M.%S'),sta,fchans(ic));

% overwrite out file query?
ofile=dir(sprintf('%s/%s',odir,dfile));
if isempty(ofile)==0 
    if strcmp(yn,'yall')==0
    yn=input('Do you want to overwrite file? (y/n/yall) ','s');
    end
    if strcmp(yn,'y')==1 || strcmp(yn,'yall')==1 
        delete(sprintf('%s/%s',odir,dfile));
    elseif strcmp(yn,'n')==1
        continue
    end
end

% Put the waveform into the trace-object:
tr.record=dbaddv(tr,'net', 'ROT','sta',sta,'chan',wfchan,...
'nsamp',nsamps,'samprate',samprate,'time',tt(1),'endtime', tt(end));%,'wfid',owfids(ic));
trinsert_data(tr,dat(:,ic));
trsave_wf(tr,dbow,'sc',sprintf('%s/%s',odir,dfile));
trdestroy( tr ) % Clean up
% rename the wfids!
unid = sprintf('sta=="%s" && chan=="%s" && time==%f',sta,wfchan,tt(1)); % unique id of this trace
dbow.record = dbfind(dbow,unid);
dbputv(dbow,'wfid',owfids(ic),'chanid',chanids(ic),'calib',calib(ic,:),'instype',instype(ic,:),'segtype',segtype(ic,:)); %change to apt wfid
dbow.record = -501; % reset to all
end % loop on channels
dbclose(dbo)


end %have both horiz; rotate
% pause
end % loop on records
end % loop on stations

try dbclose(db); end

% for is=1:length(stas)
% sta=char(stas(is));
% system(sprintf('dbsubset -v %s.wfdisc "sta == ''%s''" | dbdelete -v -',dbnam,sta))
% end
% system(sprintf('dbmerge %s %s',dbout, dbnam))
cd(wd)