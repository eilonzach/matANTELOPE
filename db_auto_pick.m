% This script will do some preliminaries for the tomography: 
% 1) go through events and stations and use TauP to put in predicted
% arrival times
% 2) remove picks that do not have a critical SNR in some window about the
% predicted arrival
% actually does this in one step: for each eq-sta pair, finds predicted P,
% S etc. times, calculate arrival time and SNR in the window around it, and
% put in a pick in the arrival table if the SNR meets a critical value
% Writes arrival and assoc tables. 

%% prelims
dbnam = 'pngbodydb';
dbdir = '/Volumes/Portadb/PORTAPNGBODY/';
phases = {'P','S','SKS','PP','SP'};
snrcrit = 2;
pretime = 20;
posttime = 20;
nZchan = 3; % number of vertical channel in "chans"

%% paranoia
if exist(sprintf('%s%s.arrival',dbdir,dbnam),'file'); system(sprintf('cp %s%s.arrival %s%s.arrival_old',dbdir,dbnam,dbdir,dbnam)); end
if exist(sprintf('%s%s.assoc',  dbdir,dbnam),'file'); system(sprintf('cp %s%s.assoc %s%s.assoc_old',    dbdir,dbnam,dbdir,dbnam)); end
try dbclose(db); end
%% start
db = dbopen(strcat(dbdir,dbnam),'r+');
dbwf = dblookup_table(db,'wfdisc');
dbor = dblookup_table(db,'origin');
dbsi = dblookup_table(db,'site');
dbarr = dblookup_table(db,'arrival');
dbass = dblookup_table(db,'assoc');
dbsich = dblookup_table(db,'sitechan');

norids = dbnrecs(dbor);
nstas = dbnrecs(dbsi);

for ie = 111:norids
    fprintf('DOING ORID %u/%u\n',ie,norids)
    dbor.record = ie-1; % ie_th event
    elat = dbgetv(dbor,'lat');
    elon = dbgetv(dbor,'lon');
    edep = dbgetv(dbor,'depth');
    orid = dbgetv(dbor,'orid');
    evtime = dbgetv(dbor,'time');

for ip = 1:length(phases)
    phase = char(phases(ip)); 
    fprintf('phase %s\n',phase)

for is = 1:nstas
    dbsi.record = is-1; % is_th sta
    slat = dbgetv(dbsi,'lat');
    slon = dbgetv(dbsi,'lon');
    sta = dbgetv(dbsi,'sta');
    
    tpt = tauptime('dep',edep,'evt',[elat elon],'sta',[slat slon],'p',phase);
    if isempty(tpt), continue; end
    phtime = evtime + tpt(1).time;
    
    % get data
    try
    [tt, dat, chans, nsamps, delta, seaz, depth,dbptrs,tref] = dbgetwf3(db, orid, phase, sta, pretime, posttime, 1);
    catch
        continue
    end
    samprate = round((unique(nsamps)-1)/(tt(end,1)-tt(1,1)));
    dat = detrend(dat(:,1:3));
    
    if length(unique(nsamps))>1, continue; end
    if any(nsamps < (pretime+posttime)*unique(samprate)), continue; end %skip if not enough data
    
    noind = 1 + [(pretime-10)*samprate:pretime*samprate]'; % noise: -10 to 0 seconds from arrival
    sigind = 1 + [(pretime)*samprate:(pretime+10)*samprate]'; % signal: 0 to 10 seconds from arrival
    snr = 0.5*max(abs(dat(sigind,:)))./std(dat(noind,:));% define as max_signal/2*std_noise
%     % plotting - can delete for real thing
%     figure(1)
%     plot(tt,dat)
%     title(sprintf('%s orid %u phase %s',sta,orid,phase))   
%     xlabel(sprintf('SNR = %.2f',snr))

    if max(snr) >= snrcrit
        if strcmp(phase(end),'P') % p-wave incident
            ic = nZchan;
        else % s-wave incident
            ic = find(snr==max(snr(setdiff(1:3,nZchan)))); % use horiz chan with max snr 
        end
        chan = char(chans(ic));
        %info for tables:
        arid = dbnrecs(dbarr)+1;
        dbsich.record = dbfind(dbsich,sprintf('sta == "%s" && chan == "%s"',sta,chan));
        chanid = dbgetv(dbsich,'chanid');
        seaz = azimuth(slat,slon,elat,elon,[6378.1,0.0033528]);
        esaz = azimuth(elat,elon,slat,slon,[6378.1,0.0033528]);
        
        % Arrival row
        dbarr = dblookup(db,'','arrival','','dbSCRATCH');
        dbputv(dbarr,'sta',sta,'time',phtime,'arid',arid,'chanid',1,'chan',chan,'iphase',phase,'azimuth',seaz,'snr',snr(ic));
        dbarr.record=dbadd(dbarr,'dbSCRATCH');
        % Assoc row
        dbass = dblookup(db,'','assoc','','dbSCRATCH');
        dbputv(dbass,'arid',arid,'orid',orid,'sta',sta,'phase',phase,'delta',tpt(1).distance,'seaz',seaz,'esaz',esaz);
        dbass.record=dbadd(dbass,'dbSCRATCH');
    end % if snr>crit
end % loop on stas
end % loop on phases
end % loop on events
    