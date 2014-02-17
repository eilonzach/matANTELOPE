% Uses an existing database (will make use of site, sitechan, and origin
% tables) to make a trxin.txt file for use with awk and trexcerpt to extract
% data snippets from a master database
% 
% output will be relevant station, channel, start time, duration
%
% Z. Eilon      1 July 2013

dbnam = 'pngbodydb';
dbdir = '/Volumes/Portadb/PORTAPNGBODY/'; % don't forget final "/"
ofile = 'trxin_plus.txt';
odir = '/Volumes/Portadb/PORTAPNGBODY/dbbuilding_files/'; % don't forget final "/"

%% Subsetting database
orids = [18,23,38,49,56,68,74,78,80,90:95,102,109,113,115,117,121,138,...
         143,159,164,169,175,180,181,185,190,204,216,233,238,245,247,253,...
         254,257,270,273,275,285,287,322,327,328,334,335]';     % e.g. [1,19,22:27,89]';
% if empty or zero then does all apart from...
orids_skip = []';

stas = {}';         % e.g. {'KIR','BASI','AGAN'}';
% if empty then does all apart from...
stas_skip = {'I'}';    % e.g. {'I'};

chans={}';          % e.g. {'BH0','BHE','BH1','BHN','BHZ'}';
% if empty then does all it can find for each station apart from...
chans_skip = {'LOG','BDH'}';   % e.g. {'LOG'}'

%Windowing options
%1) Use phases = {'phase_1','phase_2'} NB uses first occurence by default
phases={}; % can give just one phase to window around just one arrival
runup = 100; % secs before predicted phase_1 arrival to start window
rundown = 100; % secs after predicted phase_2 arrival to end window
%2) Use abs. times after evtime. DEFAULTS to this if phases is empty/fails.
windstart = 0*60; % time in seconds after evtime to start window
windend = 30*60; % time in seconds after evtime to end window

%% ################## NO EDITS NEEDED BELOW HERE ################## %%
%% prep file and write header
fout = fopen(strcat(odir,ofile),'w'); % or 'a' for append
fprintf(fout, '%-4s   %-4s   %-21s   %-12s   %-6s   %-6s\n','sta','chan','yr-jday  hr:mi:sec','windl','ph1','ph2');
fclose(fout);

if isempty(phases), phs = {''}; else phs = phases; end

% Open db and get going
try dbclose(db); end
db = dbopen(strcat(dbdir,dbnam),'r');

dbor = dblookup_table(db,'origin');
if isempty(orids)
    orids = dbgetv(dbor,'orid');
    orids = setdiff(orids,orids_skip);
end

dbsi = dblookup_table(db,'site');
if isempty(stas)
    stas = dbgetv(dbsi,'sta');
    stas = setdiff(stas,stas_skip);
end

dbsich = dblookup_table(db,'sitechan');

%% Orid info and loop
for ie = 1:length(orids)
orid = orids(ie);
dbor.record = dbfind(dbor,sprintf('orid == %u',orid));
evdeet = db2struct(dbor);

%% Station and chan info and loop
for is = 1:length(stas)
sta = char(stas(is));
dbsi.record = dbfind(dbsi,sprintf('sta == "%s"',sta));
stdeet = db2struct(dbsi);
    if isempty(chans)
    dbsich1 = dbsubset(dbsich,sprintf('sta == "%s"',sta));
    chs = dbgetv(dbsich1,'chan');
    chs = setdiff(chs,chans_skip);
    else 
    chs = chans;
    end

%% 2)Window by absolute time w.r.t. evtime
% this is the default, and the window is defined above
tstart = windstart + evdeet.time;
tend = windend + evdeet.time;
%% 1)Window by phase arrival
if ~isempty(phases)
    try % some phases may not work because wrong distances
        tt = tauptime('ph',phases{1},...
                      'dep',evdeet.depth,...
                      'sta',[stdeet.lat stdeet.lon],...
                      'evt',[evdeet.lat evdeet.lon]);
        tstart = tt(1).time - runup + evdeet.time;% s after evtime to start
        tt = tauptime('ph',phases{end},...
                      'dep',evdeet.depth,...
                      'sta',[stdeet.lat stdeet.lon],...
                      'evt',[evdeet.lat evdeet.lon]);
        tend = tt(1).time + rundown + evdeet.time; % s after evtime to end
    catch me
        tstart = windstart + evdeet.time;
        tend = windend + evdeet.time;
    end
end


%% Window 
    winddur=tend-tstart;
    timstring = epoch2str(tstart,'%Y-%j %H:%M:%S.%s');
    durstring = epoch2str(winddur,'%H:%M:%S.%s');

 %% write to outfile
fout = fopen(strcat(odir,ofile),'a'); % or 'a' for append
for ic=1:length(chs)
    chan = char(chs(ic));
    
% line of output
fprintf(fout,'%-4s   %-4s   %21s   %12s   %-6s   %-6s\n',sta,chan,timstring,durstring,char(phs(1)),char(phs(end)));

end %loop on chans 
fclose(fout);
end %loop on stas
end %loop on orids
dbclose(db)
fprintf('%u events and %u stations were written to file\n',length(orids),length(stas))



%% Bear in mind the trexconv4splt build script
% #!/bin/bash
% # Updated 11/1/2012 10:12
% # This script requires a file (trxin.txt) of stations, channels, and start times for the waveform window 
% ## statA chanA timeA durA
% ## statB chanB timeB durB
% ## statC chanC timeC durC
% ## etc. 
% 
% # stat* is the name of the station (must match the antelope tables)
% # chan* is the name of the channel (must match the antelope tables)
% # time* is the time for the beginning of the window around the relevant arrival,
% #       in the format YYYY-DDD HH:MM:SS
% # The trexcerpt command builds sac files from the window accordingly
% 
% ## Required inputs - the name of the database, the channel names,
% ## and the file with the correction angles
% db=/Volumes/zeilon/PNG/cdpapuall	  #for example
% chanE=BH0	#for example
% chanN=BH1	#for example
% chanZ=BHZ	#for example
% echo "Give the name of the file containing the time excerpt data"
% read INFILE
% i=1
% ## go through line by line getting information
% while [ $i -le  $(awk 'END { print NR }' < $INFILE) ]; do
% stat=$(awk -v ii=$i 'NR==ii {print $1}'< $INFILE)
% chan=$(awk -v ii=$i 'NR==ii {print $2}'< $INFILE)
% timestart=$(awk -v ii=$i 'NR==ii {printf("%s %s",$3,$4)}'< $INFILE)
% windl=$(awk -v ii=$i 'NR==ii {print $5}'< $INFILE)
% 
% if [ "$chan" = "$chanE" ]; then
% CHAN=e
% elif [ "$chan" = "$chanN" ]; then
% CHAN=n
% elif [ "$chan" = "$chanZ" ]; then
% CHAN=z
% fi
% 
% # do the excerpting
% trexcerpt -vvo sc -w EXCERPT/%{sta}/%Y.%j.%H.%M.%S.%{sta}.sac.$CHAN\
%  -c "chan=='$chan' && sta=='$stat'" $db trialwf "$timestart" $windl
% # N.B. If a verbose output is not required, remove the "-v" in the two lines above
% # N.B. If confirmations for each change are not required, remove the "-c" in the two lines above
% let i=i+1
% done