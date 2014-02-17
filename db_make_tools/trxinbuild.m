% Builds a -.wfmeas file for the antelope "trexcerpt" function, so as to
% excerpt waveforms in explicit mode according to the rubric:
%     Use  the  times from table trxin.txt
%     to excerpt explicit waveform segments from  the  dbwf  database.
%     Each  row  in trxin.txt specifies a station, channel, and start time
%     for the excerpt.
%Inputs: a HarvardCMT.mat file, with cmt structure
%        a station file with rows e.g. (copy precisely; tabs/spaces matter)
%       turnontime      drift       latitude    longitude   depth   name
%     yyyy:jjj:hh:mm:ss	#.######	##.#### 	###.####	####	NAME 
% writes to outfile

%will be used with awk for the simpler trexcerpt convert mode 
%output will be relecant station, channel, starttime, duration

stafile='statinfo.txt';
outfile='trxin.txt';

%% Set criteria for accepted events
% NB THESE VALUES ARE INCLUSIVE!! ( <= or >= )
% MAGNITUDE
MSmin = 0;
Mwmin = 5.7;
% TIME
startdate = [2010,03,01]; %beginning of time window, in form [yyyy,mm,dd]
enddate = [2011,08,01]; %end of time window, in form [yyyy,mm,dd]
% DEPTH
mindepth = 0; %minimum depth (km)
maxdepth = 1000; %maximum depth (km)
% DISTANCE 
refpt=[-9.8,150.5]; % [lat,long] of reference point for distance constraint
mindeg = 30; %minimum angular distance (0?degrees?180)
maxdeg = 140; %maximum angular distance (0?degrees?180)
% BACKAZ
minaz=0;
maxaz=360;

%Windowing options
%1) Use phases - {'phase_1','phase_2'} NB uses first occurence by default
phases={'P','PS'}; 
runup = 100; % secs before predicted phase_1 arrival to start window
rundown = 100; % secs after predicted phase_2 arrival to end window
%2) Use abs. times after evtime. DEFAULTS to this if phases fails.
windstart = 12*60; % time after evtime to start window
windend = 32*60; % time after evtime to end window

chans={'BH0','BHE','BH1','BHN','BHZ'}; % e.g. {'BH0','BHE','BH1','BHN','BHZ'};

%% ################## NO EDITS NEEDED BELOW HERE ################## %%

%% Get event info
try
[goodev]=gdevts_fn(MSmin,Mwmin,startdate,enddate,mindepth,maxdepth,refpt,mindeg,maxdeg,minaz,maxaz);
catch me
    error('No events fit criteria')
end

%% Get station information
fid=fopen(stafile,'r');
% [time,drift,lat,long,depth,#name]...
D=textscan(fid,'%s %f %f %f %f %s','delimiter','\t');
fclose(fid);
stas=D{6};
infostats={'lat','long','depth'}; infocolumns=[3,4,5];
for j=1:length(infostats)
    assignin('base',sprintf('sta%s',char(infostats(j))),D{infocolumns(j)});
end

% NB need to get land station information - ultimately need to get out of
% the .txt environment and start using the antelope-matlab
% interface

fout = fopen(outfile,'w'); % or 'a' for append
% Header
fprintf(fout, '%s \t%s  %s  %s  \t%s %s\n','sta','chan','yr-jday  hr:mi:sec','windl','ph1','ph2');
fclose(fout);

nos=length(stas);
noe=length(goodev.evtime);
noc=length(chans);
fprintf('%u events and %u stations being written to file with %u chan(s)\n',noe,nos,noc);
for is=1:nos
for ie=1:noe

%% 1)Window by absolute time w.r.t. evtime
% this is the default, and the window is defined above
tstart = windstart;
tend = windend;
%% 2)Window by certain phase arrival...
if length(phases{1}) > 0
    phs = phases;
    try % some phases may not work because wrong distances
    for ip = 1:2;
    tt=taupTime([],goodev.depth(ie),phases{ip},...
        'sta',[stalat(is) stalong(is)],...
        'evt',[goodev.lat(ie) goodev.long(ie)]);
    assignin('base',sprintf('ph%ut',ip),tt(1).time);
    end    
    tstart = ph1t-runup; % s after evtime to start
    tend = ph2t+rundown; % s after evtime to end
    catch me
        phs = {'',''};
        tstart = windstart;
        tend = windend;
    end
end

%% Window length
    winddur=tend-tstart;
    
%% parse starttime into yyy-jjj hh:mm:ss
    ws=rem(tstart,60);
    wmi=floor(tstart/60);
    
    SEC=rem(goodev.sec(ie)+ws,60); MIpl=floor((goodev.sec(ie)+ws)/60);
    MI=rem(goodev.minute(ie)+wmi+MIpl,60); HRpl=floor((goodev.minute(ie)+wmi+MIpl)/60);
    HR=rem(goodev.hour(ie)+HRpl,24); DAYpl=floor((goodev.hour(ie)+HRpl)/24);
    DAY=goodev.day(ie)+DAYpl;
%     Jday=(datenum(eval(num2str(evyear(ie))),eval(num2str(evmonth(ie))),eval(num2str(DAY)))-datenum(eval(num2str(evyear(ie))),0,0));
    Jday=goodev.jjj(ie)+DAYpl;
    
    % parse window length into hh:mm:ss
    wls=ceil(rem(winddur,60));
    wlmi=floor(rem(winddur,3600)/60);
    wlh=floor(winddur/3600);
    

    %% annoying things to get all time strings into right formats: 
    %annoying things to get starttime d.o.y. in right format (3-digit string)
    if Jday >= 100
        assignin('base','Jday',sprintf('%s',num2str(Jday)));
    elseif Jday >= 10 && Jday < 100
        assignin('base','Jday',sprintf('0%s',num2str(Jday)));
    elseif Jday < 10
        assignin('base','Jday',sprintf('00%s',num2str(Jday)));
    end
    %annoying things to get starttime HR in right format (2-digit string)
    if HR >= 10
        assignin('base','HR',sprintf('%s',num2str(HR)));
    elseif HR < 10
        assignin('base','HR',sprintf('0%s',num2str(HR)));
    end
    %annoying things to get starttime MI in right format (2-digit string)
    if MI >= 10
        assignin('base','MI',sprintf('%s',num2str(MI)));
    elseif MI < 10
        assignin('base','MI',sprintf('0%s',num2str(MI)));
    end
    %annoying things to get starttime SEC in right format (2-digit string)
    SEC=round(SEC);
    if SEC >= 10
        assignin('base','SEC',sprintf('%s',num2str(SEC)));
    elseif SEC < 10
        assignin('base','SEC',sprintf('0%s',num2str(SEC)));
    end
    %annoying things to get windl wlh in right format (2-digit string)
    if wlh >= 10
        assignin('base','wlh',sprintf('%s',num2str(wlh)));
    elseif wlh < 10
        assignin('base','wlh',sprintf('0%s',num2str(wlh)));
    end
    %annoying things to get windl wlmi in right format (2-digit string)
    if wlmi >= 10
        assignin('base','wlmi',sprintf('%s',num2str(wlmi)));
    elseif wlmi < 10
        assignin('base','wlmi',sprintf('0%s',num2str(wlmi)));
    end
    %annoying things to get windl wls in right format (2-digit string)
    if wls >= 10
        assignin('base','wls',sprintf('%s',num2str(wls)));
    elseif wls < 10
        assignin('base','wls',sprintf('0%s',num2str(wls)));
    end

    
 %% write to outfile
for ic=1:noc
fout = fopen(outfile,'a'); % or 'a' for append
% line of output
fprintf(fout,'%s  \t%s   %s\t %s\t %s\t %s\n',char(stas(is)),char(chans(ic)),sprintf('%s-%s %s:%s:%s',num2str(goodev.year(ie)),Jday,HR,MI,SEC),...
    sprintf('%s:%s:%s',wlh,wlmi,wls),...
    char(phs(1)),char(phs(2)));
fclose(fout);
end %loop on chans 
end %loop on evs
end %loop on stas

fprintf('%u events and %u stations were written to file with %u chan(s)\n',noe,nos,noc);




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