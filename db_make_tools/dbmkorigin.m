% This script takes an existing origin table with ?1 rows and adds all
% "good" events from a HarvardCMT.mat file (which has cmt structure) that
% are not already in the origin database. One must specify a time and
% magnitude range for "good" events.
clear('Mmin','startdate','enddate','mindepth','maxdepth','refpt','mindeg','maxdeg');
db='pngbodydb_or';

%% Set criteria for accepted events
% MAGNITUDE
% Mmin=6.5; 
Mwmin=5.7;
MSmin=0;
Mbmin=0;
M0min=0;
% TIME
startdate=[2010,03,01]; %beginning of time window, in form [yyy,mm,dd]
enddate=[2011,08,01]; %end of time window, in form [yyy,mm,dd]
% DEPTH
mindepth=0; %minimum depth (km)
maxdepth=1000; %maximum depth (km)
% DISTANCE - comment out to ignore, as it makes it slower
refpt=[-9.8,150.5]; % [lat,long] of reference point for distance constraint
mindeg = 30; %minimum angular distance (0?degrees?180)
maxdeg = 140; %maximum angular distance (0?degrees?180)

%%  Get event information from HarvardCMT.mat file
load('harvardCMT.mat'); %loads the cmt structure

yes_Mw=find(cmt.Mw > Mwmin); %impose Mw constraint
yes_MS=find(cmt.MS >= MSmin); %impose MS constraint
yes_Mb=find(cmt.Mb >= Mbmin); %impose Mb constraint
yes_M0=find(cmt.M0 >= M0min); %impose M0 constraint
goodev=intersect(yes_Mw,yes_MS); goodev=intersect(goodev,yes_Mb); goodev=intersect(goodev,yes_M0);
% impose date window
yes_startdate=find(datenum(cmt.year,cmt.month,cmt.day) >= datenum(startdate));
yes_enddate=find(datenum(cmt.year,cmt.month,cmt.day) <= datenum(enddate));
goodev=intersect(goodev,yes_startdate);
goodev=intersect(goodev,yes_enddate);
% impose depth window
yes_mindepth=find(cmt.depth >= mindepth);
yes_maxdepth=find(cmt.depth <= maxdepth);
goodev=intersect(goodev,yes_mindepth);
goodev=intersect(goodev,yes_maxdepth);
% impose distance window
if exist('refpt','var')==1
    yes_mindeg=find(distance(refpt(1),refpt(2),cmt.lat(goodev),cmt.long(goodev)) >= mindeg);
    yes_maxdeg=find(distance(refpt(1),refpt(2),cmt.lat(goodev),cmt.long(goodev)) <= maxdeg);
    goodgood=intersect(yes_mindeg,yes_maxdeg);
else
    goodgood=[1:length(goodev)];
end

infoevs={'year','month','day','jjj','hour','minute','sec','lat','long','depth','Mb','MS','M0','Mw'};
for j=1:length(infoevs)
    assignin('base',char(infoevs(j)),eval(sprintf('cmt.%s(goodev(goodgood))',char(infoevs(j))))); %assign variables
    eval(sprintf('ev%s=%s;',char(infoevs(j)),char(infoevs(j)))); %append "ev" prefix to variables
end
%% Add rows to origin table
db=dbopen(db,'r+');
dborigin=dblookup_table(db,'origin');
dbrecord=dbaddnull(dborigin);
try
orids=unique(dbgetv(dborigin,'orid'));
evtimes=unique(dbgetv(dborigin,'time'));
catch
    orids=0; evtimes=0;
end
orid=max(orids);
if orid < 1; orid=0; end;
for ie=1:length(evdepth);
    noreps=any(abs(str2epoch(sprintf('%d/%d/%d %d:%d:%f',evmonth(ie),evday(ie),evyear(ie),evhour(ie),evminute(ie),evsec(ie)))-evtimes) < 50);
    if noreps==0 %only satisfied if no events within 100s of 'this' event
    orid=orid+1;
    db.record=dbaddv(dborigin,'lat',evlat(ie),'lon',evlong(ie),'depth',evdepth(ie),...
    'time',str2epoch(sprintf('%d/%d/%d %d:%d:%f',evmonth(ie),evday(ie),evyear(ie),evhour(ie),evminute(ie),evsec(ie))),...
    'orid',orid,'mb',evMb(ie),'ms',evMS(ie),'ml',evMw(ie),'nass',1,'ndef',1);
    end
end
dbcrunch(dborigin);
dbclose(db);
