% makes a list of orids from an origin database that fit certain criteria

%% Set criteria for accepted events
% NB THESE VALUES ARE INCLUSIVE!! ( <= or >= )
% MAGNITUDE
MSmin=6.5;
Mbmin=0;
% TIME
startdate=[2010,03,01]; %beginning of time window, in form [yyyy,mm,dd]
enddate=[2011,08,01]; %end of time window, in form [yyyy,mm,dd]
% DEPTH
mindepth=00; %minimum depth (km)
maxdepth=1000; %maximum depth (km)
% DISTANCE 
refpt=[-10,150]; % [lat,long] of reference point for distance constraint
mindeg=25; %minimum angular distance (0?degrees?180)
maxdeg=135; %maximum angular distance (0?degrees?180)
% BACKAZ
minaz=0;
maxaz=360;
% ISDATA
check_data=0;   % 0 (default) don't check
                % 1 do check        >> requires wfdisc table <<
datawind=[1201 1800]; % seconds of data after event
if ~exist('sta','var')==1
    sta='J'; % only needed if checking data
end

%% database and other information
dbdir='/Volumes/Portadb/PORTAPNG/';
dbnam='portapng';

%% read in info and make variables
db=dbopen(strcat(dbdir,dbnam),'r');
dbor=dblookup_table(db,'origin');
orids=dbgetv(dbor,'orid');
MS=dbgetv(dbor,'ms');
Mb=dbgetv(dbor,'mb');
times=dbgetv(dbor,'time');
depths=dbgetv(dbor,'depth');
lats=dbgetv(dbor,'lat');
lons=dbgetv(dbor,'lon');
dbclose(db)

%% apply criteria
yes_MS=find(MS >= MSmin); %impose MS constraint
yes_Mb=find(Mb >= Mbmin); %impose Mb constraint
goodev=intersect(yes_Mb,yes_MS);
% impose date window
yes_startdate=find(times >= str2epoch(sprintf('%.0f/%.0f/%.0f 00:00',startdate)));
yes_enddate=find(times <= str2epoch(sprintf('%.0f/%.0f/%.0f 00:00',enddate)));
goodev=intersect(goodev,yes_startdate);
goodev=intersect(goodev,yes_enddate);
% impose depth window
yes_mindepth=find(depths >= mindepth);
yes_maxdepth=find(depths <= maxdepth);
goodev=intersect(goodev,yes_mindepth);
goodev=intersect(goodev,yes_maxdepth);
% impose distance window
if exist('refpt','var')==1
    yes_mindeg=find(distance(refpt(1),refpt(2),lats(goodev),lons(goodev)) >= mindeg);
    yes_maxdeg=find(distance(refpt(1),refpt(2),lats(goodev),lons(goodev)) <= maxdeg);
    yes_dist=intersect(yes_mindeg,yes_maxdeg);
    goodev=goodev(yes_dist);
end
% impose azimuth window
if exist('refpt','var')==1
    yes_minaz=find(azimuth(refpt(1),refpt(2),lats(goodev),lons(goodev)) >= minaz);
    yes_maxaz=find(azimuth(refpt(1),refpt(2),lats(goodev),lons(goodev)) <= maxaz);
    yes_az=intersect(yes_minaz,yes_maxaz);
    goodev=goodev(yes_az);
end    

goodorids=orids(goodev);
% impose data check
if check_data==1;
    db=dbopen(dbnam,'r');
    yes_data=zeros(length(goodorids),1);
    for ie=1:length(goodorids)
        try
        [tt, dat, chans, nsamps, samprate, wfids]...
            = dbgetwfz(db,sta,datawind(1)+times(goodorids(ie)),datawind(2)+times(goodorids(ie)),'epoch');
        if length(chans)>=3 && nsamps > 0.5*samprate*(datawind(2)-datawind(1)) % the half is there just in case there is a sample or so on either side missing
        yes_data(ie)=ie;
        end
        catch me
        end
    end
    dbclose(db)
    yes_data(yes_data==0)=[];
    goodorids=goodorids(yes_data);
end

if isempty(goodorids)==1, fprintf('NO GOOD ORIDS\n'), end

% While testing
% disp(goodorids)

%% clean up
clear yes_*
clear('orids','MS','Mb','times','depths','lats','lons');
clear('MSmin','Mbmin','startdate','enddate','mindepth','maxdepth','refpt','mindeg','maxdeg');
clear('goodgood','goodev');