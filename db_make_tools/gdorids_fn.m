function [goodorids, goodevts]=gdorids_fn(db,MSmin,Mbmin,startdate,enddate,mindepth,maxdepth,refpt,mindeg,maxdeg,minaz,maxaz,sta,datawind)
% [goodorids]=gdorids_fn(db,MSmin,Mbmin,startdate,enddate,mindepth,maxdepth,refpt,mindeg,maxdeg,minaz,maxaz,sta,datawind)
% makes a list of orids from an origin database that fit certain criteria
%
%
%% Set criteria for accepted events
% NB THESE VALUES ARE INCLUSIVE!! ( <= or >= )
% MAGNITUDE
% MSmin=6.5;
% Mbmin=0;
% TIME
% startdate=[start_y,start_m,start_d]; %beginning of time window, in form [yyyy,mm,dd]
% enddate=[end_y,end_m,end_d]; %end of time window, in form [yyyy,mm,dd]
% DEPTH
% mindepth=0; %minimum depth (km)
% maxdepth=1000; %maximum depth (km)
% DISTANCE 
% refpt=[reflat,reflon]; % [lat,long] of reference point for distance constraint
% mindeg=85; %minimum angular distance (0?degrees?180)
% maxdeg=140; %maximum angular distance (0?degrees?180)
% BACKAZ
% minaz=0;
% maxaz=360;
% ISDATA
% check_data=0;   % 0 (default) don't check
%                 % 1 do check        >> requires wfdisc table <<
% datawind=[1201 1800]; % seconds of data after event
% if ~exist('sta','var')==1
%     sta='J'; % only needed if checking data
% end

%% read in info and make variables
dbor=dblookup_table(db,'origin');
orids=dbgetv(dbor,'orid');
MS=dbgetv(dbor,'ms');
Mb=dbgetv(dbor,'mb');
times=dbgetv(dbor,'time');
depths=dbgetv(dbor,'depth');
lats=dbgetv(dbor,'lat');
lons=dbgetv(dbor,'lon');

%% apply criteria
yes_MS=find(MS >= MSmin | MS == -999); %impose MS constraint or null
yes_Mb=find(Mb >= Mbmin | Mb == -999); %impose Mb constraint or null
goodev=intersect(yes_Mb,yes_MS);
% impose date window
yes_startdate=find(times >= str2epoch(sprintf('%.0f/%.0f/%.0f 00:00',startdate)));
yes_enddate=find(times <= str2epoch(sprintf('%.0f/%.0f/%.0f 23:59',enddate)));
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
if nargin>17;
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

%% Make output structure

goodevts = struct([]);
for j = 1:length(goodorids)
    goodevts(1).year(j,1) = str2num(epoch2str(times(orids==goodorids(j)),'%Y'));
	goodevts(1).month(j,1) = str2num(epoch2str(times(orids==goodorids(j)),'%m'));
    goodevts(1).day(j,1) = str2num(epoch2str(times(orids==goodorids(j)),'%d'));
    goodevts(1).jjj(j,1) = str2num(epoch2str(times(orids==goodorids(j)),'%j'));
    goodevts(1).hour(j,1) = str2num(epoch2str(times(orids==goodorids(j)),'%H'));
    goodevts(1).minute(j,1) = str2num(epoch2str(times(orids==goodorids(j)),'%M'));
    goodevts(1).sec(j,1) = str2num(epoch2str(times(orids==goodorids(j)),'%S.%s'));
    goodevts(1).evtime = times(orids==goodorids);
    goodevts(1).lat(j,1) = lats(orids==goodorids(j));
    goodevts(1).long(j,1) = lons(orids==goodorids(j));
    goodevts(1).depth(j,1) = depths(orids==goodorids(j));
    goodevts(1).Mb(j,1) = Mb(orids==goodorids(j));
    goodevts(1).MS(j,1) = MS(orids==goodorids(j));
    goodevts(1).M0 = zeros(length(goodorids),1);
    goodevts(1).Mw = zeros(length(goodorids),1);
end
    
if isempty(goodorids)==1, fprintf('NO GOOD ORIDS\n'), end


end