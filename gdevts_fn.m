function [goodevts]=gdevts_fn(MSmin,Mwmin,startdate,enddate,mindepth,maxdepth,refpt,mindeg,maxdeg,minaz,maxaz)
%[goodevts]=gdevts_fn(MSmin,Mwmin,startdate,enddate,mindepth,maxdepth,refpt,mindeg,maxdeg,minaz,maxaz)
%   like gdorids_fn, but instead of getting events that satisfy criteria
%   from a database, it gets it from the cmt catalogue
%
%% Set criteria for accepted events
% NB THESE VALUES ARE INCLUSIVE!! ( <= or >= )
% MAGNITUDE
% MSmin=6.5;
% Mwmin=0;
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


%%  Get event information from HarvardCMT.mat file
load('/Users/Zach/Documents/MATLAB/SplitLab1.0.5/harvardCMT.mat'); %loads the cmt structure

yes_MS=find(cmt.MS > MSmin); %impose MS constraint
yes_Mw=find(cmt.Mw > Mwmin); %impose Mb constraint

goodev=intersect(yes_Mw,yes_MS); 
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
    yes_mindeg=find(angsep(refpt(1),refpt(2),cmt.lat(goodev),cmt.long(goodev)) >= mindeg);
    yes_maxdeg=find(angsep(refpt(1),refpt(2),cmt.lat(goodev),cmt.long(goodev)) <= maxdeg);
    goodgood=intersect(yes_mindeg,yes_maxdeg);
else
    goodgood=[1:length(goodev)];
end

infoevs={'year','month','day','jjj','hour','minute','sec','lat','long','depth','Mb','MS','M0','Mw'};
goodevts = struct([]);
for j=1:length(infoevs)
    eval(sprintf('goodevts(1).%s=cmt.%s(goodev(goodgood));',char(infoevs(j)),char(infoevs(j)))); %append "ev" prefix to variables
end
for j=1:length(goodgood)
goodevts.evtime(j) = str2epoch(sprintf('%u/%u/%u %u:%u:%.2f',goodevts.year(j),goodevts.month(j),goodevts.day(j),goodevts.hour(j),goodevts.minute(j),goodevts.sec(j)));
end
goodevts.evtime = reshape(goodevts.evtime,length(goodgood),1);

if isempty(goodevts)==1, fprintf('NO GOOD EVENTS\n'), end


end
