% This script makes a large structure of all the data across the array for
% a given arrival
% For a given database, and a ref location, this script will filter events
% by distance and then output all possible orids
% For a chosen orid, the script will then loop through stations and load
% all the waveforms associated with the orid and put them into the
% cell "alldata" which has nsta + 1 structures:
% the first nsta structures contain the info and data for each station
% the last structure contains the evinfo

% Check to see if "goodorids" needs to be run
if exist('goodorids','var')==1
yn=input('Regenerate orids list? (e.g. change parameters) (y/n) ','s');
if yn=='y'
run('/Users/Zach/Documents/MATLAB/seis_tools/gdorids');
end
else
run('/Users/Zach/Documents/MATLAB/seis_tools/gdorids');
end

%% starting info
dbdir='/Volumes/Portadb/PORTAPNG/'; % nb need the / at the end
dbnam='portapng';
rlat=-10;
rlon=150;
buffer=1; % for dbgetwfz wfdisc join
samprate=50;
filtfreqs=[0.02 0.125];
ndec=4;   % Decimation factor?

%% windowing: 
windopt='time'; % 'time' or 'arrival'
% excerpts by time
windstart = 50; %seconds after EQ to open window
windend = 1750;   %seconds after EQ to close window
%OR around arrival
phase='PS'; % phase to window around
pretime=20;
posttime=40;

if strcmp(windopt,'time')==1, windl=windend-windstart-2*buffer;
elseif strcmp(windopt,'arrival')==1, windl=posttime+pretime;
end

%% inputs
if exist('filtfreqs','var')==0;
yn='y';
else
% yn='n';
fprintf('Current filtfreqs:\n[%.3f %.3f]\ti.e. %.1fs to %.1fs\n',...
filtfreqs(1),filtfreqs(2),1/filtfreqs(1),1/filtfreqs(2));
yn=input('New filter? (y/n): ','s');
end
if yn=='y'
fprintf('Suggested filter: [0.02 0.125]\n');   
minfreq=input('Minimum frequency: ');
maxfreq=input('Maximum frequency: ');
filtfreqs=[minfreq maxfreq];
end


%% get all event info

db=dbopen(strcat(dbdir,dbnam),'r');

dbsite=dblookup_table(db,'site');
dbor=dblookup_table(db,'origin');
dbwf=dblookup_table(db,'wfdisc');

% nstas=dbquery(dbsite,'dbRECORD_COUNT');
% stas=dbgetv(dbsite,'sta');
% slats=dbgetv(dbsite,'lat');
% slons=dbgetv(dbsite,'lon');

norids=length(goodorids);
orids=goodorids;
elats=zeros(size(orids));
elons=zeros(size(orids));
depths=zeros(size(orids));
evtimes=zeros(size(orids));
seazs=zeros(size(orids));
gcarcs=zeros(size(orids));

dbjos=dbjoin(dbor,dbsite);
dbjows=dbjoin(dbjos, dbwf,{'sta', 'sarrival()'}, ...
	{'sta','time::endtime'});

fprintf(' n  orid    date       time        ~delta   ~seaz  depth  mb\n');
  for ie=1:norids
    ii=dbfind(dbor,sprintf('orid==%d',orids(ie)));
    dboo=dbor;
    dboo.record=ii;
    elats(ie)=dbgetv(dboo,'lat');
    elons(ie)=dbgetv(dboo,'lon');
	depths(ie)=dbgetv(dboo,'depth');
    evtimes(ie)=unique(dbgetv(dboo,'time'));
    seazs(ie)=azimuth(rlat,rlon,elats(ie),elons(ie));
    gcarcs(ie)=distance(rlat,rlon,elats(ie),elons(ie));
    fprintf('%2d  %2d   %s   %5.1f   %5.1f   %3.0f   %3.1f\n',ie,orids(ie),...
		 strtime(evtimes(ie)),...
		 gcarcs(ie),seazs(ie),depths(ie),dbgetv(dboo,'mb'));
  end 
%% set orid and loop through stations getting data
orid=input('Choose event orid to get data for: ');
elat=elats(orids==orid);
elon=elons(orids==orid);
edep=depths(orids==orid);
evtime=evtimes(orids==orid);

dbjowso=dbsubset(dbjows,sprintf('orid==%d',orid));

stas=unique(dbgetv(dbjowso,'sta'));
nstas=length(stas);

%% set up structure for all data - station data/info goes in the first nstas rows and
% the last row is for a structure with the orid information
data=struct('orid',orid,...
            'elat',elat,...
            'elon',elon,...
            'edep',edep,...
            'evtime',evtime,...
          'sta',cell(nstas,1),...
          'slat',zeros(nstas,1),...
          'slon',zeros(nstas,1),...
          'selev',zeros(nstas,1),...
          'seaz',zeros(nstas,1),...
          'foraz',zeros(nstas,1),...
          'gcarc',zeros(nstas,1),...
          'tt0',zeros(nstas,1),...
          'tte',zeros(nstas,samprate*windl),...
          'samprate',zeros(nstas,1),...
          'nsamps',zeros(nstas,1),...
          'dof',zeros(nstas,1),...
          'dt',zeros(nstas,1),...
          'datN',zeros(nstas,samprate*windl/ndec),...
          'datE',zeros(nstas,samprate*windl/ndec),...
          'datR',zeros(nstas,samprate*windl/ndec),...
          'datT',zeros(nstas,samprate*windl/ndec),...
          'datZ',zeros(nstas,samprate*windl/ndec));
data(1).filtfreqs = filtfreqs;
      
% Making these this big seems unecessary - CHECK
for is=1:nstas
sta=char(stas(is));
dbsta=dbsubset(dbsite,sprintf('sta=="%s"',sta));
slat=dbgetv(dbsta,'lat');
slon=dbgetv(dbsta,'lon');
selev=dbgetv(dbsta,'elev');
seaz=azimuth(slat,slon,elat,elon);
foraz=azimuth(elat,elon,slat,slon);
gcarc=distance(slat,slon,elat,elon);
fprintf('Processing sta %s orid %d\n',sta,orid);

if strcmp(windopt,'time')==1, 
    t0=evtime+windstart+buffer;
    t1=evtime+windend-buffer;
elseif strcmp(windopt,'arrival')==1, 
    tt=taupTime([],edep,phase,'sta',[slat,slon],'evt',[elat elon]);
    arrtime=tt(1).time;
    t0=evtime+arrtime-pretime;
    t1=evtime+arrtime+posttime;
end

%2.  Get records, rotate, preview
% %% get details for name of file for event
% yyyy=epoch2str(t0,'%Y');
% jjj=epoch2str(t0,'%j');
% hh=epoch2str(t0,'%H');
% mm=epoch2str(t0,'%M');
% ss=epoch2str(t0,'%S');
%% get data
try
[tt, dat,chans,nsamps,samprate, wfids] = dbgetwfz(db,sta,t0,t1,'epoch');
catch me
chans=0;
end
chanend=char({'e','n','z'});
nchan=size(chans,1);
% only continue if all  channels exist
if nchan==3; 
else data(is).sta='mark'; % mark for deletion
    continue
end % if station exists

%% Data info
dat=detrend(dat);
ampl=max(max(dat));
dat=dat ./ ampl;
%%%% taper - window...
dt = 1/samprate;
tte = tt - evtime; % time since event
% nchan=length(chans);
saz=sin(pi-seaz*pi/180); % pi shift to rotate so radial is in foraz
caz=cos(pi-seaz*pi/180);

% Things I want in the structure...
data(is).sta=sta;
data(is).slat=slat;
data(is).slon=slon;
data(is).selev=selev;
data(is).seaz=seaz;
data(is).foraz=foraz;
data(is).gcarc=gcarc;
cd /Users/Zach/Documents/MATLAB/PNG_swsplit/

%% Data processing on three channels
for k=1:nchan
eval(sprintf('data_%s=dat(:,k);',chanend(k)));
end

% FILTER %%% filt filt?
[b,a]=butter(3,filtfreqs.*2.*dt);
% data_e=filter(b,a,data_e);
% data_n=filter(b,a,data_n);
% data_z=filter(b,a,data_z);
data_e=filtfilt(b,a,data_e);
data_n=filtfilt(b,a,data_n);
data_z=filtfilt(b,a,data_z);
amp2=max([max(abs(data_e)),max(abs(data_n)),max(abs(data_z))]);
data_e=data_e./amp2;
data_n=data_n./amp2;
data_z=data_z./amp2;
tt0=tt;

% WINDOW - taper
we=window(@tukeywin,length(data_e));
wn=window(@tukeywin,length(data_n));
wz=window(@tukeywin,length(data_z));
data_e=data_e.*we;
data_n=data_n.*wn;
data_z=data_z.*wz;

%  DECIMATE here after filtering
if (ndec > 1)
  nold=length(data_z);
  rsmp=1:ndec:nold;
  tt0=tt0(rsmp);
  dtold=dt;
  dt=mean(diff(tt0));
  data_e=data_e(rsmp);
  data_n=data_n(rsmp);
  data_z=data_z(rsmp);
  nsamps=nsamps/ndec;
end

% ROTATE (T positive to the left!) (recall caz+saz have been rotated so
% angle is the forward azimuth of the arrival
data_r =  caz.*data_n + saz.*data_e;
data_t =  saz.*data_n - caz.*data_e;

chan=char({'N','E','R','T','Z'});
chanend=char({'n','e','r','t','z'});
nchan=size(chan,1);
dof=zeros(nchan,2); % work out degrees of freedom
for k=1:nchan
eval(sprintf('data(is).dat%s=data_%s;',chan(k),chanend(k)));
eval(sprintf('[dof(k,1),dof(k,2)]=scdofcalc(data_%s);',chanend(k)));% work out degrees of freedom
end
%assign to output structure
% eval(sprintf('%s.twind=[windstart-buffer windend+buffer];',sta));

data(is).nsamps=nsamps;
data(is).tt0=tt0;
data(is).tte=tt0-evtime;
data(is).dt=dt;
data(is).dof=dof;
% %reset channels
% chan=char({'e','n','z'});
% nchan=size(chan,1);
% assign to alldata
% eval(sprintf('alldata{j}=%s;',sta));
% eval(sprintf('clear(''%s'');',sta));
end % loop on stations
%% delete null rows
marked=0;
for is=1:nstas
    if strcmp(data(is).sta,'mark')==1, marked(length(marked)+1)=is; end
end
data(marked(2:end))=[]; clear('marked');
dbclose(db)

