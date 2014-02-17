function [data,traces,tts,nsamps,samprates,stas,chans,gcarcs,seazs,esazs,trefs,arids,evdata] = dbgetdata( db, varargin )
%   [data,traces,tts,nsamps,samprates,stas,chans,gcarcs,seazs,esazs,trefs] = dbgetdata( db, varargin )
%
%   Master function to extract data from a database - a synthesis of
%   previous dbgetwf, dbgetwf3, dbgetwfz functions, but doing a better job
%   using trload_css, trextract, and db2struct to get more information,
%   while also dealing with crossovers of wfdisc rows and data duplication
%
% There are 3 modes for data extraction, defined by the second argument:
%  1) Simple windowing by time, if the 2nd argument is ''time|t|epoch'
%       - in this case, you simply specify a start and end time and all
%       data within this window is extracted (+ subsetting by sta, chan)
%  2) Event windowing, if the 2nd argument is 'orid|evt|event|e'
%       - in this case, data associated with a given event is extracted.
%       This can be windowed with respect to the event time, or to a given
%       phase arrival (either using associated picks or the iasp91
%       prediction). Again, subsetting by sta, chan, supported. In this
%       mode, can only handle one event at a time.
%  3) Station windowing, if the 2nd argument is 'sta|s|station'
%       - in this case, data from all orids recorded at the given station
%       is extracted. If a list of events is given, then only data from
%       these events is extracted - in either case there is an option to
%       window about the event time or a given phase arrival ((either using
%       associated picks or the iasp91 prediction). This is useful for e.g.
%       active source record sections. Subsetting by chan is supported.
% 
% N.B. output data are detrended - i.e. best fit slope removed.
% 
% SIMPLE WINDOWING
% [data,traces,tts,nsamps,samprates] = dbgetdata(db,'time',[tstart tend],...)
% [...] = dbgetdata(db,'time',[tstart tend],...,'s|sta|station','sta_name',...)
% [...] = dbgetdata(db,'time',[tstart tend],...,'c|chan|channel','chan_name',...)
% EVENT WINDOWING
% [data,traces,tts,nsamps,samprates,stas,chans,gcarcs,seazs,esazs,trefs,arids,evdata] = dbgetdata(db,'evt',orid,...)
% [...] = dbgetdata(db,'evt',orid,...,'w|wind|window',[pretime posttime],...)
% [...] = dbgetdata(db,'evt',orid,...,'p|ph|phs|phase','phase_name',...)
% [...] = dbgetdata(db,'evt',orid,...,'mode|imode','pred/pick',...)
% [...] = dbgetdata(db,'evt',orid,...,'s|sta|station','sta_name',...)
% [...] = dbgetdata(db,'evt',orid,...,'c|chan|channel','chan_name',...)
% STATION WINDOWING
% [data,traces,tts,nsamps,samprates,stas,chans,gcarcs,seazs,esazs,trefs,arids,evdata] = dbgetdata(db,'sta',STATION,...)
% [...] = dbgetdata(db,'sta',STATION,...,'e|evt|orid',[orids],...)
% [...] = dbgetdata(db,'sta',STATION,...,'w|wind|window',[pretime posttime],...)
% [...] = dbgetdata(db,'sta',STATION,...,'p|ph|phs|phase','phase_name',...)
% [...] = dbgetdata(db,'sta',STATION,...,'mode|imode','pred/pick',...)
% [...] = dbgetdata(db,'sta',STATION,...,'c|chan|channel','chan_name',...)
%
% Output information:
% - "data" actually has all the information in its structure other than
% event and station details
% - "traces" and "tts" are nsamp x ntr matrices, with each trace in a
% different column. These are related to the stations and channels by
% unfolding the nsta x nchan values so that ntr = nsta*nchan:
% i.e.
%   [sta1_chan1 sta1_chan2 sta1_chan3 sta2_chan1 sta2_chan2 ...]
%   [    .          .          .          .          .         ]
%   [    .          .          .          .          .         ]
%   [    .          .          .          .          .         ]
% 
% - by default, the stations will be in order found in origin table
% - "trefs" are predicted or picked phase times (depending on mode)
%
% Defaults:
%  sta - all
%  chan - all
%  window - 100s before and after time
%  imode - 0 (0 is predicted, 1 is picked, requires arrival + assoc tables)
%  phase - none, so picks around event time
%  events - all (when relevant, i.e. in station windowing mode) 
% 
% Written by  Zach Eilon, July 19th 2013

%% Prelim checks
% check nargin
if(mod(nargin-1,2))
    error('Unpaired option(s)!');
end
if (~exist('db', 'var'))
   disp('database db not yet existant; must be open as db=dbopen(dbname, ''r'') ');
   disp('needs at least assoc arrival wfdisc tables');
	return
end
if (isempty(dbquery(db,'dbTABLE_COUNT')))
	disp('no data in current database db or not open')
    cd(wd);
	return
end

% temporarily move to database root directory, so data paths are correct
wd = pwd;
[dbdir, dbnam] = fileparts(dbquery(db,'dbDATABASE_NAME'));
if isempty(dbdir), dbdir = '.'; end
cd(dbdir);

%% Parse inputs
% defaults
tstart = 0;
tend = 0;
sta = '';
chan = '';
pretime = 100;
posttime = 100;
phase = '';
imode = 0; % 0 means predicted, 1 means picked
orids_do = 0;

for i = 4:2:nargin
if(isempty(varargin{i})); continue; end
    switch lower(varargin{i-1})
        case {'sta' 'station' 's'}
            sta = varargin{i};
        case {'chan' 'channel' 'c'}
            chan = varargin{i};
        case {'e' 'evt' 'orid' 'event' 'orids' 'evts' 'events'}
            orids_do = varargin{i};
        case {'wind' 'window' 'w'}
            pretime = varargin{i}(1);
            posttime = varargin{i}(end);
        case {'ph' 'phase' 'p' 'phs' 'phas'}
            phase = varargin{i};
        case {'imode','mode'}
            imode = varargin{i};
    end
end

%% absolute time mode
if any(strcmp(varargin{1},{'time','t','epoch'}))
wmode = 0; % absolute time mode

if length(varargin{2})~=2
	error('Must provide start AND end times')
else
    tstart = varargin{2}(1);
    tend = varargin{2}(end);
end

%load necessary tables
dbwf = dblookup_table(db,'wfdisc');

%subset to inputs
if ~isempty(sta)
    dbwf = dbsubset(dbwf,sprintf('sta == "%s"',sta));
end
if ~isempty(chan)
    dbwf = dbsubset(dbwf,sprintf('chan == "%s"',chan));
end
% stas = unique(dbgetv(dbwf,'sta'));
chans = unique(dbgetv(dbwf,'chan'));

% load trace object
try
    tr = trload_css(dbwf,tstart,tend);
catch
    fprintf('No data in this window for %s/n',sta);
    return
end

% calibrate
trapply_calib(tr);

samprate = unique(dbgetv(tr,'samprate')); 
if length(samprate)>1, error('Samprate non unique across selection.'); end

% splice to deal with overlaps
trsplice(tr,1/samprate)
dat = tr2struct(tr);

if nargout>1 % only do extra stuff if needed...
nsamp = max(dbgetv(tr,'nsamp')); 
ntr = dbnrecs(tr);
traces = NaN(nsamp,ntr);
tts = NaN(nsamp,ntr);
nsamps = zeros(ntr,1);
samprates = zeros(ntr,1);
for it = 1:ntr
tr.record = it-1;
trace = detrend(trextract_data(tr));
nsamps(it) = dbgetv(tr,'nsamp');
samprates(it) = dbgetv(tr,'samprate');
stas{it} = dbgetv(tr,'sta');
t0 = dbgetv(tr,'time');
% deal with offset in case less data in this trace
if nsamps(it)<nsamp
    samplag = round((t0-tstart)*samprate);
else
    samplag=0;
end

traces(1+samplag:samplag+length(trace),it)=trace;
tts(1+samplag:samplag+length(trace),it) = t0+([0:nsamps(it)-1]'./samprates(it));
end
end
trdestroy(tr);

stas = unique(stas);

if nargout>7
fprintf('Simple time window mode does not sustain gcarc,seaz,esaz,tref output arguments\n')
end
%% event mode
elseif any(strcmp(varargin{1},{'orid','evt','event','e'}))
wmode = 1; % event mode

orid = varargin{2};
if ischar(orid), orid = eval(orid); end % make sure it's a number

%load necessary tables
dbwf = dblookup_table(db,'wfdisc');
dbor = dblookup_table(db,'origin');
dbsi = dblookup_table(db,'site');
dbar = dblookup_table(db,'arrival');
dbas = dblookup_table(db,'assoc');

%event details
dbors = dbsubset(dbor,sprintf('orid == "%u"',orid));
evdata = struct([]);
[evdata(1).elat,evdata(1).elon,evdata(1).edep,evdata(1).evtime,evdata(1).mag] = dbgetv(dbors,'lat','lon','depth','time','ms');

%subset to inputs
if ~isempty(sta)
    dbsi = dbsubset(dbsi,sprintf('sta == "%s"',sta));
end
if ~isempty(chan)
    dbwf = dbsubset(dbwf,sprintf('chan == "%s"',chan));
end

% get station info and loop
[stas,slats,slons] = dbgetv(dbsi,'sta','lat','lon');
dat = struct([]);
trefs = zeros(size(stas,1),1);
arids = zeros(size(stas,1),1);
for is = 1:size(stas,1)
    sta = char(stas(is,:));
    dbwfs = dbsubset(dbwf,sprintf('sta == "%s"',sta));
% Get time window    
    if isempty(phase) % window from event time
        tstart = evdata(1).evtime - pretime;
        tend   = evdata(1).evtime + posttime;
    else % window from 1st arrival of named phase
        if ~imode % if predicted
        taupt = tauptime('ph',phase,'dep',evdata(1).edep,'sta',[slats(is),slons(is)],'evt',[evdata(1).elat evdata(1).elon]);
        if isempty(taupt)
            fprintf('No predicted %s arrival at %s, orid %u\n',phase,sta,orid)
            continue
        end
        tstart = taupt(1).time + evdata(1).evtime - pretime;
        tend   = taupt(1).time + evdata(1).evtime + posttime;
        elseif imode % if picked
        dbsas = dbsubset(dbas,sprintf('sta == "%s" && orid == %u && phase == "%s"',sta,orid,phase));
        
        if ~dbnrecs(dbsas)
            fprintf('No %s assoc for %s, orid %u\n',phase,sta,orid)
            continue
        end
        
        dbjsasr = dbjoin(dbar,dbsas);
        
        if dbnrecs(dbjsasr)>1
            fprintf('Multiple %s assocs for %s, orid %u... using last one\n',phase,sta,orid)
            dbjsasr.record=dbnrecs(dbjsasr)-1;
        end
        
        tstart = dbgetv(dbjsasr,'time') - pretime;
        tend   = dbgetv(dbjsasr,'time') + posttime;
        arids(is) = dbgetv(dbjsasr,'arid');
        end
    end

% load trace object
try
    tr = trload_css(dbwfs,tstart,tend);
catch
    fprintf('%s has no data in this window\n',sta);
    continue
end

% calibrate
trapply_calib(tr);

samprate = unique(dbgetv(tr,'samprate')); 
if length(samprate)>1, error('Samprate non unique across selection.'); end

% splice to deal with overlaps
trsplice(tr,1/samprate)

if dbgetv(tr,'nsamp')~=samprate*(pretime+posttime)
    fprintf('tr from %s has only %u samples\n',sta,dbgetv(tr,'nsamp'))
end

if isempty(dat)
dat = tr2struct(tr); 
else
dat(size(dat,1)+1,:) = tr2struct(tr); 
end
trdestroy(tr);

% Now can find tref, as a trace has definitely been loaded
trefs(is) = tstart+pretime;

end %loop on stas

if isempty(dat),
    fprintf('\nNo data exists with your specifications\n')
    cd(wd);
return
end
    
if nargout>1 % only do extra stuff if needed...
    nsamp = samprate*(tend-tstart); 
    ntr = numel(dat);
    nsta = size(dat,1);
    nchan = size(dat,2);
    
    traces = NaN(nsamp,ntr);
    tts = NaN(nsamp,ntr);
    nsamps = zeros(nsta,nchan);
    samprates = zeros(nsta,nchan);
    stas = cell(nsta,1);
    chans = cell(nchan,1);
    for is = 1:nsta
        stas{is} = dat(is,1).sta;
    for ic = 1:nchan
        chans{ic} = dat(1,ic).chan;
        nsamps(is,ic) = dat(is,ic).nsamp;
        samprates(is,ic) = dat(is,ic).samprate;
        traces(1:nsamps(is,ic),ic + (is-1)*size(dat,2)) = detrend(dat(is,ic).data);
        tts(1:nsamps(is,ic),ic + (is-1)*size(dat,2)) = [dat(is,ic).time:1/dat(is,ic).samprate:dat(is,ic).endtime]';
    end
    end
    [gcarcs,seazs] = distance(slats,slons,evdata(1).elat,evdata(1).elon);
    [~,esazs]      = distance(evdata(1).elat,evdata(1).elon,slats,slons);
    gcarcs = gcarcs(trefs~=0);
    seazs = seazs(trefs~=0);
    esazs = esazs(trefs~=0);
    arids = arids(trefs~=0);
end % if nargout>1
    trefs = trefs(trefs~=0);  
%% station mode
elseif any(strcmp(varargin{1},{'s','sta','station'}))
wmode = 2; % station mode

sta = varargin{2};
stas=sta;

%load necessary tables
dbwf = dblookup_table(db,'wfdisc');
dbor = dblookup_table(db,'origin');
dbsi = dblookup_table(db,'site');
dbar = dblookup_table(db,'arrival');
dbas = dblookup_table(db,'assoc');

%sta details
dbsis = dbsubset(dbsi,sprintf('sta == "%s"',sta));
dbwf = dbsubset(dbwf,sprintf('sta == "%s"',sta));
dbar = dbsubset(dbar,sprintf('sta == "%s"',sta));
dbas = dbsubset(dbas,sprintf('sta == "%s"',sta));
if dbnrecs(dbsis)~=1, error('No matching station in database'), end
[slat,slon] = dbgetv(dbsis,'lat','lon');

% get event info and loop
try
    [orids,elats,elons,edeps,evtimes,emags] = dbgetv(dbor,'orid','lat','lon','depth','time','mw');
catch
	[orids,elats,elons,edeps,evtimes,emags] = dbgetv(dbor,'orid','lat','lon','depth','time','mb');
end

%subset to inputs
if ~isempty(chan)
    dbwf = dbsubset(dbwf,sprintf('chan == "%s"',chan));
end

% if no orids specified, do all
if ~orids_do
    orids_do = orids;
end

norid = length(orids_do);

dat = struct([]);
trefs = zeros(norid,1);
arids = zeros(norid,1);

for ie = 1:norid
    orid = orids_do(ie);
    iorid = find(orids==orid);
%     dbors = dbsubset(dbor,sprintf('orid == "%u"',orid));
%     evdata = struct([]);
%     [evdata(1).elat,evdata(1).elon,evdata(1).edep,evdata(1).evtime] = dbgetv(dbors,'lat','lon','depth','time');

%     dbwfs = dbsubset(dbwf,sprintf('sta == "%s"',sta));
    
    if isempty(phase) % window from event time
        tstart = evtimes(iorid) - pretime;
        tend   = evtimes(iorid) + posttime;
    else % window from 1st arrival of named phase
        if ~imode % if predicted
        taupt = tauptime('ph',phase,'dep',edeps(iorid),'sta',[slat,slon],'evt',[elats(iorid) elons(iorid)]);
        if isempty(taupt)
            fprintf('No predicted %s arrival for orid %u\n',phase,orid)
            continue
        end
        tstart = taupt(1).time + evtimes(iorid) - pretime;
        tend   = taupt(1).time + evtimes(iorid) + posttime;
        elseif imode % if picked
        dbsas = dbsubset(dbas,sprintf('orid == %u && phase == "%s"',orid,phase));
        if ~dbnrecs(dbsas)
            fprintf('No %s assoc for %s, orid %u\n',phase,sta,orid)
            continue
        end
        dbjsasr = dbjoin(dbar,dbsas);
        tstart = dbgetv(dbjsasr,'time') - pretime;
        tend   = dbgetv(dbjsasr,'time') + posttime;
        arids(ie) = dbgetv(dbjsasr,'arid');
        end
    end

% load trace object
try
    tr = trload_css(dbwf,tstart,tend);
catch
    fprintf('%s has no data in this window for orid %u\n',sta,orid);
    continue
end


% calibrate
trapply_calib(tr);

samprate = unique(dbgetv(tr,'samprate')); 
if length(samprate)>1, error('Samprate non unique across selection.'); end

% splice to deal with overlaps
trsplice(tr,1/samprate)

if dbgetv(tr,'nsamp')~=samprate*(pretime+posttime)
    fprintf('tr from orid %u has only %u samples\n',orid,dbgetv(tr,'nsamp'))
end

if isempty(dat)
dat = tr2struct(tr); 
else
dat(size(dat,1)+1,:) = tr2struct(tr); 
end
trdestroy(tr);

% Now can find tref, as a trace has definitely been loaded
trefs(ie) = tstart+pretime;
    
end %loop on orids

if isempty(dat),
    fprintf('\nNo data exists with your specifications\n')
    cd(wd);
    return
end
    
if nargout>1 % only do extra stuff if needed...
    nsamp = samprate*(tend-tstart); 
    ntr = numel(dat);
    norid = size(dat,1);
    nchan = size(dat,2);
    
    traces = NaN(nsamp,ntr);
    tts = NaN(nsamp,ntr);
    nsamps = zeros(norid,nchan);
    samprates = zeros(norid,nchan);
    chans = cell(nchan,1);
    for ie = 1:norid
    for ic = 1:nchan
        chans{ic} = dat(1,ic).chan;
        nsamps(ie,ic) = dat(ie,ic).nsamp;
        samprates(ie,ic) = dat(ie,ic).samprate;
        traces(1:nsamps(ie,ic),ic + (ie-1)*size(dat,2)) = detrend(dat(ie,ic).data);
        tts(1:nsamps(ie,ic),ic + (ie-1)*size(dat,2)) = [dat(ie,ic).time:1/dat(ie,ic).samprate:dat(ie,ic).endtime]';
    end
    end
    [gcarcs,seazs] = distance(slat,slon,elats,elons);
    [~,esazs]      = distance(elats,elons,slat,slon);
    gcarcs = gcarcs(trefs~=0);
    seazs = seazs(trefs~=0);
    esazs = esazs(trefs~=0);
    arids = arids(trefs~=0);
    evdata = struct('orid',orids(trefs~=0),'elat',elats(trefs~=0),...
                    'elon',elons(trefs~=0),'edep',edeps(trefs~=0),...
                    'evtime',evtimes(trefs~=0),'mag',emags(trefs~=0));
end % if nargout>1
    trefs = trefs(trefs~=0);
end % excerpt mode - time, orid, or station based.

data = dat;

cd(wd);


end % on function

function disp(x)
    return
end

