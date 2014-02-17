% Script to make a wfdisc table from a directory with sac files from an
% IRIS data request.
% 
% Requires any mseed files to have been converted to SAC via mseed2sac
%
% Expects sac files in the form:
% net.sta.chan.M.yyyy,jjj,hh:mm:ss.SAC
% NB odd punctuation is important!! see textscan line below

datadir='/Users/Zach/Documents/MATLAB/AK_swsplit/data/'; % directory with sac files
dbdir = '/Users/Zach/Documents/MATLAB/AK_swsplit/'; % database directory
dbnam = 'AKcheck'; % database name

D=dir(datadir);
nfiles=length(D);
x=0;
for ii=1:nfiles
    if isempty(regexp(D(ii).name,'.SAC','matchcase'))~=1 % if '.SAC' is found in the filename
        %% Get all details from file name
        junk = regexprep(D(ii).name,'\.',' ');
        C=textscan(junk,'%s %s %s %s %u,%u,%u:%u:%u %s');
        network = char(C{1});
        sta     = char(C{2});
        chan    = char(C{3}); ochan=lower(chan(end));
        year    = C{5};
        jjj     = C{6};
        hh      = C{7};
        mm      = C{8};
        ss      = C{9};  
        evtime=str2epoch(sprintf('%u/%u %u:%u:%u',C{5},C{6},C{7},C{8},C{9}));
%% Change sacfile name
        dfile=sprintf('%s.%s.sac.%s',epoch2str(evtime,'%Y.%j.%H.%M.%S'),sta,ochan);
        ofile=strcat(datadir,dfile);
        movefile(strcat(datadir,D(ii).name),strcat(datadir,dfile))
%% Read data
try
    dat=rsac(ofile);
catch meh
    dat=rsacsun(ofile);
end
delete(ofile) % delete sac file, so antelope can rewrite it. 

tt = dat(:,1);
dat=dat(:,2);
nsamps=length(dat);
dt=mean(diff(tt));
dt=0.00001*round(100000*dt); % round off to nearest 0.01 millisecs
samprate = 1./dt;
%% Make wfdisc entry
dbo=dbopen(strcat(dbdir,dbnam),'r+');
dbow=dblookup_table(dbo,'wfdisc');
dbr=dbow;
% Set up trace and its parameters
tr=trnew;
tr=dblookup_table(tr,'trace');
endtime= tr_endtime(evtime,samprate,nsamps);
%find wfdisc row corresponding to the dfile, if there is one
try
    dbr.record=dbfind(dbr,sprintf('dfile=="%s"',dfile));
if dbr.record~=-102; dbdelete(dbr); end % if null pointer no need to delete
catch meh
end

% Put the waveform into the trace-object:
tr.record=dbaddv(tr,'net', network,'sta',sta,'chan',chan,...
'nsamp',nsamps,'samprate',samprate,'time',evtime,'endtime', endtime);
trinsert_data(tr,dat);
% Save the trace data in a new database, with the underlying file in sac format:
trsave_wf(tr,dbow,'sc',ofile);
trdestroy( tr ) % Clean up
dbclose(dbo) 
    end
end
