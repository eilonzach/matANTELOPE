% This script goes through an origin table and an arrival table and seeks
% to build an association table - NB. not sure why this isn't working in
% dbpick!

%% Starts
dbnam = 'pngbodydb';
dbdir = '/Users/Zach/Documents/MATLAB/PNG_tomog/PNGbodydb/';% remember final "/"
dbtabdir = '/Users/Zach/Documents/MATLAB/PNG_tomog/PNGbodydb/'; % remember final "/"
maxTdiff = 25; % time in seconds to search around predicted arrival
phases = {'P','S','PP','SP','SKS','SKP','PKP','ScS','sS','pS'}';
ph_cutoff = 0; 

db = dbopen(strcat(dbdir,dbnam),'r');
dbsite = dblookup_table(db,'site');
dbor = dblookup_table(db,'origin');
dbarr = dblookup_table(db,'arrival');
dbass = dblookup_table(db,'assoc');

% see if assoc table already exists, if so move it to assoc_old
if dbnrecs(dbass) > 0
    yn = 'n'; %by default
    yn = input('Assoc table seems to exist. Do you want to overwrite? (y/n) ','s');
    if strcmp(yn,'y')==1,   movefile(strcat(dbtabdir,dbnam,'.assoc'),strcat(dbtabdir,dbnam,'.assoc_old')); 
    else return; end
end
% station details
nstas = dbnrecs(dbsite);
stas = dbgetv(dbsite,'sta');
slats=dbgetv(dbsite,'lat');
slons=dbgetv(dbsite,'lon');
% orid details
norids = dbnrecs(dbor);
orids = dbgetv(dbor,'orid');
elats=dbgetv(dbor,'lat');
elons=dbgetv(dbor,'lon');
edeps=dbgetv(dbor,'depth');
evtimes=dbgetv(dbor,'time');

% 3D matrix of predicted phase arrival times, "artimes"
% each row has different orid, each column is an different phase, each
% layer is a different station. Similarly, distance and seaz, but no
% dimension for the different phases
if exist('artimes','var')==0
    
artimes = zeros(norids,length(phases),nstas);
deltas = zeros(norids,nstas);
seazs = zeros(norids,nstas);
esazs = zeros(norids,nstas);
for is = 1:nstas
for ie = 1:norids
    [deltas(ie,is),seazs(ie,is)] = distance(slats(is),slons(is),elats(ie),elons(ie));
    [~,            esazs(ie,is)] = distance(elats(ie),elons(ie),slats(is),slons(is));
for ip = 1:length(phases)
%     if deltas(ie,is) > ph_cutoff && strcmp(char(phases(ip)),'S')==1, continue; end % don't do S above cutoff
%     if deltas(ie,is) < ph_cutoff && strcmp(char(phases(ip)),'S')~=1, continue; end % don't do core phases below cutoff
    try % get predicted arrival times - try because some phases don't exist
%   evalc('tauptt = tauptime(edeps(ie),char(phases(ip)),''deg'',deltas(ie,is));');
%     evalc('tauptt = tauptime(''dep'',edeps(ie),''phases'',char(phases(ip)),''sta'',[slats(is) slons(is)],''evt'',[elats(ie) elons(ie)]);'); % suppress output
    evalc('tauptt = taupant(''dep'',edeps(ie),''phases'',char(phases(ip)),''sta'',[slats(is) slons(is)],''evt'',[elats(ie) elons(ie)]);'); % suppress output
    artimes(ie,ip,is) = tauptt(1).time + evtimes(ie); % epochal
    catch 
        if strcmp(char(phases(ip)),'P')
        try
%             evalc('tauptt = tauptime(''dep'',edeps(ie),''phases'',''Pdiff'',''sta'',[slats(is) slons(is)],''evt'',[elats(ie) elons(ie)]);'); % suppress output
            evalc('tauptt = taupant(''dep'',edeps(ie),''phases'',''Pdiff'',''sta'',[slats(is) slons(is)],''evt'',[elats(ie) elons(ie)]);'); % suppress output
            artimes(ie,ip,is) = tauptt(1).time + evtimes(ie); % epochal
        end
        elseif strcmp(char(phases(ip)),'S')
        try
%             evalc('tauptt = tauptime(''dep'',edeps(ie),''phases'',''Sdiff'',''sta'',[slats(is) slons(is)],''evt'',[elats(ie) elons(ie)]);'); % suppress output
            evalc('tauptt = taupant(''dep'',edeps(ie),''phases'',''Sdiff'',''sta'',[slats(is) slons(is)],''evt'',[elats(ie) elons(ie)]);'); % suppress output
            artimes(ie,ip,is) = tauptt(1).time + evtimes(ie); % epochal
        end
        end
    end
end % loop on phases
end % loop on orids
fprintf('Station %s: done - %.1f%% \n',char(stas(is)),100*(is/nstas));
end % loop on stas
end % if no artimes variable - don't want to do this more than once!

dbo = dbopen(strcat(dbdir,dbnam,'_new'),'r+');
dboassoc = dblookup_table(dbo,'assoc');
dbrecord = dbaddnull(dboassoc);
% now loop through arrivals, trying to associate
narrs = dbnrecs(dbarr);
arids = dbgetv(dbarr,'arid');
for ia = 1:narrs
    if any(ia==0:round(narrs/20):narrs), fprintf('%.0f%% done\n',100*ia/narrs); end
    dbarr.record = dbfind(dbarr,sprintf('arid == %u',arids(ia)));
    if any(strcmp(phases,dbgetv(dbarr,'iphase'))) % only if phase is one of those to look for
    is = find(strcmp(dbgetv(dbarr,'sta'),stas)==1); % station index
    ip = find(strcmp(dbgetv(dbarr,'iphase'),phases)==1); % phase index
    artime = dbgetv(dbarr,'time');
    [tdiff, ie] = min(abs((artimes(:,ip,is)-artime)));
    if tdiff < maxTdiff % if time discrepancy is smaller than crit
    % put new row into assoc table.
    dbo.record=dbaddv(dboassoc,'arid',arids(ia),'orid',orids(ie),'sta',char(stas(is)),...
    'phase',char(phases(ip)),'delta',deltas(ie,is),'seaz',seazs(ie,is),...
    'esaz',esazs(ie,is),'timeres',artimes(ie,ip,is)-artime,'azres',0);
    else fprintf('No assoc arid #%u, orid %3u %4s %4s %s (tdiff %.2f) \n',arids(ia),...
    orids(ie),dbgetv(dbarr,'iphase'),dbgetv(dbarr,'sta'),strydtime(artime),tdiff);
%         if tdiff > 100     
%         dbclose(dbo)
%         dbclose(db)
%         return
%         end
    end
    end 
end % loop on arrivals
dbcrunch(dboassoc);
       
dbclose(dbo)
dbclose(db)

movefile(strcat(dbdir,dbnam,'_new','.assoc'),strcat(dbdir,dbnam,'.assoc'));


