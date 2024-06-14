function [traces,tt,status,cmdout] = run_propmat(LAYmodel,ID,ph,samprate,inc,synthperiod,nsamps,cutf,sourc)
% Most argument checks are done in run_propmat_or_telewavesim.m 
if isempty(samprate); 
    samprate = 5; 
end
if isempty(nsamps); 
    nsamps = 2^11; 
end
if isempty(inc); 
    inc = 5; 
end 
if ~all(unique(factor(nsamps))==2)
    error('Nsamps must be some power of 2')
end
if isempty(cutf)
    cutf = ceil((4 / synthperiod) * (nsamps / samprate));
end

% Function to run the propagator matrix code for a given layerised model. 
% demoPlot = false; % whether to plot the propmat results
demoPlot = strcmp(ph, 'Ps'); 

paths = getPaths(); 

% remaining parms
obsdist = 0;
ocomps = 2; % 1 is [x,y,z], 2 is [r,t,z]

%% filenames
if ~ischar(ID), ID = num2str(ID);end
modfile =    [ ID,'.mod'];
execfile =   [ ID,'.csh']; % brb2022.03.22 trying .csh seems to be a faster extension when runing from csh shell (csh shebang.). 
odatfile =   [ ID,'.rtz'];
ifile =      [ ID,'_synth.in'];
ofile0 =     [ ID,'_synth.out0'];
ofile1 =     [ ID,'_synth.out1'];
ofile2 =     [ ID,'_synth.out2'];
oimagout =   [ ID,'_imag.out'];

%% =======================================================================

if strcmp(ph,'Ps')
    Vbot = LAYmodel.Vp(end);
elseif strcmp(ph,'Sp')
    Vbot = LAYmodel.Vs(end);
end
nlay = LAYmodel.nlay;


%% write to PropMatrix format
writePROPMATmodfile( LAYmodel,modfile)
writePROPMATparmfile(ifile, Vbot, nlay+1,nsamps,samprate,cutf) % add one layer for the halfspace
if strcmp(sourc,'gauss')
    writePROPMATexecfile_gauss( execfile,modfile,ifile,ofile0,ofile1,ofile2,oimagout,odatfile,inc,ph,synthperiod,obsdist,ocomps)
elseif strcmp(sourc,'sine')
    writePROPMATexecfile(       execfile,modfile,ifile,ofile0,ofile1,ofile2,oimagout,odatfile,inc,ph,synthperiod,obsdist,ocomps)
end
fileattrib(execfile, '+x'); 

%% do PropMatrix on it
% warning('This is where things get slow')
[status ,cmdout ] = system([paths.timeout ' 5 ./',execfile]); % now not deleting imag file. 
[errorInfo]=assess_timeout_results(status, cmdout); 
% [status2,cmdout2] = system(sprintf('rm %s &',oimagout)); % Try deleting file in background? brb2022.03.25

% [status2,cmdout2] = system([paths.timeout ' 0.3 echo hi world']); % brb2022.03.22 just for testing things
% % system(sprintf('cp %s slow.%s',execfile, execfile)); 
% % fprintf('\nrun propmat: copying propmat exefile to slow*, and only tieout .3 s\n')

% [res ,cmdout ,err ] = jsystem([paths.timeout ' 0.3 ./',execfile], '/bin/bash'); 
% [res ,cmdout ,err ] = jsystem(['./',execfile], '/bin/bash -c'); 
% [a,b]=system('thing that no work')
% [res2,cmdout2,err2] = jsystem([paths.timeout ' 0.3 echo hi world']); % brb2022.03.22 just for testing slowness
% [res2,cmdout2,err2] = jsystem(['echo $PATH']); % brb2022.03.22 just for testing slowness

% p = java.lang.ProcessBuilder({paths.timeout,...
%     '0.3',...
%     ['./' execfile] }).start();
% reader = java.io.BufferedReader(java.io.InputStreamReader(p.getInputStream()));
% str = char(java.util.Scanner(reader).useDelimiter('\A').next());

% system(sprintf('cp %s slow.%s',execfile, execfile)); 
% fprintf('\nrun propmat: copying propmat exefile to slow*, and only tieout .3 s\n')

%% read PropMatrix output
try 
    [traces,tt] = readPROPMATtr(odatfile); %brb2023/05/26 Having problems with this. I've had problems compiling Fortran lately on this Mac, maybe that is the problem? 
    delete(modfile,execfile,odatfile,ifile,ofile0,ofile1,ofile2);%% delete files
catch e 
    delete(modfile,execfile,odatfile,ifile,ofile0,ofile1,ofile2);%% delete files
    warning('Cant read propmat files. Did timeout propmat not finish?')
    error(getReport(e))
end




if demoPlot; 
    % plot
    figure(2); clf, hold on
    comps = {'VERTICAL','RADIAL','TRANSVERSE'}; 
    traces2 = traces(:,[3,1,2]);
    for ip = 1:3
        subplot(3,1,ip)
        plot(tt,traces2(:,ip),'Linewidth',1.5)
        xlim([0 max(tt)]);
        ylabel(comps{ip},'fontsize',19,'fontweight','bold')
    end
    set(gcf,'position',[680         273        1058         825])
end
% fprintf('Propmat %s%s took %.5f s\n',ID,ph,toc)

