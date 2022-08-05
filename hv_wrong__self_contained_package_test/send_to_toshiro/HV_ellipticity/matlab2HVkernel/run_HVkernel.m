function [HVr,HVK,phV,grV] = run_HVkernel(model,swperiods,ID,ifdelete,ifplot,ifverbose)
% [HVrs,HVKs,phVs] = run_HVkernel(model,swperiods,ID,ifdelete,ifplot,ifverbose)
% 
% Function to run HV kernel script for a given model and extract the phase
% velocities and kernels at a bunch of input periods. 
% 
% Outputs are horizontal / vertical values, positive for retrograde motion.
% 

tic

if nargin < 3 || isempty(ID)
    ID = 'eg';
end
if nargin < 4 || isempty(ifdelete)
    ifdelete = true;
end
if nargin < 5 || isempty(ifplot)
    ifplot = false;
end
if nargin < 6 || isempty(ifverbose)
    ifverbose = true;
end



%% filenames
if ~ischar(ID), ID = num2str(ID);end
execfile = [ID,'.run_HVker'];
modfile = [ID,'.model'];
ofile = [ID,'.out'];
logfile = [ID,'.log'];


%% =======================================================================
wd = pwd;
cd('/Users/zeilon/Work/codes/HV_Tanimoto/matlab_to_HVkernel');


%% model is an input file or a matlab structure?
if ischar(model) && exist(model,'file')==2
    modfile = model;
    delcard = false;
else
	writeHVkernel_modfile(modfile,model.Z,model.Vp,model.Vs,model.rho,[],[],1./swperiods);
    delcard =true;
end

writeHVkernel_execfile( execfile,modfile,ofile,logfile);
system(['chmod u+x ' execfile]);


%% do MINEOS on it
if ifverbose
    fprintf('    > Running HV kernel computation code. \n    > Will take some time...')
end
[status,cmdout] = system(['/opt/local/bin/gtimeout 100 ./',execfile]);
if ifverbose
     fprintf(' success!\n')
end
%% read modes output
if status~=124
    try
    [HVr,HVK,phV,grV] = readHVkernel_ofile(ofile,swperiods,ifplot);
    HVr = -1./HVr;
    catch
        fprintf('some error - check model file layers not too thin!\n')
        error
    end
    phV = phV(:);
    grV = grV(:);
else 
    error('KVkernels did not finish in 100s')
end
    


%% delete files
if ifdelete
    delete(execfile,ofile);
    if exist(logfile,'file')==2, delete(logfile); end
    if delcard, delete(modfile); end
end
cd(wd);

% %% plot
% if ifplot
%     figure(88), clf; set(gcf,'pos',[331 385 848 613]);
%     ax1 = axes; hold on;
%     % dispersion curves
%     hd(1)=plot(ax1,swperiods,phV,'o-','linewidth',2);
%     hd(2)=plot(ax1,swperiods,grV,'o-','linewidth',2);
%     hl = legend(ax1,hd,{'Phase (c)','Group (U)'},'location','southeast');
%     set(hl,'fontsize',16,'fontweight','bold')
%     set(ax1,'fontsize',16)
%     xlabel(ax1,'Period (s)','interpreter','latex','fontsize',22)
%     ylabel(ax1,'Velocity (km/s)','interpreter','latex','fontsize',22)
% end
if ifverbose
    fprintf('HV   %s took %.5f s\n',ID,toc)
end

 