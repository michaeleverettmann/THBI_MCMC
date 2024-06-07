function [phV,grV,eigfiles_fix] = run_mineos(model,swperiods,par_mineos,ifdelete,ifplot,ifverbose,maxrunN)
% Note: Cannot actually use the arguments - options approach since the
% function is already built to work using nargin. .... Matlab...
%     arguments 
%         model
%         swperiods
%         par_mineos
%         ifdelete
%         ifplot
%         ifverbose
%         options.maxrunN = 5e2 % This was default value Zach/Jon had. But when we are trying to initiate a model and it fails, no need to try 500 times before we decide it fails. 
%     end
% [phV,grV] = run_mineos(model,swperiods,par_mineos,ifdelete,ifplot,ifverbose)
% 
% Function to run the MINEOS for a given model and extract the phase
% velocities at a bunch of input periods. If you keep the output files
% around (ifdelete==false) then they can be used to calculate perturbation
% kernels with the complementary run_kernelcalc.m script

%%%%% Previous. Temporary. 
paths = getPaths(); 

tic1 = now;

if nargin < 3 || isempty(par_mineos)
    par_mineos = [];
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
if nargin < 7 || isempty(maxrunN) 
    maxrunN = 100; 
end

%% parameters
% default parameters
parm = struct('R_or_L','R',...
              'phV_or_grV','ph',...
              'ID','example',...
              'lmin',0,...            % minimum angular order
              'lmax',3500,...         % expected max angular order
              'fmin',0.05,...         % min frequency (mHz)
              'fmax',200.05,...       % max frequency (mHz) - gets reset by min period 
              'l_increment_standard',2,... % 
              'l_increment_failed',5,...
              'maxrunN',maxrunN,... % bb2021.09.21 changing temporarily from 5e2 to 5e1. Why would we want 500 failed iterations? That takes forever. 
              'qmodpath',[paths.THBIpath '/matlab_to_mineos/safekeeping/qmod']); % bb2021.09.14 if matlab_to_mineos gets moved this will cause a problem. 
          % replace default values with user values, where appropriate. 
fns = fieldnames(par_mineos);
for ii = 1:length(fns)
    parm.(fns{ii}) = par_mineos.(fns{ii});
end

%% filenames
if ~ischar(parm.ID)
    parm.ID = num2str(parm.ID);
end
ID = [parm.ID,parm.R_or_L(1)];
cardfile = [parm.ID,'.card']; % this might be overwritten later, but is default card file name

switch parm.R_or_L(1)
    case 'R'
        modetype = 'S';
    case 'L'
        modetype = 'T';
end

qmod= parm.qmodpath;

%% =======================================================================
wd = pwd;


%% write MINEOS executable and input files format
if ischar(model) && java.io.File([pwd '/' model]).exists; % exist(model,'file')==2 % input model is a card file, not a matlab structure %TODOEXIST bb2021.11.22 exist is SUPER slow
    cardfile = model; % 'model' is just the path to the cardfile
    delcard = false; % do not delete the card file you fed in
else
%% WRITE CARD FILE
    if ~isfield(model,'Sanis')
        model.Sanis = zeros(size(model.z));
    end
    if ~isfield(model,'Panis')
        model.Panis = zeros(size(model.z));
    end
    xi = 1 + model.Sanis/100;  % assumes Sanis is a percentage of anis about zero
    phi = 1 + model.Panis/100;  % assumes Panis is a percentage of anis about zero
    [ vsv,vsh ] = VsvVsh_from_VsXi( model.VS,xi );
    [ vpv,vph ] = VpvVph_from_VpPhi( model.VP,phi );

	write_cardfile(cardfile,model.z,vpv,vsv,model.rho,[],[],vph,vsh);
    delcard = true; % delete this card file afterwards (you still have the model to re-make it if you want)
end

% count lines in cardfile
[ model_info ] = read_cardfile( cardfile );
skiplines = model_info.nlay + 5; % can skip at least this many lines at the beginning of the .asc output file(s)
save([ID,'vel_profile'],'model_info');

% compute max frequency (mHz) - no need to compute past the minimum period desired
parm.fmax = 1000./min(swperiods)+1; % need to go a bit beyond ideal min period..

%% do MINEOS on it
if ifverbose
    fprintf('    > Running MINEOS normal mode summation code. \n    > Will take some time...')
end



%%% Now we have card file. Go to Josh's code to run it. 
% Test 1. Can we simply execute his script? To start. Will be silly. 
disp('Starting to use alternate mineos.'); 
mkdir('CARDS'); 
copyfile(cardfile, sprintf('./CARDS/%s', cardfile)); 
copyfile('/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/MINEOS_synthetics/run_MINEOS/CARDS/prem_35.qmod', sprintf('./CARDS/%s.qmod', parm.ID)); % Temporary


% parameter_FRECHET
% a1_run_mineos_check
% a2_mk_kernels
% a3_pull_dispersion
% %brb20240605 You can add a4 onward for eigenfunctions or synthetic
% %waveforms, but we were not using them. 
% % a4_pull_eigenfuncs
% % b1_run_idagrn
% % % a5_run_idagrn_branch % brb2024/05/22 I don't see this script in the Github repository. Maybe it is no longer needed. 
% % b2_run_idagrn_excitation
% 
% 
% % brb20240605 TODO next order of business is putting Josh's codes in their
% % own functions so we don't clear the workspace. 
% phV = interp1(mode.Tq,mode.phvq,swperiods); % Get the q corrected phase velocities, or I guess the phase velocity that would be measured if there were no q? brb20240605
% grV = interp1(mode.T ,mode.grv ,swperiods); % Group velocity. I don't think this is corrected for Q. 
% % TODO might have to transpose.
% eigfiles_fix = []; % What is this?  
% %%%

% Try putting JOshs stuff in one function 
[mode] = f_mineos_run_phv_grv(); 
% brb20240605 TODO next order of business is putting Josh's codes in their
% own functions so we don't clear the workspace. 
phV = interp1(mode.Tq,mode.phvq,swperiods); % Get the q corrected phase velocities, or I guess the phase velocity that would be measured if there were no q? brb20240605
grV = interp1(mode.T ,mode.grv ,swperiods); % Group velocity. I don't think this is corrected for Q. 
% TODO might have to transpose.
eigfiles_fix = [mode.fname(1:end-1) 'eig']; % Looks like eig files should be the ones that end in .eig. Since, for now, Josh's code has the q file with the same name, just replace q with eig at the end of the q filename. Might not be sustainable.   
%%%

% Below, comment out the original code. 

% % %%% Now we have card file. Go to Josh's code to run it. 
% % % Test 1. Can we simply execute his script? To start. Will be silly. 
% % mkdir('CARDS'); 
% % 
% % copyfile('/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/MINEOS_synthetics/run_MINEOS/CARDS/prem_35.card', sprintf('./CARDS/%s', cardfile)); 
% % copyfile('/Users/brennanbrunsvik/Documents/repositories/Peoples_codes/MINEOS_synthetics/run_MINEOS/CARDS/prem_35.qmod', sprintf('./CARDS/%s.qmod', parm.ID)); % Temporary
% % parameter_FRECHET
% % a1_run_mineos_check
% % %%%
% 
% %% Run mineos once by default
% lrun = 0; lrunstr = num2str(lrun);
% 
% execfile = [ID,'_',lrunstr,'.run_mineos'];
% ascfile =  [ID,'_',lrunstr,'.asc'];
% eigfile =  [ID,'_',lrunstr,'.eig'];
% modefile = [ID,'_',lrunstr,'.mode'];
% 
% writeMINEOSmodefile( modefile, modetype,parm.lmin,parm.lmax,parm.fmin,parm.fmax )
% writeMINEOSexecfile( execfile,cardfile,modefile,eigfile,ascfile,[ID,'.log']);
% 
% fileattrib(execfile, '+x'); 
% [status,cmdout] = system([paths.timeout ' 100 ./',execfile]); % make sure your system can see your timeout function
% [errorInfo]=assess_timeout_results(status, cmdout); 
% %brb2022.04.05 Changing from 100s timeout to 15s. 
% 
% delete(execfile,modefile); % kill files we don't need
% 
% % read prelim output
% [~,llast,lfirst,Tmin] = readMINEOS_ascfile(ascfile,0,skiplines);
%  
% ascfiles = {ascfile};
% eigfiles = {eigfile};
% llasts = llast; lrunstrs = {lrunstr};
% 
% 
% %% Re-start iteratively on higher modes if necessary
% % Tmin = max(swperiods);
% tFact = (24 * 60 * 60); % Factor to convert from serial time to seconds. 
% t_start_rerun = now() * tFact(); 
% max_time_rerun = 100; % Don't run mineos for more than this many more seconds.  Otherwise, we could take 100 s per run 100 times!
% while Tmin > min(swperiods)
% 
%     lrun = lrun + 1; lrunstr = num2str(lrun);
%     lmin = llast + parm.l_increment_standard;
%     modes_done = max(round(llast-lfirst+1)); 
%     
%     if (ifverbose) || (lrun > 10); 
%         fprintf('\n        %4u modes done, failed after mode %u... restarting at %u. Run=%1.0f',modes_done,llast,lmin,lrun)
%     end
%     
%     time_ran = (now() * tFact() - t_start_rerun); 
%     if time_ran > max_time_rerun; % Time since starting to rerun mineos, in seconds.  
%         error('bb2022.07.28 error: Taking %1.2f s to rerun failed Mineos modes. %1.0f modes recalculated.',time_ran, modes_done); 
%     end
%     
%     if lrun > parm.maxrunN
%     %     fprintf('More than %u tries; breaking mineos eig loop\n',parm.maxrunN);
%         error('bb2021.09.21 error: More than %u tries; breaking mineos eig loop\n',parm.maxrunN);
%     end
% 
% 
%     execfile = [ID,'_',lrunstr,'.run_mineos'];
%     ascfile =  [ID,'_',lrunstr,'.asc'];
%     eigfile =  [ID,'_',lrunstr,'.eig'];
%     modefile = [ID,'_',lrunstr,'.mode'];
% 
%     writeMINEOSmodefile( modefile, modetype,lmin,parm.lmax,parm.fmin,parm.fmax )
%     writeMINEOSexecfile( execfile,cardfile,modefile,eigfile,ascfile,[ID,'.log']);
% 
%     fileattrib(execfile, '+x'); 
%     [status,cmdout] = system([paths.timeout ' 100 ./',execfile]); % run execfile 
%     [errorInfo]=assess_timeout_results(status, cmdout); 
%     %brb2022.04.05 Changing from 100s timeout to 15s. 
% 
%     delete(execfile,modefile); % kill files we don't need
% 
%     % read asc output
%     [~,llast,lfirst,Tmin] = readMINEOS_ascfile(ascfile,0,skiplines);
% 
%     if isempty(llast) % what if that run produced nothing?
%         llast=lmin + parm.l_increment_failed;
%         Tmin = max(swperiods);
%         lfirst = llast+1;
%         delete(ascfile,eigfile); % kill files we don't need
%         continue
%     end
% 
%     % save eigfiles and ascfiles for stitching together later
%     ascfiles{length(ascfiles)+1} = ascfile;
%     eigfiles{length(eigfiles)+1} = eigfile;
%     llasts(length(eigfiles)) = llast;
%     lrunstrs{length(eigfiles)} = lrunstr;
% 
% 
% end % on while not reached low period
% if ifplot
%     readMINEOS_ascfile(ascfiles,ifplot,skiplines);
% end
% 
% if (ifverbose) || (lrun > 10)
%      fprintf(' success!\n')
% end
% 
% %% Do eig_file fixing
% eigfiles_fix = eigfiles;
% for ief = 1:length(eigfiles)-1
%     execfile = [ID,'_',lrunstrs{ief},'.eig_recover'];
%     writeMINEOSeig_recover( execfile,eigfiles{ief},llasts(ief) )
%     
%     fileattrib(execfile, '+x'); 
%     [status,cmdout] = system([paths.timeout ' 100 ./',execfile]); % run execfile 
%     [errorInfo]=assess_timeout_results(status, cmdout); 
%     %brb2022.04.05 Changing from 100s timeout to 15s. 
% 
%     eigfiles_fix{ief} = [eigfiles{ief},'_fix'];
%     delete(execfile);
% end
% 
% %% Do Q-correction
% qexecfile = [ID,'.run_mineosq'];
% writeMINEOS_Qexecfile( qexecfile,eigfiles_fix,qmod,[ID,'.q'],[ID,'.log'] )
% 
% fileattrib(qexecfile, '+x'); 
% [status,cmdout] = system([paths.timeout ' 100 ./',qexecfile]); % run qexecfile 
% [errorInfo]=assess_timeout_results(status, cmdout); 
% %brb2022.04.05 Changing from 100s timeout to 15s. 
% 
% 
% delete(qexecfile);
% 
% %% Read phase and group velocities
% try
% %     disp([ID,'.q'])
% %     disp(swperiods)
%     [phV,grV] = readMINEOS_qfile([ID,'.q'],swperiods);
% catch e 
%     warning('Cant read Mineos files. Did timeout not finish?')
%     fprintf('\n%s\n',getReport(e)); 
%     error('some error with extracting phV and grV from q-file')        
% end
%   
% phV = phV(:);
% grV = grV(:);
% if any(isnan(phV)) || any(isnan(grV))         
% 	error('Some NaN data for %s',cardfile)
% end
% 
% %% delete files
% if ifdelete
%     delete([ID,'_*.asc'])
%     delete([ID,'.q'])
%     delete([ID,'_*.eig'])
%     delete([ID,'_*.eig_fix'])
%     if java.io.File([pwd '/' ID '.log']).exists, delete([ID,'.log']); end 
%     if delcard, delete(cardfile); end
% end
% cd(wd);
% 
% %% plot
% if ifplot
%     figure(88), clf; set(gcf,'pos',[331 385 848 613]);
%     set(gcf, 'color', 'white'); 
%     ax1 = axes; hold on;
%     % dispersion curves
%     hd(1)=plot(ax1,swperiods,phV,'o-','linewidth',2);
%     hd(2)=plot(ax1,swperiods,grV,'o-','linewidth',2);
%     hl = legend(ax1,hd,{'Phase (c)','Group (U)'},'location','southeast');
%     set(hl,'fontsize',16,'fontweight','bold')
%     set(ax1,'fontsize',16)
%     xlabel(ax1,'Period (s)','interpreter','latex','fontsize',22)
%     ylabel(ax1,'Velocity (km/s)','interpreter','latex','fontsize',22)
%     box on; 
% end
% 
% if ifverbose
% 	fprintf('Kernels %s took %.5f s\n',ID,(now-tic1)*86400)
% end
 