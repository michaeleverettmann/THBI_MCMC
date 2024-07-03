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
% 
% Modified from old matlab_to_mineos folder to work with Josh's Mineos. 

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

% compute max frequency (mHz) - no need to compute past the minimum period desired
parm.fmax = 1000./min(swperiods)+1; % need to go a bit beyond ideal min period..
% parm.fmin = 1000./max(swperiods)-1; % TODO Should we limit the min frequency?  


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

% Try putting Joshs stuff in one function 
[mode] = f_mineos_run_phv_grv(parm, swperiods); 
phV = interp1(mode.Tq,mode.phvq,swperiods); % Get the q corrected phase velocities, or I guess the phase velocity that would be measured if there were no q? brb20240605
grV = interp1(mode.T ,mode.grv ,swperiods); % Group velocity. I don't think this is corrected for Q. 
eigfiles_fix = [mode.fname(1:end-1) 'eig']; % Looks like eig files should be the ones that end in .eig. Since, for now, Josh's code has the q file with the same name, just replace q with eig at the end of the q filename. Might not be sustainable.   
% figure(1); clf; hold on; plot(mode.Tq, mode.phvq); 

if any(isnan(phV)); 
    warning('Nan phase velocities')
end 