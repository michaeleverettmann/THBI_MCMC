function [SW_V_kernels] = run_kernels(swperiods,par_mineos,eigfiles,ifdelete,ifplot,ifverbose,ifanis)
% [SWV_kernels] = run_kernels(swperiods,par_mineos,eigfiles,ifdelete,ifplot,ifverbose,ifanis)
% 
% Function to calculate perturbational phase velocity kernels, having
% previously run MINEOS

paths = getPaths(); 

tic1 = now;

if nargin < 2 || isempty(par_mineos)
    par_mineos = [];
end
if nargin < 3 || isempty(eigfiles)
    eigfiles = {[par_mineos.ID,'_0.eig']}; % brb2021.11.03 TODO This will break if you do not provide par_mineos! Previously, it just utilized the (undefined) variable ID, which is supposed to come from par_mineos. 
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
if nargin < 9 || isempty(ifanis)
    ifanis = false;
end

%% parameters
% This needs to be exactly the same as in run_mineos.m. TODO make a function
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
              'qmodpath',[paths.THBIpath '/matlab_to_mineos/safekeeping/qmod']); % bb2021.09.14 making properly dynamic paths again

fns = fieldnames(par_mineos);
for ii = 1:length(fns)
    parm.(fns{ii}) = par_mineos.(fns{ii});
end

% compute max frequency (mHz) - no need to compute past the minimum period desired
parm.fmax = 1000./min(swperiods)+1; % need to go a bit beyond ideal min period..
% parm.fmin = 1000./max(swperiods)-1; % TODO Should we limit the min frequency?  


% phase or group or both ([1 0] or [0 1] or [1 1] respectively)
ph_gr = [0 0];
if ~isempty(regexp(parm.phV_or_grV,'ph','once')) ||...
    ~isempty(regexp(parm.phV_or_grV,'phase','once')) ||...
    ~isempty(regexp(parm.phV_or_grV,'c','once'))
   ph_gr(1)= true; 
end
if ~isempty(regexp(parm.phV_or_grV,'gr','once')) ||...
    ~isempty(regexp(parm.phV_or_grV,'group','once')) ||...
    ~isempty(regexp(parm.phV_or_grV,'U','once'))
   ph_gr(2)= true; 
end

%% filenames
if ~ischar(parm.ID), parm.ID = num2str(parm.ID);end
ID = [parm.ID,parm.R_or_L(1)];
logfile = [ID,'.log'];
execfile_k = [ID,'.run_kernels'];
stripfile = [ID,'.strip'];
tabfile = [ID,'.table'];
qfile = [ID,'.q'];
kernelfile = [ID,'.frechet'];


% standard inputs, don't get re-written
qmod = parm.qmodpath; % This should resolve any path problems finding qmod. bb2021.09.14

%% =======================================================================
wd = pwd;

%% CALCULATE AND READ IN PERTURBATION KERNELS 
%(frechet derivatves of parm perturbation)

SW_V_kernels = f_mineos_kernels(parm, swperiods); 
