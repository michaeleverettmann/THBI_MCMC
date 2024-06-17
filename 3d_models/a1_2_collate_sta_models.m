%%% Before running this script, may need to run a1_0_collate_results.m
%%% to get results combined from each station. 

run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 

%% load project & station details
stations = load([infodir 'stations.mat']); 
stainfo = stations.stainfo; 
stainfo.overallQ = ones(size(stainfo.slons)); 

% % % % retrieve comparison models
% % % semPath = '~/Documents/repositories/data/models_seismic/SEMum2_avg_VS'
% % % addpath(semPath); 
% % % a = SEMum2_avgprofiles(0,[semPath '/']);

mdls = struct('nwk', {}, 'sta', {}, 'lat', [], 'lon', [], 'model', {},...
    'misfit', {}, 'fposterior', {}, 'dir', {});  % Fill out structure with models we have. 


%% find stas with good fits overall and good fit to Sp data;
gdstas = zeros(stainfo.nstas,1);
for is = 1:stainfo.nstas 

    clear sdtyp

    %%%
    sta = stainfo.stas{is}; 
    nwk = stainfo.nwk {is}; 
    
    resdir = sprintf('%s%s_%s_dat%.0f/%s',proj.STAinversions,sta,nwk,generation,STAMP);
    fdir = [resdir, '/final_model.mat']; 
    fallmod = [resdir, '/allmodels_perchain_orig.mat']; 
    fgoodchains = [resdir,'/goodchains.mat']; 
    fposterior = [resdir,'/posterior.mat']; 
    misfitdir = sprintf('%s/final_misfit.mat',resdir); 
    gdstas(is) = logical(exist(fdir, 'file'          ) ...
        * logical(exist(misfitdir, 'file')           ) ...
        * logical(exist([resdir,'/par.mat'], 'file') ) ...
        * logical(exist(fallmod, 'file')             ) ...
        * logical(exist(fgoodchains, 'file')         )     ); 
    if ~gdstas(is); continue; end; % can't load par if it doesn't exist. Just do continue. 
       
    par=load([resdir,'/par.mat']); par=par.par; 
    gdstas(is) = gdstas(is) * ((par.inv.niter==desired_iter) && (par.inv.nchains==desired_chains)); % Don't plot results if we didn't run the inversion with enough chains or iterations. Might have been test runs..
    if ~gdstas(is); continue; end; % can't load par if it doesn't exist. Just do continue. 

    if gdstas(is) && ~ exist(fposterior, 'file'); 
        warning('Station %s %s has the other files but not posterior?',nwk,sta) 
        continue; 
    end

    model = load(fdir).final_model; 
    final_misfit = load(misfitdir).final_misfit; 
%     posterior = load(fposterior).posterior; 


    igdsta = sum(gdstas); 

    %%% Save useful things to .mat structure. 
    mdls(1).nwk{igdsta,1} = nwk; 
    mdls.sta{igdsta,1} = sta; 
    mdls.lat(igdsta,1) = stainfo.slats(is); 
    mdls.lon(igdsta,1) = stainfo.slons(is); 
    mdls.model{igdsta,1} = model; 
    mdls.fposterior{igdsta,1} = fposterior; % Takes too much space to duplicate the posterior here. Just put the file position of it. 
    mdls.misfit{igdsta,1} = final_misfit; 
    mdls.dir{igdsta,1} = resdir; 
    %%%


% % %     %%% Trying to make PDF of velocity at each depth. 
% % %     allmod = load(fallmod).allmodels_perchain; 
% % %     goodchains = load(fgoodchains).goodchains; 
% % %     %%%

end


fresults = sprintf('%s/compiled_results_%s.mat',out_dir,STAMP); 
save(fresults,"mdls"); 