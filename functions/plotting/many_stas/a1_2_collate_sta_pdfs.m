%%% brb2022.09.27. In progress. Might come back to this. But the
%%% posterior.mat files might be all I need, so I'll try those first. 

run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 

%%
nkernel = 100; % Number of points in histogram kernels. 
widthkernel = 0.02; % FOund this looks ok in t1_make_pdf.m

fresults = sprintf('%s/compiled_results_%s.mat',out_dir,STAMP); 
fpdfs    = sprintf('%scompiled_pdfs_%s.mat',out_dir,STAMP); 
mdls = load(fresults).mdls; 

pdfs = struct('nwk', {}, 'sta', {}, 'mm', {}, 'pm', {}); 

for is = 1:length(mdls.lon); 
% for is = 1:20; 
    posterior = load(mdls.fposterior{is}).posterior;
    pdfs(is).nwk = mdls.nwk{is}; 
    pdfs(is).sta = mdls.sta{is}; 
    for iz = 40; 
        depth = posterior.zatdep(iz); 
        samps = posterior.VSmantle(:,iz); 
        [pdfm, mm] = ksdensity(samps, 'width', widthkernel, 'NumPoints', nkernel); % pdf of model parameter. And ... model parameters. 
%         figure(1); clf; hold on; 
%         plot(mm, pdfm); 
        pdfs(is).mm = mm'; 
        pdfs(is).pm = pdfm'; 
        
    end
end

save(fpdfs, 'pdfs', 'is', 'iz'); % Save iz and is as a reminder that this might not be complete.  

% %% find stas with good fits overall and good fit to Sp data;
% gdstas = zeros(stainfo.nstas,1);
% for is = 1:stainfo.nstas 
% 
%     clear sdtyp
% 
%     %%%
%     sta = stainfo.stas{is}; 
%     nwk = stainfo.nwk {is}; 
%     
%     resdir = sprintf('%s%s_%s_dat%.0f/%s',proj.STAinversions,sta,nwk,generation,STAMP);
%     fdir = [resdir, '/final_model.mat']; 
%     fallmod = [resdir, '/allmodels_perchain_orig.mat']; 
%     fgoodchains = [resdir,'/goodchains.mat']; 
%     misfitdir = sprintf('%s/final_misfit.mat',resdir); 
%     gdstas(is) = logical(exist(fdir, 'file'          ) ...
%         * logical(exist(misfitdir, 'file')           ) ...
%         * logical(exist([resdir,'/par.mat'], 'file') ) ...
%         * logical(exist(fallmod, 'file')             ) ...
%         * logical(exist(fgoodchains, 'file')         )     );  
%     if ~gdstas(is); continue; end; % can't load par if it doesn't exist. Just do continue. 
%     par=load([resdir,'/par.mat']); par=par.par; 
%     gdstas(is) = gdstas(is) * ((par.inv.niter==desired_iter) && (par.inv.nchains==desired_chains)); % Don't plot results if we didn't run the inversion with enough chains or iterations. Might have been test runs..
%     if ~gdstas(is); continue; end; % can't load par if it doesn't exist. Just do continue. 
% 
%     model = load(fdir).final_model; 
%     final_misfit = load(misfitdir).final_misfit; 
% 
%     igdsta = sum(gdstas); 
% 
%     %%% Save useful things to .mat structure. 
%     mdls(1).nwk{igdsta,1} = nwk; 
%     mdls.sta{igdsta,1} = sta; 
%     mdls.lat(igdsta,1) = stainfo.slats(is); 
%     mdls.lon(igdsta,1) = stainfo.slons(is); 
%     mdls.model{igdsta,1} = model; 
%     mdls.misfit{igdsta,1} = final_misfit; 
%     mdls.dir{igdsta,1} = resdir; 
%     %%%
% 
% 
%     %%% Trying to make PDF of velocity at each depth. 
%     allmod = load(fallmod).allmodels_perchain; 
%     goodchains = load(fgoodchains).goodchains; 
%     %%%
% 
% end


