
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 

%% Parameters to set. 
nkernel = 100; % Number of points in histogram kernels. 
widthkernel = 0.02; % Found this looks ok in t1_make_pdf.m

%% Loading
fresults = sprintf('%s/compiled_results_%s.mat',out_dir,STAMP); 
fpdfs    = sprintf('%scompiled_pdfs_%s.mat',out_dir,STAMP); 
mdls = load(fresults).mdls; 

%% Initiate pdf structures. 
pdfs = struct('nwk', {}, 'sta', {}, 'mm', {}, 'pm', {}); 
pdfs_empty = pdfs; 

%%% Fill structures, pdfs for individual parameters. 
indiv_parameters = ["zsed","zmoh","kcrust","kmantle","VSsedtop","VSsedbot",...
    "VScrusttop","VScrustbot","VSmanttop","fdVSsed",...
    "fdVSmoh","vpvs","xicrust","ximant"]'; % Individual parameters that are m-models x 1 vectors. Same math applies to each. 

pdfs_allparm = struct(); 
for iparam=1:length(indiv_parameters); 
    fn = indiv_parameters(iparam); % field name. 
    pdfs_allparm(1).(fn) = pdfs_empty; 
end

%%% Fill structures, pdfs for different depths velocities. 
posterior_temp = load(mdls.fposterior{1}).posterior; 
depth_not_changed = posterior_temp.zatdep' == ([5 10 15	20	25	30	35	40	45	50	55	60	65	70	75	80	85	90	95	100	105	110	115	120	125	130	135	140	145	150	155	160	165	170	175	180	185	190	195	200	205	210	215	220	225	230	235	240	245	250	255	260	265	270	275	280	285	290	295	300]); 
if ~ all(depth_not_changed); 
    error('Double check that each posterior has the same zatdepth for getting pdfs of that depths velocity. As of 2022.10.11 it was the long series of numbers typed above. ')
end
zatdep = posterior_temp.zatdep; 
pdfs_allparm.vs = {}; 
for iz = 1:length(zatdep); 
    pdfs_allparm(1).vs{iz} = pdfs_empty; 
end
% % % pdfs_allparm.vs = cell2struct(pdfs_allparm.vs, fieldnames(pdfs_empty)', length(zatdep) )'; 
% % % pdfs_allparm.zatdep = zatdep; 

%% Loop over posteriors and make pdfs. 
for is = 1:length(mdls.lon); 

    posterior = load(mdls.fposterior{is}).posterior;
%     pdfs(is).nwk = mdls.nwk{is}; 
%     pdfs(is).sta = mdls.sta{is}; 

    % Loop over most individual parameters. 
    for iparam = 1:length(indiv_parameters); 
        fn = indiv_parameters(iparam); 
        samps = posterior.(fn); 

        [pdfm, mm] = ksdensity(samps, 'width', widthkernel, 'NumPoints', nkernel); % pdf of model parameter. And ... model parameters. 
%         figure(1); clf; hold on; 
%         plot(mm, pdfm); 
        pdfs_allparm(is).(fn)(1).mm = mm'; 
        pdfs_allparm(is).(fn)(1).pm = pdfm'; 
        pdfs_allparm(is).(fn)(1).nwk = mdls.nwk{is}; 
        pdfs_allparm(is).(fn)(1).sta = mdls.sta{is}; 
    end

    % Look over each depth. DIfferent structure organization than for other parameters. 
    pdfs_allparm(is).zatdep = zatdep; 
    for iz = 1:length(zatdep); 
        depth = posterior.zatdep(iz); 
        samps = posterior.VSmantle(:,iz); 
        [pdfm, mm] = ksdensity(samps, 'width', widthkernel, 'NumPoints', nkernel); % pdf of model parameter. And ... model parameters. 
%         figure(1); clf; hold on; 
%         plot(mm, pdfm); 
        pdfs_allparm(is).vs{iz}(1).mm = mm'; 
        pdfs_allparm(is).vs{iz}(1).pm = pdfm'; 
        pdfs_allparm(is).vs{iz}(1).nwk = mdls.nwk{is}; 
        pdfs_allparm(is).vs{iz}(1).sta = mdls.sta{is}; 
    end

    fprintf('%1.2f%% Done\n', is/length(mdls.lon)*100 )

end

% save(fpdfs, 'pdfs', 'is', 'iz'); % Save iz and is as a reminder that this might not be complete.  
save(fpdfs, 'pdfs_allparm', 'is', 'iz'); 
