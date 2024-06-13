function [misfits,allmodels,goodchains,misfits_orig,allmodels_orig,allmodels_collated]...
               = c1_PROCESS_RESULTS( misfits,allmodels,par,ifsave,ofile, goodChainsManual )
% [misfits,allmodels,goodchains] = c1_PROCESS_RESULTS( misfits,allmodels,par,ifsave,ofile )
% 
% Script to process the results and make some plots of the misfit + the
% likelihood changing with iteration

if nargin < 4 || isempty(ifsave)
    ifsave=false;
end
if nargin < 5 || isempty(ofile)
    ofile = 'figs/misfits_vs_iter';
end
if nargin < 6 || isempty(goodChainsManual)
    goodChainsManual = []; 
end

figure(88);clf, set(gcf,'pos',[850 198 900 888])
% title
htit = title_custom([par.data.stadeets.sta,' ',par.data.stadeets.nwk],0.5,'fontweight','bold','fontsize',25);

ax1 = subplot(2,1,1); ax1pos = get(ax1,'pos');
ax3 = subplot(2,1,2); 
ax2 = axes('pos',ax1pos); 

set(ax1,'Color','none');
set(ax2,'color','none','YAxisLocation','right'); 
hold(ax1,'on');hold(ax2,'on');hold(ax3,'on');

pdm_min = inf;
pdm_max = -inf;
ch2_min = inf;
ch2_max = -inf;

downsampfac = 2;

%% loop through each chain
nchains = length(misfits);
goodchains = true(nchains,1);

for iii = 1:nchains
    
%     if nchains>1, mf = misfits{iii};   else mf = misfits; end
%     if nchains>1, am = allmodels{iii}; else am = allmodels; end
    mf = misfits{iii};
    am = allmodels{iii};
    if isempty(am), continue; end
    try am(1).Nstored; catch, am = am{1}; end
    basecol = colour_get(iii,nchains+1,0,parula); basecol = basecol(:)';
    
    %% TRIM MISFITS AND ALLMODELS STRUCTURES

    N = am(1).Nstored;

    % kill allmodels beyond N
    am(N+1:end) = [];

    % kill misfits beyond N
    fns = fieldnames(mf);
    for ifn = 1:length(fns)
        if length(mf.(fns{ifn})) > N
            mf.(fns{ifn})(N+1:end) = [];
        end
    end
    
    %% move on if not enough good results
    if length(am) < par.inv.bestNmod2keep && ~isinf(par.inv.bestNmod2keep)
        mf.bestmods = false(mf.Nstored,1);
        am = dealto(am,'bestmods',false);
        goodchains(iii) = false;
        if nchains>1, misfits{iii} = mf;   else misfits = mf; end
        if nchains>1, allmodels{iii} = am; else allmodels = am; end
        continue
    end
    
    %% identify "best" models 
%     from those in the top bestNmod2keep based on summed ranking of their fits to each datatype

    bestind = [1:length(mf.iter)]';
    bestind(mf.iter <= par.inv.burnin) = []; % kill all before burnin
    
    %%% Remove models that had particularly low liklihood. These could come up if resetting kernels causes a huge drop in likelihood. Bad models might be accepted around those iterations, so just ignore them. 
    % Currently done on per-chain basis. Might be better to do across all chains combined, but we would risk eliminating chains who sampled local minima (even if the local minima is the true earth model). 
    logL_dtype = nan(length(bestind),length(par.inv.datatypes));
    for id = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{id};
        logL = [mf.logL_indivdat.(dtype)]'; % all logL, including outside of bestind. 
        logL_dtype(:,id) = logL(bestind); 
    end
    logL_mean = mean(logL_dtype, 1); % Mean of logl for each data type. 
    logL_std = std(logL_dtype,1); 
    cutoff_logL = logL_mean - 3*logL_std; % Cutoff values, too low of logL to keep. 
    ignore_logL = any(logL_dtype < cutoff_logL, 2); % For each iteration, remove it if any of the data types had log L below cutoff. 
    bestind(ignore_logL) = []; 
    perc_removed = 100 * sum(ignore_logL) / length(ignore_logL); 
    fprintf(['\nRemoved %1.2f%% of models from chain %1.0f due to logL'...
        '< mean - 3(?) std of logL for some datype.\n'], perc_removed, iii); 
    %%% END remove models that had particularly low likelihood 

    if ~isinf(par.inv.bestNmod2keep) % subset if not inf to keep. Else keep all
        if par.inv.bestNmod2keep>0 % if specifying how many to keep based on low error
    
            score_overall = zeros(length(bestind),length(par.inv.datatypes));
            for id = 1:length(par.inv.datatypes)
                dtype = par.inv.datatypes{id};
                chi2 = [mf.chi2.(dtype)]';
                [~,irank_mf] = sort(sum(chi2(bestind,:),2));
                [~,score_overall(:,id)] = sort(irank_mf);
            end
            score_overall = sum(score_overall,2);   
            sort_score_overall = sort(score_overall);
            min_score_overall = sort_score_overall(min([par.inv.bestNmod2keep,length(bestind)]));
            bestind(score_overall>min_score_overall) = [];
        elseif par.inv.bestNmod2keep<0
            bestind = bestind(randperm(length(bestind),min([-par.inv.bestNmod2keep,length(bestind)])));
        end
    end
        
    mf.bestmods = false(mf.Nstored,1);
    mf.bestmods(bestind) = true;
    
%     subplot(211),semilogy(mf.iter,mf.chi2,'bo',mf.iter(mf.bestmods),mf.chi2(mf.bestmods),'ro')
%     bestind = [1:mf.Nstored]';
%     bestind(misfits.iter <= par.inv.burnin) = [];
%     [~,rank_mf] = sort(mf.chi2(bestind));
%     max_rank_ind = min([par.inv.bestNmod2keep,length(bestind)]);
%     bestind = bestind(rank_mf(1:max_rank_ind));
%     mf.bestmods = false(mf.Nstored,1);
%     mf.bestmods(bestind) = true;
%     subplot(212),semilogy(mf.iter,mf.chi2,'bo',mf.iter(mf.bestmods),mf.chi2(mf.bestmods),'ro')

    am = dealto(am,'bestmods',mf.bestmods);

    %% plug back into structures
    if nchains>1, misfits{iii} = mf;   else misfits = mf; end
    if nchains>1, allmodels{iii} = am; else allmodels = am; end

    
    %% PLOT IMPROVEMENT IN MODEL FIT
    downsamp = (1:downsampfac:length(am));

    h1 = plot(ax1,mf.iter(downsamp),mf.chi2sum(downsamp));
    h2 = plot(ax2,mf.iter(downsamp),mf.logLike(downsamp));
    
    pdm_min = min([pdm_min,min(mf.logLike)-1]);
    pdm_max = max([pdm_max,max(mf.logLike)+10]);
    ch2_min = min([ch2_min,min(mf.chi2sum)/2]);
    ch2_max = max([ch2_max,max(mf.chi2sum)*2]);

    col1 = [0 0.447 0.741];
    col2 = [0.85 0.325 0.098];

    ytk1 = unique(round_level(linspace(ch2_min,ch2_max,5),5));
    set(ax1,'Yscale','log','ylim',[ch2_min ch2_max],'ytick',ytk1,'Fontsize',16,'ycolor',col1);
    ytk2 = unique(round_level(linspace(pdm_min+1,pdm_max-10,5),5));
    set(ax2,'Yscale','log','ylim',[pdm_min pdm_max],'ytick',ytk2,'yticklabel',ytk2,'xtick',[],'Fontsize',16,'ycolor',col2);
    
    set(h1,'marker','o','linestyle','none','markersize',5,...
           'markerfacecolor',basecol,'markeredgecolor',col1);
    set(h2,'marker','o','linestyle','none','markersize',5,...
           'markerfacecolor',basecol,'markeredgecolor',col2);
    set(get(ax1,'ylabel'),'String','$\chi^2$ misfit','Fontsize',20,'interpreter','latex','color',col1)
    set(get(ax2,'ylabel'),'String','$\log_{10}{\,p(m|d)}$','Fontsize',20,'interpreter','latex','color',col2)
        
    Nd = length(par.inv.datatypes);
    cols = colour_get([1:Nd],Nd,1,[parula;flipud(spring)]);
    
    for idt = 1: Nd
        dtype = par.inv.datatypes{idt};
        rms = [mf.rms.(dtype)]';
        hrms(idt)=plot(ax3,mf.iter(downsamp),sum(rms(downsamp,:),2),'o'); hold on
        set(hrms(idt),'marker','o','linestyle','none','markersize',5,...
           'markerfacecolor',basecol,'markeredgecolor',cols(idt,:));
    end
    hrms(idt+1) = hrms(1); set(hrms(idt+1),'markerfacecolor','none');
    set(hrms(2:idt+1),'linewidth',2)
    set(ax3,'yscale','log','Fontsize',16 )
    ylabel(ax3,'RMS misfit','Fontsize',20,'interpreter','latex')
    xlabel(ax3,'Iteration','Fontsize',20)
    hl = legend(hrms(1:idt),strrep(par.inv.datatypes,'_','-'));
end % loop on chains

if ifsave
    fprintf('saving, may take a while\n')
    save2png(88,ofile,'/');
end
pause(0.05)

%% PARSE GOOD AND BAD CHAINS
chi2_alldata = nan(nchains,length(par.inv.datatypes));

% gather average rms errors of each chain
nchainsRej = 4;
nIterRej = 500; 
if (nchains <= nchainsRej) || (par.inv.niter < nIterRej); % bb2021.09.17 If you have too few chains, they often fail here. This makes debugging hard. So, if there are too few models, let's assume it is for debugging and not reject them.  
    warning('BRENNAN WARNING: nchains is less than %1.0f or nIter is less than %1.0f. Thus, NOT rejecting chains. Only use this few of chains if you are debugging.', nchainsRej, nIterRej)
    goodchains = 1:nchains;
else
    for iii=1:nchains
        if goodchains(iii)==false, chi2_alldata(iii,:) = nan; continue; end % already know it's bad
        for id = 1:length(par.inv.datatypes)
            % assign structures
            if nchains>1, mf = misfits{iii}; else mf = misfits; end
            if isempty(mf), goodchains(iii)=false; continue; end
            ind = mf.iter > par.inv.burnin & mf.iter~=0;
            dtype = par.inv.datatypes{id};
            % work out chi2 for this dtype (average across all data streams for this dtype)
            chi2 = [mf.chi2.(dtype)]';
            if any(isnan(chi2(ind))); 
                fprintf('Nan in chi2, %1.0f cases. Should not happen... Ignoring chi2 for this datatype for this iteration when determining whether this chain was good overall. dtype = %s. \n', sum(isnan(chi2)), dtype); 
            end
            chi2_alldata(iii,id) = nanmean(sum(chi2(ind,:),2)); 
            % work out if the chain got stuck - if there is no change to the
            % data over many iterations - must be stuck for 600 iterations to
            % signify
            Nstuck = 600;
            if any(any(diff(chi2(ind,:),ceil(Nstuck./par.inv.saveperN),1)==0)) % if differences between chi2 at index i and i+Nstuck is 0...
                if strcmp(dtype,'HKstack_P'), continue; end % don't do for HK stack - may stick for ages!
                fprintf('Chain %s stuck\n',mkchainstr(iii));
                chi2_alldata(iii,id) = nan; 
            end
        end
    end

    % goodchains = true(nchains,1);
    fprintf(['\nStarting with %1.0f chains. ',...
        'Removing some if high chi2 for any dtype'],...
        sum(goodchains)); 
    for id = 1:length(par.inv.datatypes) % brb2022.07.18. This looks across each data type. It finds the mean chi2 across each chain for that datatype. Then it finds the std of chi2 for that datatype (for some reason only considering the chains that had lower than average chi2...?). Then for chains with chi2 greater than the mean + 5 * the standard deviation, for any data type(?), that chain is removed. 
        mean_chi2_dtp = nanmean(chi2_alldata(:,id));
        std_chi2_gdtp = nanstd(chi2_alldata(chi2_alldata(:,id)<mean_chi2_dtp,id)); % bb2021.09.17 Might get debugging problems here if nanstd(f) calculates on f where f has only 1 non nan value. ::: Why are we calculating the standard deviation only for chains that did better than average? That's not exactly a standard deviation...
        keep_this_d = (chi2_alldata(:,id) < mean_chi2_dtp + 5*std_chi2_gdtp); % bb2021.09.17 chains might fail here if you are using small iteration chains for debugging. 
%         if any(~keep_this_d); 
%             fprintf(['\nRemoving chain %1.0f for %s, chi2=%1.3f, ',...
%                 'mean chi2 across chains=%1.3f, std=%1.3f\n\n'],...
%                 find(~keep_this_d),par.inv.datatypes{id},...
%                 chi2_alldata(~keep_this_d,id),...
%                 mean_chi2_dtp,std_chi2_gdtp ); 
%         end
        goodchains = goodchains & keep_this_d; 
    end
    fprintf(['\n -> Ended with %1.0f chains. \n'],sum(goodchains) ); 
    goodchains=find(goodchains);
end

%% Option to specify good chains 
if ~isempty(goodChainsManual)
    goodchains = goodChainsManual; 
end

%% GET TARGET MODEL for comparison
global TRUEmodel
if ~isempty(TRUEmodel)
    Z = TRUEmodel.Z;
    vs = TRUEmodel.vs;
    vp = TRUEmodel.vp;
    rho = TRUEmodel.rho;
    fprintf('TRUE sed thickness = %.1f km\n',Z(find(vs>=3.2,1,'first')))
    fprintf('TRUE moho depth = %.1f km\n',Z(find(vs>4.0,1,'first')))
    fprintf('TRUE Vs seds top = %.1f km/s\n',vs(1))
    fprintf('TRUE Vs seds bot = %.1f km/s\n',vs(find(vs<3.45,1,'last')))
    fprintf('TRUE Vs crust top = %.1f km/s\n',vs(find(vs>3.45,1,'first')))
    fprintf('TRUE Vs crust bot = %.1f km/s\n',vs(find(rho<3.2,1,'last')))
    fprintf('TRUE fractional dVs sed/crust = %.1f %% \n',-100*(vs(find(vs<3.5,1,'last')) - vs(find(vs>3.5,1,'first')))/vs(find(vs<3.5,1,'last')))
    fprintf('TRUE fractional dVs crust/mantle = %.1f %% \n',-100*(vs(find(rho<3.2,1,'last')) - vs(find(rho>3.2,1,'first')))/vs(find(rho<3.2,1,'last')))
    for ii = linspace(par.mod.sed.hmax+par.mod.crust.hmax,par.mod.maxz,6)
        try fprintf('TRUE Vs at %.0f km = %.2f km/s\n',ii,linterp(Z,vs,ii));end
    end
end

%% collate all
misfits_orig = misfits;    % misfits_perchain = misfits_perchain_original;
allmodels_orig = allmodels;% allmodels_perchain = allmodels_perchain_original;
if par.inv.nchains > 1
    allmodels = allmodels(goodchains);
end
misfits = misfits(goodchains);

[ allmodels_collated ] = collate_allmodels_perchain( allmodels,par );