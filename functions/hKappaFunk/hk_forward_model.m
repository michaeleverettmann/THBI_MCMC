function [predata, parReturn] = hk_forward_model(...
    par, model, predataOrig, pdtyps, predataPrev, options)
    arguments
        par
        model
        predataOrig
        pdtyps
        predataPrev = [] % Will need this any time you don't do full HK recalculation. 
        options.posterior = []
        options.showPlot = false; % Leave as true... this can only plot on iterations where doing HK stack anyway. 
        options.insistRerun = false; % Turn this to true only if you want full HK stack plot every iteration 
    end
    
predata = predataOrig; 
parReturn = par; 

% Par needs some default values that are assigned in the parfor loop but
% are removed after the par for loop. 
if ~isfield(par, 'hkResetInfo'); 
    par.hkResetInfo = struct('timesReset', 0, 'timesSWKernelsReset',0); 
    par.res.chainstr = ''; 
end
if ~isfield(par, 'ii'); 
    par.ii = 0; 
end
% end default par values. 

HKdat = par.inv.datatypes{find(strcmp(pdtyps(:,1),'HKstack'),1,'first')};

if (model.vpvs > max(predata.(HKdat).K)) || (model.vpvs < min(predata.(HKdat).K)) ... % Give maximum error if model is outside the bounds of the actual HK stack. 
        || (model.zmoh > max(predata.(HKdat).H)) || (model.zmoh < min(predata.(HKdat).H)); 
    predata.(HKdat).H = model.zmoh;
    predata.(HKdat).K = model.vpvs;
    predata.(HKdat).E_by_Emax = min(min(predata.(HKdat).Esum));
    warning('brb2022.03.02 HK outside bounds! THis should not happen. Thus, I am simply giving maximum HK stack error')
else % When model is within HK bounds, run this. 
    nSurfRes = 1; % Multiple of surface wave kernel resets that we will also reset HK stack
    
    % Kind of complicated to figure out if we do full HK or just one h and k
    runFullHK = ...
        ((par.hkResetInfo.timesReset==0)                       || ...
        (par.hkResetInfo.timesSWKernelsReset == nSurfRes) || ...
        (par.ii==par.inv.niter))                          || ...
        (options.insistRerun); % Could put this first to avoid the if isfield code, but this will almost certainly cause me to mess up the code in a few months. 

    
    if runFullHK; 
        
        par.hkResetInfo.timesReset = par.hkResetInfo.timesReset + 1; 
        par.hkResetInfo.timesSWKernelsReset = 0; % Reset our surface wave kernel counter. This only is analyzed for the HK stacks. 
        
        fprintf('\nResetting HK stack on iteration %1.0f (%1.0fx as often as surface wave kernels)\n', par.ii, nSurfRes);  
        HK_new = zeros(size(predata.HKstack_P.Esum_orig)); 
        for iWave = [1:size(predata.HKstack_P.waves.rf,2)]; 
            RF   = predata.HKstack_P.waves.rf(:,iWave); 
            tt   = predata.HKstack_P.waves.tt; 
            rayp = predata.HKstack_P.waves.rayParmSecDeg(iWave); 
            [HK_A, HK_H, HK_K] = HKstack_anis_wrapper(...
                par, model, RF, tt, rayp, 'ifplot', false, ...
                'hBounds', [par.mod.crust.hmin, par.mod.crust.hmax], ...
                'kBounds', [par.mod.crust.vpvsmin, par.mod.crust.vpvsmax], ...
                'kNum', size(HK_new,1), 'hNum', size(HK_new,2),...
                'posterior', options.posterior);   
            HK_new = HK_new + HK_A; 
        end
        predata.(HKdat).H    = HK_H; 
        predata.(HKdat).K    = HK_K; 
        predata.(HKdat).Esum = HK_new; 
                
        % Extract some values to make things easier to work with. 
        h          = predata.(HKdat).H; 
        k          = predata.(HKdat).K; 
        Esum       = predata.(HKdat).Esum; 
        kTrial     = model.vpvs; 
        hTrial     = model.zmoh; 
        
        E_by_Emax = interpn(k',h,Esum,kTrial,hTrial) ...
            / maxgrid(Esum);
        
        plot_HK_stack(h, k, Esum, ...
            'model_vpvs', model.vpvs, 'model_zmoh', model.zmoh, ...
            'title', sprintf('Iteration = %1.0f',par.ii),'figNum', 198) 
        
        if options.showPlot; 
            plot_HK_stack(HK_H, HK_K, HK_new, ...
                'model_vpvs', model.vpvs, 'model_zmoh', model.zmoh, ...
                'title', sprintf('Iteration = %1.0f',par.ii),'figNum', 198) ; 

            exportgraphics(gcf, sprintf('./HK_%s_at_iter_%07d.pdf',...
                par.res.chainstr, par.ii) ) ; 
        end
    else % If not remaking whole HK stack, just sample waveforms directly. This is what usually runs. 
        RF   = predata.HKstack_P.waves.rf(:,:); 
        tt   = predata.HKstack_P.waves.tt; 
        rayp = predata.HKstack_P.waves.rayParmSecDeg(:); 
        [E_not_by_Emax] = HK_stack_anis_wrapper_no_grid(...
                par, model, RF, tt, rayp); 
            
        % Extract some values to make things easier to work with. 
        predata.(HKdat).Hgrid = predataPrev.(HKdat).Hgrid; % Store old HK values in new structure for propogating between iterations. The stack is only approximation now, but the specific h-k we sampled is correct.   
        predata.(HKdat).Kgrid = predataPrev.(HKdat).Kgrid; 
        predata.(HKdat).Esum  = predataPrev.(HKdat).Esum; 
        h          = predataPrev.(HKdat).Hgrid; 
        k          = predataPrev.(HKdat).Kgrid; 
        Esum       = predataPrev.(HKdat).Esum; 
        kTrial     = model.vpvs; 
        hTrial     = model.zmoh; 
        
        E_by_Emax = E_not_by_Emax / maxgrid(Esum); % Dividing over maxgrid of the old HK stack. This induces a small amount of error, only in so much as the maximum value in the HK stack changes between iterations. I think this is small... We recalculate the full HK stack every so often. 
        
        E_by_Emax = min(E_by_Emax, 1); % brb2022.03.06 Just in case an un-updated Esum in maxgrid causes us to have higher than 1 normalized energy. That could cause problems when converting energy to mismatch (negative mismatch). 
        
%         plot_HK_stack(h, k, Esum, ...
%             'model_vpvs', model.vpvs, 'model_zmoh', model.zmoh, ...
%             'title', sprintf('Iteration = %1.0f',par.ii),'figNum', 198) 
        
        displayHKGridChange = false; 
        if displayHKGridChange; % If we want to see how much the value of energy at h,k has changed from the old grid versus the value sampled at only h and k just now. 
            E_by_Emax2 = interpn(k',h,Esum,kTrial,hTrial) ...
                / maxgrid(Esum)
            fprintf('\nHK change: %1.4f%%\n',(E_by_Emax-E_by_Emax2)./E_by_Emax * 100)
        end

    end

    % brb2022.03.04 The distance approach below is mostly irrelevant. 
    % Right now I always set all weight toward the actual HK stack values. 
    
    % Find the portion of the h_kappa stack that gives maximum energy. 
    [rowMin, colMin] = find(max(max(Esum))==Esum); 
    kBest            = k(rowMin); % TODO verify that these are in correct order. Find might return collumn then row. 
    hBest            = h(colMin); 

    % Get a "Euclidian" distance for h and kappa. 
    % Will scale to penalize h-kappa with penalty of 1 for a scaled
    % radius of 1 from the best h and kappa value. 
    kDistScale = 0.5 ; % Can do this: max(k) - min(k). However, this means that any time priors change, the entire inversion will work differently. This value is artificial anyway. I used 0.5, as that is what is in use as of 2022.01.18. 
    hDistScale = 45  ; % max(h) - min(h)
    EDist = sqrt(   ((kBest - kTrial) / kDistScale) ^ 2 + ...
                    ((hBest - hTrial) / hDistScale) ^ 2       ); 

    % Convert EDist to a penalty. 1 is perfect fit (0 distance). 0 is worst fit, or max
    % distance. It's ok if it becomes negative, but that will be higher
    % weight than E_by_Emax will generally give. 
    EDist = 1 - EDist; 

% % The commented code here is an alternative way to weight between E_by_Emax and Edist.                 
% %         % Now give weighted average of EDist and E_by_Emax
% %         % Fraction of the way through burnin. 
% %         % Will use to weight E, so don't let go beneath 0. 
% %         % When value is 0 use only E_by_Emax
% %         fracBurn = max( ...
% %             (par.inv.burnin - par.ii)/par.inv.burnin, ...
% %             0); 
% %         
% %         % Relative weighting for EDist and E_by_Emax
% %         % When value is 0 use only E_by_Emax
% %         % When value is 1 use only EDist
% %         wDistEmax = 0.5; 
% %         wDist = (    fracBurn) * (    wDistEmax); 
% %         wEmax = (1 - fracBurn) * (1 - wDistEmax); 
% %         
% %         % Scale combination of weights to make them sum to 1. 
% %         % This should help prevent h-kappa stacks sigma (error) 
% %         % value from being inverted as totally wrong during burnin. 
% %         wTot = wDist + wEmax; 
% %         wDist = wDist / wTot; 
% %         wEmax = wEmax / wTot; 

% %         fracBurn = max( ...
% %             (par.inv.burnin - par.ii)/par.inv.burnin, ...
% %             0); 

% %         % Assign the penalty value. 
% %         predata.(HKdat).E_by_Emax = E_by_Emax * wEmax + EDist * wDist; 

    % Give relative weighting to EDist
    % Start with max val. 
    % Scale toward min val by end of burnin. 
    % NOTE I'm not using this anymore, so use 0 for dist and 1 for stack. brb2022.02.08
    wDistMax = par.datprocess.HKappa.weightDistanceMax; 
    wDistMin = par.datprocess.HKappa.weightDistanceMin; 
    wDistWithSlope = wDistMax - par.ii * (wDistMax-wDistMin)/par.inv.burnin; % Intercept slope formula with par.ii is independent variable
    wDist = max(wDistMin, wDistWithSlope); 
%         sprintf('wDist = %1.2f', wDist)
    wEmax = 1 - wDist; 

    % Assign the penalty value. 
    predata.(HKdat).E_by_Emax = E_by_Emax * wEmax + EDist * wDist; 

    % Assign other value for compatibility with remaining code. 
    predata.(HKdat).Hgrid = predata.(HKdat).H;
    predata.(HKdat).Kgrid = predata.(HKdat).K; % For some reason these are getting replaced with single values, making it hard to plot later. So keep track of the grid values here. 
    predata.(HKdat).H = model.zmoh;
    predata.(HKdat).K = model.vpvs;        
end


end