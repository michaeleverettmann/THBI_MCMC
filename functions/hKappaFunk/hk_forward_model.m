function [predata, par] = hk_forward_model(...
    par, model, predata, pdtyps, options)
    arguments
        par
        model
        predata
        pdtyps
        options.posterior = []
        options.showPlot = false; % Leave as true... this can only plot on iterations where doing HK stack anyway. 
        options.insistRerun = false; % Turn this to true only if you want full HK stack plot every iteration 
    end

HKdat = par.inv.datatypes{find(strcmp(pdtyps(:,1),'HKstack'),1,'first')};

if (model.vpvs > max(predata.(HKdat).K)) || (model.vpvs < min(predata.(HKdat).K)) ... % Give maximum error if model is outside the bounds of the actual HK stack. 
        || (model.zmoh > max(predata.(HKdat).H)) || (model.zmoh < min(predata.(HKdat).H)); 
    predata.(HKdat).H = model.zmoh;
    predata.(HKdat).K = model.vpvs;
    predata.(HKdat).E_by_Emax = min(min(predata.(HKdat).Esum));
    warning('brb2022.03.02 HK outside bounds! THis should not happen. Thus, I am simply giving maximum HK stack error')
else       
    % TODO brb2022.03.01 temporarily remake HK stack here. 
    % Eventually, start only remaking HK sometimes...

    nSurfRes = 1; % Multiple of surface wave kernel resets that we will also reset HK stack
    if ((par.hkResetInfo.timesReset==0)                || ...
            (par.hkResetInfo.timesSWKernelsReset == nSurfRes) || ...
            (par.ii==par.inv.niter)) || ...
            (options.insistRerun); 
        par.hkResetInfo.timesReset = par.hkResetInfo.timesReset + 1; 
        par.hkResetInfo.timesSWKernelsReset = 0; % Reset our surface wave kernel counter. This only is analyzed for the HK stacks. 

        fprintf('\nResetting HK stack on iteration %1.0f (%1.0fx as often as surface wave kernels)\n', par.ii, nSurfRes);  
        HK_new = zeros(201, 200); % temporary. TODO brb2022.03.01
        for iWave = [1:size(predata.HKstack_P.waves.rf,2)]; 
            RF   = predata.HKstack_P.waves.rf(:,iWave); 
            tt   = predata.HKstack_P.waves.tt; 
            rayp = predata.HKstack_P.waves.rayParmSecDeg(iWave); 
            [HK_A, HK_H, HK_K] = HKstack_anis_wrapper(...
                par, model, RF, tt, rayp, 'ifplot', false, ...
                'hBounds', [par.mod.crust.hmin, par.mod.crust.hmax], ...
                'kBounds', [par.mod.crust.vpvsmin, par.mod.crust.vpvsmax], ...
                'posterior', options.posterior);   
            HK_new = HK_new + HK_A; 
        end
        predata.(HKdat).H    = HK_H; 
        predata.(HKdat).K    = HK_K; 
        predata.(HKdat).Esum = HK_new; 

%             update_HK_plot = mod(par.ii, 50) == 0; 
        if options.showPlot; 
            plot_HK_stack(HK_H, HK_K, HK_new, ...
                'model_vpvs', model.vpvs, 'model_zmoh', model.zmoh, ...
                'title', sprintf('Iteration = %1.0f',par.ii),'figNum', 198) ; 
%             paths = getPaths(); 
%             exportgraphics(gcf, sprintf('%s/HK_%s_at_iter_%07d.pdf',...
%                 par.res.resdir, par.res.chainstr, par.ii) ); 
            exportgraphics(gcf, sprintf('./HK_%s_at_iter_%07d.pdf',...
                par.res.chainstr, par.ii) ) ; 
             % If we want to use HD
        end
    end

    % Extract some values to make things easier to work with. 
    h          = predata.(HKdat).H; 
    k          = predata.(HKdat).K; 
    Esum       = predata.(HKdat).Esum; 
    kTrial     = model.vpvs; 
    hTrial     = model.zmoh; 
    ii         = par.ii; 

    % Energy at models part of h_kappa stack. Interpolate in case of small perterbations to model. 
    E_by_Emax = interpn(k',h,Esum,kTrial,hTrial) ...
                / maxgrid(Esum);


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