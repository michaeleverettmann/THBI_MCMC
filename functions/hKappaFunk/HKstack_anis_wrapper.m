function [HK_A, HK_H, HK_K, t_pred] = HKstack_anis_wrapper(par, model, ...
                                    RF, tt, rayp, options)  
    arguments
        par % Eventually contain phase weights
        model % Model structure
        RF,
        tt,
        rayp, 
        options.ifplot = false; 
        options.phase_wts = [.5, .3, .2]; % brb2022.02.08 1/3 for each gave h-kappa stacks that don't look like others... They had way to tight of constraint on moho depth (for synthetics). I'm making a compromise. 
        options.hBounds = [10, 70]; 
        options.kBounds = [1.5, 2.1]; 
        options.kNum = 201; % kNum ~= hNum to help catch errors. 
        options.hNum = 200; 
        options.posterior = []
    end

% MCMC/THBI wrapper to get relevant model/data and use it for radially 
% anisotropic HK stack. 
%brb2022.11.22 Modifying this to use my hk_anis repository. 

% Extract anisotropy values from the model. 
if isempty(options.posterior); % Have to calculate averages from 1-D profile. 
    xi = model.Sanis/100+1; % Radial S anisotropy. 
    phi = model.Panis/100+1; % Radial P anisotropy
else % If we have the posterior, get average values from there. 
    posterior = options.posterior; 
    xi = median(posterior.xicrust); 
    phi = median(posterior.phicrust); 
    % If there is an error about xi or phi being a scalar in stead of vector, maybe this is where things are going wrong. 
end
eta = ones(size(model.Sanis)); % TODO if introduce non 1 eta, need to change this!!!

% Get average crustal values to be consistent with Zhu and Kanamori math. This also uses earth flattening transformation. 
[rhoav, vsav, vpav, xiav, phiav, etaav] =  hk_average_crust(...
    model.rho, model.VS, model.VP, ...
    xi, phi, eta, ...
    model.z, model.zmoh);  

% Do HK stack on one receiver function. 
[HK_A, HK_H, HK_K, t_pred] = hk_anis(...
    RF, tt, rayp, vsav, rhoav, xiav, phiav, etaav,...
    'kNum', options.kNum, 'hNum', options.hNum,...
    'kBounds', options.kBounds, 'hBounds', options.hBounds,...
    'phase_wts', options.phase_wts); 

% The rest of the script expects transposed hk stack. For IRIS ears data loaded from Zach's functions, K is first dimension, H is second bb2021.12.31
HK_A = HK_A'; 
HK_H = HK_H'; 
HK_K = HK_K'; 

% options.ifplot = true; if options.ifplot; warning('Setting options.ifplot = true'); end; 
if options.ifplot; 
    figure(198); clf; hold on; set(gcf,'color','white');
    subplot(1,1,1); hold on; 
    xlabel('kappa'); ylabel('H'); title('Synthetic H-kappa stack (black dot ignores sediment kappa)'); 
    set(gca, 'ydir', 'reverse');         
    sf = pcolor(HK_K, HK_H, HK_A'); %Can't use surf. It's 3d. It always covers other plotting objects. 
    sf.EdgeAlpha = 0; 
    colorbar(); 
    xlim([min(HK_K), max(HK_K)]); 
    ylim([min(HK_H), max(HK_H)]); 


    % An unfinished attempt to interpolate to find the exact point of maximum h-k energy. 
%         F = griddedInterpolant({HK_K, HK_H},HK_A,'spline');
%         HK_K2 = linspace(min(HK_K)+.1, max(HK_K)-.1, 1000); 
%         HK_H2 = linspace(min(HK_H)+.1, max(HK_H)-.1, 1001); 
%         [HK_K2, HK_H2] = ndgrid(HK_K2, HK_H2); 
%         HK_A2 = F(HK_K2, HK_H2); 
%         Emax = max(max(HK_A2));
%         [xmax,ymax] = find(Emax==HK_A2); 
%         scatter(gca, HK_K2(xmax), HK_H2(ymax), 50, 'red')

    % Simple wasy to estimate maximum h-k energy point. 
    Emax = max(max(HK_A));
    [xmax,ymax] = find(Emax==HK_A); 
    kBest = HK_K(xmax); 
    hBest = HK_H(ymax); 

    % Plot position of max energy. 
    scatter(gca, kBest, hBest, 50, 'red')
    text(kBest, hBest, sprintf('H = %2.2f, k = %1.3f',...
        HK_H(ymax),HK_K(xmax)), 'verticalalignment', 'bottom' )

    % Plot true position of max energy. TODO add to legend. 
    try
        scatter(model.vpvs, model.zmoh, 50, 'k'); 
    catch
        warning('Could not plot model.vpvs and/or model.zmoh'); 
    end

    % Plot the receiver function and the predicted times of different phases. This is a function in HKstack.m 
    HKstack(RF, tt, ...
        rayp/111.1949, phase_wts, vsAv, ... 
        linspace(10, 70, 2000)',linspace(1.5, 2.1, 2001), ...
            'plotPicks',true,'kChoice',kBest,'hChoice',hBest, ...
            'vpvsANIS_vpvs',vpvsEffective/model.vpvs,'Vsv_av',vsvAn); 

end

end