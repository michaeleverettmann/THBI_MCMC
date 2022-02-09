function [ data ] = z1_SYNTH_DATA_h_kappa(par,trudata,trumodel,ifplot)
%% brb2022.02.08 Function to calculate synthetic h-kappa stack from
% synthetic receiver functions. 

% ifplot = true; warning('brb2022.02.08 remoe ifplot')

samprate = par.synth.samprate;
data = trudata; 

for id = 1:length(par.inv.datatypes)
    pdtyps(id,:) = parse_dtype(par.inv.datatypes{id}); 
end

%% ==================== Get H-kappa stack using P wave receiver function ==
if any(string(pdtyps(:,1))=='HKstack'); % Use the Ps receiver function to get h-kappa stack
    
    % Choose weighting for different phases. There are a few valid choices.
%     phase_wts = [.7, .2, .1]; % Zhu and Kanamori 2000 weights. Ordered as w1, w2, w3, but I'm not sure I Zach used the same ordering in HKstack.m bb2021.12.08
%     phase_wts = [1/3, 1/3, 1/3]; % Choose equal weights for each phase. Supposedly this is done in the IRIS-EARS database (Crotwell 2007 dissertation). brb2022.02.08
    phase_wts = [.5, .3, .2]; % brb2022.02.08 1/3 for each gave h-kappa stacks that don't look like others... They had way to tight of constraint on moho depth (for synthetics). I'm making a compromise. 

    % Get average crust velocity. There are again a few choices. 
    vsAv = trumodel.vsAvCrust; % brb2022.02.08 Newly calculating this based on travel time of vertically travelling S wave from Moho to elv=0. See z0_SYNTH_MODEL_splinemoduse.m
%     vsAv = trudata.RF_Ps.Vs_surf; % brb2022.02.08 Not sure how this value was calculated. 

    % Actual h k stack calculation. 
    [HK_A, HK_H, HK_K] = HKstack(trudata.RF_Ps.PSV(:,2), trudata.RF_Ps.tt, ...
        trudata.RF_Ps.rayp/111.1949, phase_wts, ...
        vsAv, ...
        linspace(10, 70, 200)',linspace(1.5, 2.1, 201)); 

    % The rest of the script expects transposed hk stack. For IRIS ears data loaded from Zach's functions, K is first dimension, H is second bb2021.12.31
    HK_A = HK_A'; 
    HK_H = HK_H'; 
    HK_K = HK_K'; 
    
    
    if ifplot; 
        figure(198); clf; hold on; 
        subplot(1,1,1); hold on; 
        xlabel('kappa'); ylabel('H'); title('Synthetic H-kappa stack'); 
        set(gca, 'ydir', 'reverse');         
        sf = pcolor(HK_K, HK_H, HK_A'); %Can't use surf. It's 3d. It always covers other plotting objects. 
        sf.EdgeAlpha = 0; 
        colorbar(); 
        
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
        scatter(trumodel.crustmparm.vpvs, trumodel.zmoh, 50, 'k'); 
        
        % Plot the receiver function and the predicted times of different phases. This is a function in HKstack.m 
        HKstack(trudata.RF_Ps.PSV(:,2), trudata.RF_Ps.tt, ...
            trudata.RF_Ps.rayp/111.1949, phase_wts, vsAv, ... 
            linspace(10, 70, 2000)',linspace(1.5, 2.1, 2001), ...
                'plotPicks',true,'kChoice',kBest,'hChoice',hBest); 
    end
    
    % Store stack info in structure. 
    HKstack_P = struct('H', HK_H, 'K', HK_K', 'Esum', HK_A, ...
        'Nobs', par.mod.data.deg_of_freedom.h_kappa, ...
        'dof' , par.mod.data.deg_of_freedom.h_kappa ); 
    data.HKstack_P = HKstack_P; 
    
    % We needed to calculate synthetic receiver function just to get h-kappa stack. Now, remove the receiver function time series from data. 
    disp('Removing RF_Ps from data, keeping only h-kappa stack.')
    data = rmfield(data, 'RF_Ps'); 
end



end