function [ data ] = z1_SYNTH_DATA_h_kappa(par,trudata,ifplot)

data = trudata; 

% synthetic data parameters
samprate = par.synth.samprate;

for id = 1:length(par.inv.datatypes)
    pdtyps(id,:) = parse_dtype(par.inv.datatypes{id}); 
end

%% ==================== Get H-kappa stack using P wave receiver function ==
if any(string(pdtyps(:,1))=='HKstack'); % Use the Ps receiver function to get h-kappa stack
    phase_wts = [.7, .2, .1]; % Zhu and Kanamori 2000 weights. Ordered as w1, w2, w3, but I'm not sure I Zach used the same ordering in HKstack.m bb2021.12.08
%     [HK_A, HK_H, HK_K] = HKstack(PsRF.PSV(:,2), PsRF.tt, PsRF.rayp/111, phase_wts, ...
%         PsRF.Vs_surf, linspace(25, 70, 200)',linspace(1.6, 2.1, 200)); 
    [HK_A, HK_H, HK_K] = HKstack(trudata.RF_Ps.PSV(:,2), trudata.RF_Ps.tt, ...
        trudata.RF_Ps.rayp/111, phase_wts, ...
        trudata.RF_Ps.Vs_surf, ...
        linspace(10, 70, 200)',linspace(1.5, 2.1, 150)); 
%     warning('10 to 60 km h kappa')

    % The rest of the script expects transposed hk stack. 
    % For IRIS ears data loaded from Zach's functions, K is first
    % dimension, H is second
    % bb2021.12.31
    HK_A = HK_A'; 
    HK_H = HK_H'; 
    HK_K = HK_K'; 
    
    if ifplot; 
        figure(1); clf; hold on; 
        subplot(1,1,1); 
        surf(HK_K, HK_H, HK_A, 'EdgeAlpha', 0)
        xlabel('kappa'); ylabel('H'); title('Synthetic H-kappa stack'); 
        set(gca, 'ydir', 'reverse'); 
        colorbar(); 
    end
    
    % TODO what is Vs_surf? Average vs? Why named surf? 
    % TODO change H vector. Is 25 to 70 in IRIS ears It seems... Might just
    % load a real HK stack and use those same values. 
    % PsRF.inc is P_inc. 
%     S_inc0 = rayp2inc(rayps,TLM.Vs(1),6371-TLM.zlayb(1))/ 
    
%     HKstack_P = struct('H', HK_H, 'K', HK_K', 'Esum', HK_A, ...
%         'Nobs', length(gcarcs), 'dof', length(gcarcs) ); 
    HKstack_P = struct('H', HK_H, 'K', HK_K', 'Esum', HK_A, ...
        'Nobs', par.mod.data.deg_of_freedom.h_kappa, ...
        'dof' , par.mod.data.deg_of_freedom.h_kappa ); 
    warning('bb2021.12.20 artifically setting 7 h-kappa observations')
    % Not sure what dof is. Not sure how to establish Nobs, or if it even
    % matters. 
    data.HKstack_P = HKstack_P; 
    
    data = rmfield(data, 'RF_Ps'); 
    disp('Removing RF_Ps from data, keeping only h-kappa stack')
end



end