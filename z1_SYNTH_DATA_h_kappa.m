function [ data ] = z1_SYNTH_DATA_h_kappa(par,trudata,trumodel,ifplot)
%% brb2022.02.08 Function to calculate synthetic h-kappa stack from
% synthetic receiver functions. 

ifplot = true; warning('brb2022.02.08 remove ifplot')

samprate = par.synth.samprate;
data = trudata; 

for id = 1:length(par.inv.datatypes)
    pdtyps(id,:) = parse_dtype(par.inv.datatypes{id}); 
end

%% ==================== Get H-kappa stack using P wave receiver function ==
if any(string(pdtyps(:,1))=='HKstack'); % Use the Ps receiver function to get h-kappa stack

    %%% This wrapper takes the receiver function waveforms we have
    % and makes HK stack while accounting for radial anisotropy. 
    waves = struct('rf',   trudata.RF_Ps.PSV(:,2), ...
                   'tt',   trudata.RF_Ps.tt, ...
                   'rayParmSecDeg', trudata.RF_Ps.rayp); % Parameters used in forward modelling HK stack on waveforms. 
    
    trumodelIso = trumodel; warning('!!!')
    trumodelIso.Panis = trumodel.Panis .* 0; 
    trumodelIso.Sanis = trumodel.Sanis .* 0; 
    [HK_A, HK_H, HK_K] = HKstack_anis_wrapper(par, trumodelIso, ...
        waves.rf, waves.tt, waves.rayParmSecDeg,...
        'ifplot', true, ...
        'hBounds', [par.mod.crust.hmin, par.mod.crust.hmax], ...
        'kBounds', [par.mod.crust.vpvsmin, par.mod.crust.vpvsmax] ); 
% 
    [HK_A, HK_H, HK_K] = HKstack_anis_wrapper(par, trumodel, ...
        waves.rf, waves.tt, waves.rayParmSecDeg,...
        'ifplot', true, ...
        'hBounds', [par.mod.crust.hmin, par.mod.crust.hmax], ...
        'kBounds', [par.mod.crust.vpvsmin, par.mod.crust.vpvsmax] ); 
    
    plot_HK_stack(HK_H, HK_K, HK_A, ...
        'model_vpvs', trumodel.vpvs, 'model_zmoh', trumodel.zmoh, ...
        'title', sprintf('Synthetic H-kappa stack (true model values)'),...
        'figNum', 198 ,...
        'saveString', sprintf('%s/HK_stack_synth_truparms.pdf',par.res.resdir) ); 

    % Store stack info in structure. 
    HKstack_P = struct('H', HK_H, 'K', HK_K', 'Esum', HK_A, ...
        'Nobs', par.mod.data.deg_of_freedom.h_kappa, ...
        'dof' , par.mod.data.deg_of_freedom.h_kappa, ...
        'waves', waves, 'RF_Ps', data.RF_Ps); 
    data.HKstack_P = HKstack_P; 
        
    % We needed to calculate synthetic receiver function just to get h-kappa stack. Now, remove the receiver function time series from data. 
    disp('Removing RF_Ps from data, keeping only h-kappa stack.')
    data = rmfield(data, 'RF_Ps'); 
end



end