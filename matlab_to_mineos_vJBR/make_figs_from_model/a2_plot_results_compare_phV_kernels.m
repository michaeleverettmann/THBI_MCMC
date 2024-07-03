clear; clc; close all; 

% versions = ["old", "new"]; 

% Where we are storing figs, and pre data and kernels.
root_dir = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/matlab_to_mineos_vJBR/make_figs_from_model/'; 
fdat = [root_dir 'fig_out_%s_mineos/%s'];
fig_comp = [root_dir, 'fig_comparison/']; 

chdir(root_dir); 

datn = load(sprintf(fdat, 'new', 'dat_and_kbase.mat')); 
dato = load(sprintf(fdat, 'old', 'dat_and_kbase.mat')); 

types = ["Ray", "Lov"]; 

%% Plot phase velocity difference. 
figure(1); clf; hold on; 
for itype = 1:2; 
    box on; 

    kn = datn.kbase.(types(itype)); 
    ko = dato.kbase.(types(itype)); 

    yyaxis('left'); 
    ylim([3, 4.5]); 
    xlabel('Period (s)'); 
    ylabel('Phase velocity (km/s)'); 
    plot(ko.periods, ko.phV, '-k', 'DisplayName', 'Zach'); 
    plot(kn.periods, kn.phV, '-b', 'DisplayName', 'Josh'); 

    yyaxis('right'); 
    set(gca(), 'YColor', 'r'); 
    ylim([-.05, 0.05]); 
    ylabel('Difference (km/s)'); 
    plot(ko.periods, kn.phV-ko.phV, '-r', 'DisplayName', 'Difference'); 
    
    if itype == 1; 
        legend(); 
    end
end
exportgraphics(figure(1), [fig_comp, 'phV_.pdf']); 


%% Plot kernel differences. 
 
for itype = 1:2; 
    figure(2); clf; hold on; 
    tiledlayout(3,2, 'TileSpacing','tight');
    sgtitle(types(itype));

    kn = datn.kbase.(types(itype)); 
    ko = dato.kbase.(types(itype)); 

    parms = ["Vsv", "Vsh", "Vpv", "Vph", "eta", "rho"]; 
    for iparm = 1:length(parms); 
        parm = parms(iparm); 
        nexttile(iparm); hold on; 
        box on; 
        ylabel(parm); 
        xlim([0, 250]); 
        for iper = round(linspace(1, length(ko.phV), 4)); % 
            plot(ko.Kph{iper}.Z/1000, ko.Kph{iper}.(parm), '-k', 'DisplayName', 'Zach', 'linewidth', 1); 
            plot(kn.Kph{iper}.Z/1000, kn.Kph{iper}.(parm), '--b', 'DisplayName', 'Josh', 'linewidth', 1.5); 
            % plot(ko.Kph{iper}.(parm), kn.Kph{iper}.(parm)); 
            % plot(ko.Kph{iper}.Z, kn.Kph{iper}.Z); 
        end
    end
    legend(); 
    xlabel('Depth (km)'); 
    exportgraphics(figure(2), [fig_comp, 'kernels_', convertStringsToChars(types(itype)), '.pdf']); 
end