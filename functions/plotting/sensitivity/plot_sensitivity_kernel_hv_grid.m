function []=plot_sensitivity_kernel_hv(k, options); 
    arguments
        k
        options.dat=[] % Structure containing vector of HVr (HV ratio) and periods. Like struct('HVr', hv_vector, 'periods', prediod_vector); 
        options.z_limits=[-2,300]
        options.model=[]; 
        options.predata=[]; 
        options.filename=[]; 
    end
% NOt quite finished. 
    
dat   = options.dat; 
model = options.model; 
Np    = length(k); 

figure(88); clf; hold on;
set(gcf,'pos',[-1347 303 541 740], 'color', 'white'); 
each_var = {'Kzh_Vs','Kzh_Vp','Kzh_rho'}; % Also have 'Kph_Vs','Kph_Vp','Kph_rho',
h=tiledlayout(length(each_var),1,'TileSpacing','tight'); 

% Plot sensitivity kernels for each variable. 
for ivar = 1:length(each_var); 
    this_var = each_var{ivar}; 
    nexttile(ivar); cla; hold on; box on; 
    set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
    title(strrep(this_var, '_', ' '), 'fontweight', 'normal'); 

    zgrid = 0.5*(k{1}.Z1+k{1}.Z2); % [0:300]'; 
    sgrid = dat.periods; 
    pgrid = zeros(length(zgrid), length(sgrid)); 
    

    for ip = 1:Np
        sens = k{ip}.(this_var)/1e3
        pgrid(:,ip) = sens; 
    end  
    colorbar(); 
    colormap(redbluecmap()); 
    contourf(sgrid, zgrid, pgrid,  'LineColor', 'none'); 
    xlabel('Period (s)'); 
    ylabel('Depth (km)'); 
        
    set(gca, 'xscale', 'log'); 
    ylim([0, 50]); 
    caxis([-1e-4, 1e-4]); 
end

if ~isempty(options.filename); 
    exportgraphics(gcf, options.filename); 
end

end