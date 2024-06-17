function []=plot_sensitivity_kernel_hv(k, options); 
    arguments
        k
        options.dat=[] % Structure containing vector of HVr (HV ratio) and periods. Like struct('HVr', hv_vector, 'periods', prediod_vector); 
        options.z_limits=[-2,300]
        options.model=[]; 
        options.predata=[]; 
        options.filename=[]; 
    end
    
dat   = options.dat; 
model = options.model; 
Np    = length(k); 

figure(88); clf; hold on;
set(gcf,'pos',[-1347 303 1323 713], 'color', 'white'); 
each_var = {'Kzh_Vs','Kzh_Vp','Kzh_rho'}; % Also have 'Kph_Vs','Kph_Vp','Kph_rho',
h=tiledlayout(1,length(each_var)+2,'TileSpacing','tight'); 

% Plot sensitivity kernels for each variable. 
for ivar = 1:length(each_var); 
    this_var = each_var{ivar}; 
    nexttile(ivar); cla; hold on; box on; 
    set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
    ylim(options.z_limits); 
    grid on; 
    title(strrep(this_var, '_', ' '), 'fontweight', 'normal'); 
    for ip = 1:Np
        plot(k{ip}.(this_var)/1e3,0.5*(k{ip}.Z1+k{ip}.Z2),...
            '-','linewidth',2, 'DisplayName',num2str(k{ip}.period));
    end  
    
    if ivar == 1; 
        leg=legend('Location', 'southeast'); 
        ylabel('Depth (km)'); 
    else
        set(gca, 'YTickLabel', []); 
    end
    
    xlim([-15e-5,15e-5]); 
end

% If model is provided, plot it. 
if ~isempty(model); 
    ivar = ivar + 1; 
    nexttile(ivar); cla; hold on; box on; set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
    title('Model', 'fontweight', 'normal'); 
    ylim(options.z_limits); 
    xlim([1, 9.5]); 
    grid on; 
    ylabel('Depth (km)'); 
    plot(model.VS , model.z, 'linewidth', 2, 'DisplayName', 'Vs' ); 
    plot(model.VP , model.z, 'linewidth', 2, 'DisplayName', 'Vp' ); 
    plot(model.rho, model.z, 'linewidth', 2, 'DisplayName', 'rho'); 
    legend('Location', 'best'); 
end

% Plot data calculated using full forward model. 
if ~isempty(dat); 
    ivar = ivar + 1; 
    nexttile(ivar); cla; hold on; box on; set(gca,'ydir', 'reverse', 'LineWidth', 1.5); 
    title('Ellipticity', 'fontweight', 'normal'); 
    ylabel('Period'); 
    xlabel('H/V'); 
    xlim([0.6, 1.2]); 
    ylim([min(dat.periods)-1.5, max(dat.periods)+5]); 
    set(gca, 'yscale', 'log');
    grid off; grid on;  % Log scale requires turning grid off before on
    plot(dat.HVr, dat.periods,...
        'linewidth', 2); 
end

if ~isempty(options.filename); 
    exportgraphics(gcf, options.filename); 
end

end