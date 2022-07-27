function []=plot_sensitivity_kernel_ray(k, options); 
    arguments
        k
        options.dat=[] % Structure containing vector of data and periods. Like struct('HVr', hv_vector, 'periods', prediod_vector); 
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
each_var = {'Vsv', 'Vsh', 'Vpv', 'Vph', 'rho', 'eta'}; % Also have 'Kph_Vs','Kph_Vp','Kph_rho',
h=tiledlayout(1,length(each_var)+2,'TileSpacing','tight'); 

for ivar = 1:length(each_var); 
    this_var = each_var{ivar}; 
    nexttile(ivar); cla; hold on; box on; 
    set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
    ylim(options.z_limits); 
    grid on; 
    title(strrep(this_var, '_', ' '), 'fontweight', 'normal'); 
    for ip = 1:Np
        plot(k{ip}.(this_var),k{ip}.Z/1000,...
            '-','linewidth',2, 'DisplayName',num2str(k{ip}.period));
    end  
    
    if ivar == 1; 
        leg=legend('Location', 'southeast'); 
        ylabel('Depth (km)'); 
    else
        set(gca, 'YTickLabel', []); 
    end
end

if ~isempty(model); 
    ivar = ivar + 1; 
    nexttile(ivar); cla; hold on; box on; set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
    title('Model', 'fontweight', 'normal'); 
    ylim(options.z_limits); 
    grid on; 
    set(gca, 'YTickLabel', []); 
    plot(model.VS , model.z, 'linewidth', 2, 'DisplayName', 'Vs' ); 
    plot(model.VP , model.z, 'linewidth', 2, 'DisplayName', 'Vp' ); 
    plot(model.rho, model.z, 'linewidth', 2, 'DisplayName', 'rho'); 
    legend('Location', 'best'); 
end

if ~isempty(dat); 
    ivar = ivar + 1; 
    nexttile(ivar); cla; hold on; box on; set(gca,'ydir', 'reverse', 'LineWidth', 1.5); 
    title('Rayleigh Phv', 'fontweight', 'normal'); 
    ylabel('Period'); 
    xlabel('V (km/s)'); 
%     xlim(
%     ylim(
    set(gca, 'YScale', 'log'); 
    grid off; grid on; 
    plot(dat.phV,dat.periods,...
        'linewidth', 2)
end

if ~isempty(options.filename); 
    exportgraphics(gcf, options.filename); 
end

end