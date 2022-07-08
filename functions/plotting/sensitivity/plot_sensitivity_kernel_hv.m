function []=plot_sensitivity_kernel_hv(k, options); 
    arguments
        k
        options.z_limits=[-2,300]
        options.model=[]; 
        options.predata=[]; 
        options.filename=[]; 
    end
model=options.model; 
predata=options.predata; 

Np = length(k); 

figure(88); clf; hold on;
set(gcf,'pos',[-1347 303 1323 713], 'color', 'white'); 
each_var = {'Kzh_Vs','Kzh_Vp','Kzh_rho'}; % Also have 'Kph_Vs','Kph_Vp','Kph_rho',
h=tiledlayout(1,length(each_var)+2,'TileSpacing','tight'); 

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
end

if ~isempty(model); 
    nexttile(ivar+1); cla; hold on; box on; set(gca, 'ydir', 'reverse', 'LineWidth', 1.5); 
    title('Model', 'fontweight', 'normal'); 
    ylim(options.z_limits); 
    set(gca, 'YTickLabel', []); 
    plot(model.VS , model.z, 'linewidth', 2, 'DisplayName', 'Vs' ); 
    plot(model.VP , model.z, 'linewidth', 2, 'DisplayName', 'Vp' ); 
    plot(model.rho, model.z, 'linewidth', 2, 'DisplayName', 'rho'); 
    legend('Location', 'best'); 
    grid on; 
end

if ~isempty(predata); 
    nexttile(ivar+2); cla; hold on; box on; set(gca,'ydir', 'reverse', 'LineWidth', 1.5); 
    title('Ellipticity', 'fontweight', 'normal'); 
    ylabel('Period'); 
    plot(predata.SW_HV.HVr, predata.SW_HV.periods,...
        'linewidth', 2, 'DisplayName', 'rho'); 
    grid on;
end

if ~isempty(options.filename); 
    exportgraphics(gcf, options.filename); 
end

end