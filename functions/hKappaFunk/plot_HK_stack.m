function []=plot_HK_stack(HK_H, HK_K, HK_A, options)
    arguments
        HK_H
        HK_K
        HK_A
        options.model_vpvs = nan; 
        options.model_zmoh = nan; 
        options.title = ''; 
        options.figNum = 198; 
        options.saveString = ''; 
    end
figure(options.figNum); clf; hold on; set(gcf,'color','white');
subplot(1,1,1); hold on; 
xlabel('kappa'); 
ylabel('H'); 
title(options.title); 
set(gca, 'ydir', 'reverse');         
sf = pcolor(HK_K, HK_H, HK_A'); %Can't use surf. It's 3d. It always covers other plotting objects. 
sf.EdgeAlpha = 0; 
colorbar(); 
xlim([min(HK_K), max(HK_K)]); 
ylim([min(HK_H), max(HK_H)]); 

% Simple wasy to estimate maximum h-k energy point. 
Emax = max(max(HK_A));
[xmax,ymax] = find(Emax==HK_A); 
kBest = HK_K(xmax); 
hBest = HK_H(ymax); 

% Plot position of max energy. 
scatter(gca, kBest, hBest, 50, 'red')
text(kBest, hBest, sprintf('H = %2.2f, k = %1.3f',...
    HK_H(ymax),HK_K(xmax)), 'verticalalignment', 'bottom' )

if ~isnan(options.model_vpvs); 
    % Plot true position of max energy. TODO add to legend. 
    scatter(options.model_vpvs, options.model_zmoh, 50, 'k'); 
end

if ~isempty(options.saveString); 
    exportgraphics(gcf, options.saveString ) 
end

end