function plot_h_kappa_progress(trudata, allmodels, resdir, chainNo)

% Extract some info
hStack = trudata.HKstack_P.H; 
kStack = trudata.HKstack_P.K;  
eStack = trudata.HKstack_P.Esum; 
hIter = [allmodels.zmoh]'; 
kIter = [allmodels.vpvs]'; 
iIter = [allmodels.iter]'; 



% start basic figure properties. 
figure(1); clf; hold on; 
set(gcf, 'pos', [1554 413 1458 1303]); 
% xlim([min(k), max(k)]); 
% ylim([min(h), max(h)]); 
set(gca, 'ydir', 'reverse'); 

% Copy basic axis. One axis for h-kappa stack, one for model locations. 
% Need two axes for different colorbars. 
ax1 = gca; 
ax2 = copyobj(ax1, gcf); % Displayed in middle order.  
ax1.Visible = 'on' ; 
ax2.Visible = 'off'; % Only show what I plot here. Not labels, etc. 
box(ax1); 

% Plot k-kappa stack
axes(ax1); 
surf(kStack, hStack, eStack', 'EdgeAlpha', 0); 
% colorbar('east'); 


% Plot iterations progress
axes(ax2);
% 

scatter(kIter, hIter, 30, iIter); 



% xlabel('kappa'); ylabel('H'); title('Synthetic H-kappa stack'); 
% set(gca, 'ydir', 'reverse'); 
% colorbar(); 
linkaxes([ax1, ax2]); 

linePlt = plot(kIter, hIter, 'k'); 

% iterLabels = round(linspace(1,length(hIter), 10)); 
iterLabels = [1, 100, 400, 1500, 3000, 8000]; 
iterLabels = iterLabels( iterLabels < max(length(hIter)) ); 
iterLabels(1) = 1; 
for iIterLabel = 1:length(iterLabels); 
    text( kIter(iterLabels(iIterLabel)),... 
          hIter(iterLabels(iIterLabel)),...
          num2str(iterLabels(iIterLabel)),...
         'FontSize', 15); 
end
% colorbar('westoutside'); 

exportgraphics(figure(1), [resdir '/h_k_progress_chain' num2str(chainNo) '.jpg']);  

end