function [ax1, ax2] = start_hk_ax(ax1, kStack, hStack, eStack)
% Function to help with other function plot_h_kappa_progress.m 


% Copy basic axis. One axis for h-kappa stack, one for model locations. 
% Need two axes for different colorbars. 


set(ax1, 'xlim', [min(min(kStack)), max(max(kStack))]); 
set(ax1, 'ylim', [min(min(hStack)), max(max(hStack))]); 
ax2 = copyobj(ax1, gcf); % Displayed in middle order.  
hold on; 
set(ax2, 'color', 'none'); 
linkaxes([ax1, ax2]); 
hold on; 

ax2.Visible = 'on'; % Only show what I plot here. Not labels, etc. 
ax1.Visible = 'off' ; 
box(ax2); 

% Plot k-kappa stack
axes(ax1); 
h = pcolor(ax1, kStack, hStack, eStack'); 
set(h, 'EdgeColor', 'none');

set(ax1, 'YDir', 'reverse'); 
set(ax2, 'YDir', 'reverse'); % 

axes(ax1); 
axes(ax2); 

end