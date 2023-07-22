%%%% DANGEROUS! I don't have time to write this to a function properly. 
% Set up colorbar stuff for the script b1_plot_paper_maps_w_litho_mld.m


% Colorbar and label
ax = gca(); 
cbar = colorbar('Location', 'south'); 
cbar.Position(3) = ax.Position(3) * .4; 
cbar.Position(1) = (ax.Position(1)+ax.Position(3)*.5) ; 
if ifig == 1; 
    cbar.Position(1) = cbar.Position(1) - 0.01; % For some reason first axis is off. 
end