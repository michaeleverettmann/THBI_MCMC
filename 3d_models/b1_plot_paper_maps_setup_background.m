% TODO this should be put into a function
% Set up map default stuff for the script b1_plot_paper_maps_w_litho_mld.m

% State lines
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end

% Coast lines
cst = m_coast('patch',[1 1 1], 'FaceAlpha', 0); 

%%%% Tectonic stuff
% Tectonic fronts. 
m_plot(app_bord(1,:), app_bord(2,:), 'linewidth', 1.5, 'color', color_front); 
% m_plot(gre_bord(1,:), gre_bord(2,:), 'linewidth', 3, 'color', color_front); 

% Load tectonic files to plot. 
plt_ftrs = ["grv_frt"]; % Which things to plot. Not showing thrusts because not pertinent. Those they are mislabeled as mislabeled as  "something_province" % , "MCR", "Reelfoot"
plt_ftrs_path = tectpath1 + plt_ftrs + "*.txt";
for iftr = 1:length(plt_ftrs); 
    fplot = ls(plt_ftrs_path(iftr)); 
    fplot = splitlines(fplot(1:end-1)); % Need to remove last line, to prevent nonsense empty extra cell index. 
    this_color = clr_tectfiles.(plt_ftrs(iftr)); 
    for ifplot = 1:length(fplot); 
        xy = load(fplot{ifplot});  
        x = xy(:,1); 
        y = xy(:,2); 
        m_plot(x, y, 'Color', this_color, 'linewidth', 1.5); % TODO color for each feature
    end
end

% Label features. 
if ifig == 1; 
    m_text(-85, 39, 'GF', 'color', color_text, 'fontsize', 12, ...
        'rotation', 75, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontweight', 'bold'); 
    m_text(-80.8, 39, 'AF', 'color', color_text, 'fontsize', 12, ...
        'rotation', 55, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontweight', 'bold'); 
end


%%%% Labels, ticks, states, etc. 
% Label the plot
m_text(-74, 30.5, label, 'Units','data', 'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom'); 

% Handle ticks across different subplots. 
xtck = [-88:6:-68]; 
ytck = [28:4:44]; 
xtckstr = xtck; 
ytckstr = ytck; 

% Ticks for 3x2
if main_paper_figure % Main paper
    xtckcel = {[], [], [], [], [], [], xtckstr, xtckstr, xtckstr}; 
    ytckcel = {ytckstr, [], [], ytckstr, [], [], ytckstr, [],[]}; 
else % Supplemental figure
    xtckcel = {xtckstr, xtckstr, xtckstr, xtckstr, xtckstr, xtckstr, xtckstr, xtckstr, xtckstr}; % Manually make 3x3. Sorry, future self. 
    ytckcel = {ytckstr, ytckstr, ytckstr, ytckstr, ytckstr, ytckstr, ytckstr, ytckstr, ytckstr}; 
end

% Grid, coastlines I think, etc. Background colors. 
m_grid('box','fancy','linestyle','none','gridcolor',.8 .*[1,1,1],...
    'backcolor',.9.*[1,1,1], 'xticklabel', xtckcel{ifig}, 'yticklabel', ytckcel{ifig},... % 'backcolor',[.3 .75 1]
    'xtick', xtck, 'ytick', ytck);

% Figure subletter labels - do this late in script so letters don't get covered. 
subplot_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P']; 
txt= text(0.04, 0.015, subplot_letters(ifig), 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment','bottom', ...
    'FontSize',13, 'Color',[0,0,0]); 
