function [ ] = plot_MISFIT_TRENDS(par,allmodels,misfits,resdir)

% Parameters you define. 
resolution = 200; % DPI. These are possibly supplemental figures. 
ylimits    = [10^(-2.5), 2]; % Set manually based on what you expect to be the highest versus lowest error and sigma. 
position   = [-1439 534 802 509]; % Position of figure. Same for both. 
position   = [2 10 863 1327]; 
ndattyps   = length(par.inv.datatypes); 
n_rows     = 6;
n_cols     = 3; 
dot_size   = 5; 
box_line_width = 0.8; 

% Increase n_rows if needed
if (n_rows * n_cols) < ndattyps; 
    n_rows = ceil(ndattyps / n_cols); 
end

% Plot RMS mismatch of each data type. 
figure(156), clf, set(gcf,'pos',position), hold on; 
h = tiledlayout(n_rows, n_cols, 'TileSpacing', 'tight', 'Padding', 'compact');% Tight, compact, and normal. 
sgtitle('RMS Mismatch')

for idat = [1:ndattyps]; 
    thisdat = par.inv.datatypes{idat}; 
    
    ax=nexttile(idat, [1,1]); 
    hold on; box on; set(gca, 'linewidth', box_line_width); 
    ylim(ylimits); 
    set(gca, 'yscale', 'log');     
    
    title(replace(thisdat, '_', ' '), 'fontweight', 'normal');

    for ichain = [1:length(misfits)]; 
        misfiti = misfits{ichain}; 
         
        thisval = [misfiti.rms.(thisdat)]'; 
        iter = misfiti.iter; 
       
        scatter(iter, thisval, dot_size, 'filled');
    end
end

% exportgraphics(gcf, [resdir '/datatype_misfits_rms.png' ], 'Resolution', 130); 
exportgraphics(gcf, [resdir '/datatype_misfits_rms.jpeg'], 'Resolution', resolution); % JPEG tends to look better for same KB size. 

% Plot sigma, or inverted standard deviation. 
figure(157), clf, set(gcf,'pos',position), hold on;  
h = tiledlayout(n_rows, n_cols, 'TileSpacing', 'tight', 'Padding', 'compact');% Tight, compact, and normal. 
sgtitle('\sigma (Inverted standard deviation)')

for idat = [1:ndattyps]; 
    thisdat = par.inv.datatypes{idat}; 
    
    ax=nexttile(idat, [1,1]); 
    hold on; box on; set(gca, 'linewidth', box_line_width); 
    ylim(ylimits); 
    set(gca, 'yscale', 'log');     
    
    title(replace(thisdat, '_', ' '), 'fontweight', 'normal');

    for ichain = [1:length(allmodels)]; 
        allmodelsi = allmodels{ichain}; 
        all_sigma = [allmodelsi.datahparm]; 
        thisval = [all_sigma.(['sig_' thisdat])]'; 
        iter = [allmodelsi.iter]';        
        scatter(iter, thisval, dot_size, 'filled');
    end
end

exportgraphics(gcf, [resdir '/datatype_misfits_sigma.jpeg'], 'Resolution', resolution);


% Plot likelihood corresponding to each datatype. 
figure(157), clf, set(gcf,'pos',position), hold on;  
h = tiledlayout(n_rows, n_cols, 'TileSpacing', 'tight', 'Padding', 'compact');% Tight, compact, and normal. 
sgtitle('Likelihood (log) contribution')

for idat = [1:ndattyps]; 
    thisdat = par.inv.datatypes{idat}; 
    
    ax=nexttile(idat, [1,1]); 
    hold on; box on; set(gca, 'linewidth', box_line_width); 
%     ylim(ylimits); 
%     set(gca, 'yscale', 'log');     
    
    title(replace(thisdat, '_', ' '), 'fontweight', 'normal');

    for ichain = [1:length(misfits)]; 
        misfiti = misfits{ichain}; 
         
        thisval = [misfiti.logL_indivdat.(thisdat)]'; 
        iter = misfiti.iter; 
       
        scatter(iter, thisval, dot_size, 'filled');
    end
end

exportgraphics(gcf, [resdir '/datatype_misfits_loglik.jpeg'], 'Resolution', resolution);
    
end