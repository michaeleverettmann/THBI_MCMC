clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
fpdfs = sprintf('%scompiled_pdfs_%s.mat',out_dir,STAMP); % File with pdfs. a2_1...m


v_at_depth = true; % Use velocity from a depth, or one of the other parameters like moho depth. 
% z_vs loaded in a0....m
% for idep = 1:length(z_vs);
for idep = [2, 13, 30, 50]; 


depth = z_vs(idep); 

this_inversion = sprintf('vs%1.0f',depth); % String name affiliated with figures and files. 
iz = find(z_vs == depth); 

mdls = load(fresults).mdls; 
% pdfs = load(fpdfs).pdfs; 
sfsmat = load('surface_out_example.mat'); xgrid = sfsmat.xgrid; ygrid = sfsmat.ygrid; llminmax = sfsmat.llminmax; latgrid = sfsmat.latgrid; longrid = sfsmat.longrid; 
sfsmat2= load(sprintf('%s/surface_values', this_inversion)); vgrid_out = sfsmat2.vgrid_out; 


pdf_file = load(fpdfs); 
pdfs = pdf_file.pdfs_allparm; 
%%% Put this in function later. 
pdfs_vs = pdfs(1).vs{1}; % Make a new structure (obnoxious). And have to start with the correct field names. Reason for new structure is that, I used a cell array for each different depth. Matlab doesn't actually access the nth stations ith cell array all in one call. 
nsta = length(pdfs); 
for ista = 1:nsta
    pdfs_vs(ista) = pdfs(ista).vs{iz}; 
end
pdfs = pdfs_vs; 
%%% Put this in function later. 


%%
lolim = [-87, -76; -86, -72; -88, -76]; 
lalim = [ 43,  35;  30,  44;  36,  32]; 
fnum = 101; 
scale_pdf = .035; 
% offsecmax = 0.65; %5%  distance off section allowed, in degrees
offsecmax = 1.5; %5%  distance off section allowed, in degrees
mlim_manual = [4.1, 5]; 

%% Prep pdf-section style figure
figure(fnum); clf; hold on; 
set(gcf, 'pos', [1053 564 767*2 329*ceil(.5*size(lolim,1))], 'color', 'white'); 
tiledlayout(ceil(.5 * size(lolim, 1 )), 2,'TileSpacing', 'Compact')

%% Prep map view figure
figure(1); clf; hold on; 
lonmin = llminmax(1); lonmax = llminmax(2); latmin = llminmax(3); latmax = llminmax(4); 
set(gcf, 'color', 'white'); 
m_proj('lambert', 'long',[lonmin, lonmax],'lat',[latmin, latmax]);
m_coast('patch',[1 1 1]); % m_coast('patch',[1 .85 .7]);
m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[.3 .75 1]);
title('Cross-sections', 'fontweight', 'normal')
contourf(xgrid, ygrid, vgrid_out, 15,...
    'LineStyle','none'); 
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end
colorbar(); turbo_map = turbo(); turbo_map = turbo_map(end:-1:1,:); colormap(turbo_map); 
% caxis(mlim_manual); 
% caxis([4.35, 4.7]); 

%% Get station data and plot it. Loop over different x sections. 
for i_xsect = 1:size(lolim, 1); 

Q1 = [lalim(i_xsect, 1), lolim(i_xsect, 1)];
Q2 = [lalim(i_xsect, 2), lolim(i_xsect, 2)]; 
[profd,profaz] = distance(Q1(1),Q1(2),Q2(1),Q2(2));


by_line = logical(zeros(size(mdls.lon))); 
d_perp  = zeros(        size(mdls.lon)) ; 
d_par   = zeros(        size(mdls.lon)) ; 
d_lon   = zeros(        size(mdls.lon)) ; 
for ista = 1:length(mdls.lon); 
    [d_perp(ista),d_par(ista)] = dist2line_geog( Q1,Q2,...
        [mdls.lat(ista),mdls.lon(ista)],0.0025 );   
    d_lon(ista) = mdls.lon(ista); 

    if d_perp(ista) <= offsecmax; 
        by_line(ista) = true; 
    end
end

section_letter = char(64+i_xsect); % Text for cross-section name. ith letter of alphabet
figure(fnum); 
nexttile(); hold on; 
set(gca, 'LineWidth', 1.5); 
xlabel('Distance along section (degree)'); 
ylabel('Vs (km/s'); 
title('P(Vs) versus surface Vs', 'FontWeight','normal'); 
grid on; 
box on; 
t1=text(0.02, .98, section_letter    , 'fontsize', 20, 'color', 'r', 'units', 'normalized', 'VerticalAlignment','top'); 
t2=text(0.98, .98, section_letter+"'", 'fontsize', 20, 'color', 'r', 'units', 'normalized', 'VerticalAlignment','top', 'HorizontalAlignment','right'); 

%% Plot station pdfs. 
% % % [junk, sta_plot_order] = sort(-d_perp); % Order stations from furthest to closest. 
% % % % plot_these_stas = sta_plot_order(by_line); % Plot furthest stations first, then closest. 
% % % plot_these_stas = plot_these_stas(d_perp(sta_plot_order) > offsecmax)
[junk, sta_plot_order] = sort(d_perp); % Order stations from furthest to closest. 
plot_these_stas = sta_plot_order(junk<offsecmax); % Plot furthest stations first, then closest. 
% plot_these_stas = plot_these_stas(d_perp(sta_plot_order) > offsecmax)

% by_line = sta_plot_order(d_perp(sta_plot_order)<offsecmax); 

for ista = plot_these_stas(end:-1:1)'; 

    mm = pdfs(ista).mm; 
    pm = pdfs(ista).pm; 

    cpm = cumtrapz(mm, pm); 
    median_mm = linterp(cpm, mm, 0.5); 

    xbase = d_par(ista) .* ones(size(pm));
    pm = pm * scale_pdf ; 

    disti = d_perp(ista); 
    dist_rat = disti / offsecmax; 
%     color = dist_rat * ones(1,3);% .* ones(2*length(mm),3); 
    face_alpha = 1-dist_rat; 
    face_alpha = face_alpha ^ 2; 
    color = zeros(1,3); 
    edge_alpha = 1 * (disti < 0.25);

    fill([xbase + pm; xbase - pm(end:-1:1)], [mm; mm(end:-1:1)], color, ...
        'FaceAlpha', face_alpha, ...
        'EdgeColor', 'k', 'EdgeAlpha', edge_alpha); % mean([face_alpha,1])  
    scatter(xbase(1), median_mm, 25, 'r', 'diamond', 'filled', 'MarkerFaceAlpha', face_alpha); 
end

%% Interpolate model surface to our line section and plot. 
n_surf_pts = 100; 
gcarc = linspace(min(d_par(by_line)), max(d_par(by_line)), n_surf_pts)'; 
azline = profaz; 
[lat_surf_line, lon_surf_line] = reckon(Q1(1), Q1(2), gcarc, azline); 
[vs_surf_line] = griddata(longrid, latgrid, vgrid_out, lon_surf_line, lat_surf_line); 

figure(fnum); % Reopen this figure. 
plot(gcarc, vs_surf_line, 'color', 'blue', 'LineWidth', 3); 
% ylim([4.1, 5.0]); 

%% Plot position of line section on the map. 
figure(1); 
[stax, stay] = m_ll2xy(mdls.lon(by_line), mdls.lat(by_line));       
scatter(stax, stay, 30, 'filled', 'MarkerEdgeColor','k'); 
[x_surf_line, y_surf_line]=m_ll2xy(lon_surf_line, lat_surf_line); 
plot(x_surf_line, y_surf_line, 'color', 'blue', 'LineWidth', 5); 
t1=text(x_surf_line(1  ), y_surf_line(1  ), section_letter    , 'fontsize', 20, 'color', 'r');
t2=text(x_surf_line(end), y_surf_line(end), section_letter+"'", 'fontsize', 20, 'color', 'r');





end % End loop on different sections. 

%% Save figures. 
exportgraphics(figure(fnum), sprintf('%s/surface_versus_pdf.pdf', this_inversion)); 
exportgraphics(figure(1   ), sprintf('%s/surface_versus_pdf_mapview.pdf', this_inversion)); 

end
