clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
fpdfs = sprintf('%scompiled_pdfs_%s.mat',out_dir,STAMP); % File with pdfs. a2_1...m

version_surf = 6; 
n_surf_pts = 100; 

lolim = [-87, -76; -86, -68; -88, -78; -87, -80.5]; 
lalim = [ 43,  35;  30,  47;  36,  33;  38,  25  ]; 
fnum = 101; 
scale_pdf = .035; 
offsecmax = 1.5; %5%  distance off section allowed, in degrees
% mlim_manual = [4.1, 5]; 


v_at_depth = true; % Use velocity from a depth, or one of the other parameters like moho depth. 
% z_vs loaded in a0....m
% for idep = 1:length(z_vs);


for idep = int16([20, 65]/5); 

depth = z_vs(idep); 

this_inversion = sprintf('vs%1.0f',depth); % String name affiliated with figures and files. 
iz = find(z_vs == depth); 

mdls = load(fresults).mdls; 
% pdfs = load(fpdfs).pdfs; 
sfsmat = load('surface_out_example.mat'); xgrid = sfsmat.xgrid; ygrid = sfsmat.ygrid; llminmax = sfsmat.llminmax; latgrid = sfsmat.latgrid; longrid = sfsmat.longrid; 
sfsmat2= load(sprintf('%s/surface_values_V%1.0f', this_inversion, version_surf)); vgrid_out = sfsmat2.mgrid_out; 


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


%% Prep pdf-section style figure
figure(fnum); clf; hold on; 
% set(gcf, 'pos', [1053 564 767*2 329*ceil(.5*size(lolim,1))], 'color', 'white'); %2x2
% tiledlayout(ceil(.5 * size(lolim, 1 )), 2,'TileSpacing', 'Compact'); %2x2
set(gcf, 'pos', [1053 564 767 329*size(lolim,1)])
tiledlayout(size(lolim, 1 ), 1,'TileSpacing', 'Compact')


%% Get station data and plot it. Loop over different x sections. 
lat_surf_line_all = zeros(n_surf_pts, size(lolim,1) ); 
lon_surf_line_all = lat_surf_line_all; 
all_ms = []; 
all_pms = []; 
all_cpms = []; 
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

%%
section_letter = char(64+i_xsect); % Text for cross-section name. ith letter of alphabet
figure(fnum); 
nexttile(); hold on; 
set(gca, 'LineWidth', 1.5); 
xlabel('Distance along section (degree)'); 
ylabel('m'); 
title('P(m) versus inverted m', 'FontWeight','normal'); 
grid on; 
box on; 
t1=text(0.02, .98, section_letter    , 'fontsize', 20, 'color', 'r', 'units', 'normalized', 'VerticalAlignment','top'); 
t2=text(0.98, .98, section_letter+"'", 'fontsize', 20, 'color', 'r', 'units', 'normalized', 'VerticalAlignment','top', 'HorizontalAlignment','right'); 

%% Plot station pdfs. 
[junk, sta_plot_order] = sort(d_perp); % Order stations from furthest to closest. 
plot_these_stas = sta_plot_order(junk<offsecmax); % Plot furthest stations first, then closest. 


for ista = plot_these_stas(end:-1:1)'; 

    mm = pdfs(ista).mm; 
    pm = pdfs(ista).pm; 

    cpm = cumtrapz(mm, pm); 
    median_mm = linterp(cpm, mm, 0.5); 

    keep_p = (cpm > 0.001) & (cpm < 0.999); % Get rid of 99th percentile to help clean up the plots. 
    mm = mm(keep_p); 
    pm = pm(keep_p); 
    cpm = cpm(keep_p); 
    mm = [mm(1), mm, mm(end)]; 
    pm = [0, pm, 0]; 
    cpm = [0, cpm, 1]; 

    all_ms = [all_ms, mm];
    all_pms = [all_pms, pm]; 
    all_cpms = [all_cpms, cpm]; 

    xbase = d_par(ista) .* ones(size(pm));
    pm = pm * scale_pdf ; 

    disti = d_perp(ista); 
    dist_rat = disti / offsecmax; 
%     color = dist_rat * ones(1,3);% .* ones(2*length(mm),3); 
    face_alpha = 1-dist_rat; 
    face_alpha = face_alpha ^ 2; 
    color = zeros(1,3); 
    edge_alpha = 1 * (disti < 0.25);

    fill([xbase + pm; xbase - pm(end:-1:1)]', [mm; mm(end:-1:1)]', color, ...
        'FaceAlpha', face_alpha, ...
        'EdgeColor', 'k', 'EdgeAlpha', edge_alpha); % mean([face_alpha,1])  
    scatter(xbase(1), median_mm, 25, 'r', 'diamond', 'filled', 'MarkerFaceAlpha', face_alpha); 
end

% ylim([])

%% Interpolate model surface to our line section and plot. 
gcarc = linspace(min(d_par(by_line)), max(d_par(by_line)), n_surf_pts)'; 
azline = profaz; 
[lat_surf_line, lon_surf_line] = reckon(Q1(1), Q1(2), gcarc, azline); 
lat_surf_line_all(:,i_xsect) = lat_surf_line; 
lon_surf_line_all(:,i_xsect) = lon_surf_line; 
[vs_surf_line] = griddata(longrid, latgrid, vgrid_out, lon_surf_line, lat_surf_line); 

figure(fnum); % Reopen this figure. 
plot(gcarc, vs_surf_line, 'color', 'blue', 'LineWidth', 3); 
% ylim([4.1, 5.0]); 

%% Plot position of line section on the map. 
% % % figure(1); 
% % % [stax, stay] = m_ll2xy(mdls.lon(by_line), mdls.lat(by_line));       
% % % scatter(stax, stay, 30, 'filled', 'MarkerEdgeColor','k'); 
% % % [x_surf_line, y_surf_line]=m_ll2xy(lon_surf_line, lat_surf_line); 
% % % plot(x_surf_line, y_surf_line, 'color', 'blue', 'LineWidth', 5); 
% % % t1=text(x_surf_line(1  ), y_surf_line(1  ), section_letter    , 'fontsize', 20, 'color', 'r');
% % % t2=text(x_surf_line(end), y_surf_line(end), section_letter+"'", 'fontsize', 20, 'color', 'r');





end % End loop on different sections. 

% min_ms = min(all_ms); 
% max_ms = max(all_ms); 
% min_y = max(all_ms(all_cpms < 0.01)); 
% max_y = min(all_ms(all_cpms > 0.99)); 
min_y = mean(all_ms(all_cpms < 0.01)); 
max_y = mean(all_ms(all_cpms > 0.99)); 
y_range = max_y - min_y; 
figure(fnum); 
% linkaxes(gcf().Children.Children, 'y'); 
linkaxes(gcf().Children.Children); 
ylim([min_y - .2 * y_range, max_y + .2 * y_range]); 


%% Save figures. 
exportgraphics(figure(fnum), sprintf('%s/surface_versus_pdf_V%1.0f.pdf', ...
    this_inversion, version_surf),'ContentType','vector'); 
a3_2_plot_surface_simple(llminmax, 'stalon', mdls.lon, 'stalat', mdls.lat, ...
    'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', vgrid_out,...
    'sectlon', lon_surf_line_all, 'sectlat', lat_surf_line_all); 
exportgraphics(figure(1   ), sprintf('%s/surface_versus_pdf_mapview_V%1.0f.pdf',...
    this_inversion, version_surf)); 
% exportgraphics(gcf, sprintf('xsections/xsections_map_V%1.0f.pdf', version_surf)); 


end
