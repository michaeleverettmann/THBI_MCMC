clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
fpdfs    = sprintf('%scompiled_pdfs_%s.mat',out_dir,STAMP); % File with pdfs. a2_1...m

mdls = load(fresults).mdls; 
pdfs = load(fpdfs).pdfs; 
load('surface_out_example.mat');


%%
fnum = 101; 
lolim = [-87, -76]; 
lalim = [35, 43]; 
Q1 = [lalim(2), lolim(1)];
Q2 = [lalim(1), lolim(2)]; 
offsecmax = 1; %5%  distance off section allowed, in degrees
[profd,profaz] = distance(Q1(1),Q1(2),Q2(1),Q2(2));


by_line = logical(zeros(size(mdls.lon))); 
d_perp  = zeros(        size(mdls.lon)) ; 
d_par   = zeros(        size(mdls.lon)) ; 
for ista = 1:length(by_line); 

    [d_perp(ista),d_par(ista)] = dist2line_geog( Q1,Q2,...
        [mdls.lat(ista),mdls.lon(ista)],0.1 );    

    if d_perp(ista) <= offsecmax; 
        by_line(ista) = true; 
    end

end


scale_pdf = .1; 
figure(fnum); clf; hold on; 
set(gcf, 'pos', [593 570 832 373], 'color', 'white'); 
xlabel('Distance along section (degree)'); 
ylabel('Vs (km/s'); 
title('P(Vs) versus surface Vs', 'FontWeight','normal'); 
for ista = 1:length(by_line); 
    if ~ by_line(ista); continue; end

    mm = pdfs(ista).mm; 
    pm = pdfs(ista).pm; 

    xbase = d_par(ista) .* ones(size(pm));
    pm = pm * scale_pdf ; 

    color = pm; 
%     plot(xbase, mm, 10)
    plot(xbase + pm, mm, 'LineWidth', 2, 'Color', 'k')
end
grid on; 
box on; 

n_surf_pts = 100; 
gcarc = linspace(min(d_perp)-1, max(d_perp)+1, n_surf_pts)'; 
azline = profaz; 
[lat_surf_line, lon_surf_line] = reckon(Q1(1), Q1(2), gcarc, azline); 
[vs_surf_line] = griddata(longrid, latgrid, vgrid_out, lon_surf_line, lat_surf_line); 

figure(fnum); % Reopen this figure. 
plot(gcarc, vs_surf_line, 'color', 'blue', 'LineWidth', 3)


a3_2_plot_surface_simple(llminmax, 'stalon', mdls.lon(by_line),...
    'stalat', mdls.lat(by_line), 'stav', zeros(sum(by_line),1)); 
[x_surf_line, y_surf_line]=m_ll2xy(lon_surf_line, lat_surf_line); 
plot(x_surf_line, y_surf_line, 'color', 'blue', 'LineWidth', 4); 
set(gcf, 'pos', [593 570 832 373], 'color', 'white'); 



exportgraphics(figure(fnum), './surface_versus_pdf.pdf'); 
exportgraphics(figure(1   ), './surface_versus_pdf_mapview.pdf'); 

