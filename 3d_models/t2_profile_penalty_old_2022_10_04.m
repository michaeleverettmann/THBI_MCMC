clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
fpdfs    = sprintf('%scompiled_pdfs_%s.mat',out_dir,STAMP); % File with pdfs. a2_1...m


mdls = load(fresults).mdls; 

%% % Reorganize structure. Could be more - or less - useful one way or other. 
%%% BORROWED FROM A2
models = struct(); 
zmoh          = zeros(length(mdls.nwk),1); 
zmohsig       = zeros(length(mdls.nwk),1);
zsed          = zeros(length(mdls.nwk),1);
zsedsig       = zeros(length(mdls.nwk),1);
xicr          = zeros(length(mdls.nwk),1); 
% xicrboundsmaybe= zeros(length(mdls.nwk),2);
vpvs          = zeros(length(mdls.nwk),1); 
% vpvssig       = zeros(length(mdls.nwk),1);


z_vs          = [0.5, 15, 50, 100, 150, 200]; 
vs            = zeros(length(mdls.nwk),length(z_vs)); 

for imd = 1:length(mdls.nwk); 

    % Make alternative format models variable. 
    models(imd).nwk   = mdls.nwk  {imd}; 
    models(imd).sta   = mdls.sta  {imd}; 
    models(imd).lat   = mdls.lat  (imd); 
    models(imd).lon   = mdls.lon  (imd); 
    models(imd).model = mdls.model{imd}; 
    models(imd).dir   = mdls.dir  {imd}; 
    
    % Extract map-view inversion results in cases where there's just one paramete for a station. 
    mdl = mdls.model{imd}; 
    zmoh              (imd) = mdl.zmohav; 
    zmohsig           (imd) = mdl.zmohsig; 
    zsed              (imd) = mdl.zsedav;     
    zsedsig           (imd) = mdl.zsedsig; 
    xicr              (imd) = mdl.xicrav;     
%     xicrboundsmaybe   (imd,:) = mdl.xicrsig2'; 
    vpvs              (imd) = mdl.vpvsav; 

    % Pseudo-tomography
    dist_from_z = (abs(z_vs - mdl.Z)); 
    [C,I ] = min(dist_from_z); 
    vs(imd,:) = mdl.VSav(I'); % Velocity at depths. 
end


%% Setting up coordinate system. 
lonmin = min(mdls.lon); 
latmin = min(mdls.lat); 
lonmax = max(mdls.lon); 
latmax = max(mdls.lat); 

llminmax = [lonmin, lonmax, latmin, latmax]; % For easy passing to plotting functions

figure(1); clf; hold on; 
m_proj('lambert', 'long',[lonmin, lonmax],'lat',[latmin, latmax]);
[stax, stay] = m_ll2xy(mdls.lon, mdls.lat);  % Get stax and stay!!!

%% Plot stations positions. 
a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, ...
    'fignum', 1, 'title', 'station positions')

%% Setting up surface. 
nx = 15; 
ny = 17;  
ivel = 4; % Temporary. Expand to a loop. 
edge_space = 0.2; % How far to go beyond max and min stax and y, in fraction. 


DX = max(stax) - min(stax); 
DY = max(stay) - min(stay); 
xline = linspace(min(stax) - edge_space * DX, max(stax) + edge_space * DX, nx)'; 
yline = linspace(min(stay) - edge_space * DY, max(stay) + edge_space * DY, ny)'; 
[xgrid, ygrid] = ndgrid(xline, yline);  
[longrid, latgrid] = m_xy2ll(xgrid, ygrid); 

vs_interp = griddata(stax, stay, vs(:,ivel), xgrid, ygrid, 'cubic'); 
% vs_interp(isnan(vs_interp)) = nanmean(vs_interp, 'all'); 

% Plot surface. 
a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', vs(:,ivel),...
    'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', vs_interp, ...
    'fignum', 2, 'title', 'Simple v interpolation');  
exportgraphics(gcf, 'surface_simple_interpolation.pdf'); 



%% Get PDF for stations. 
% Use an example pdf while developing surface inversion stuff. 
% fpdf = '~/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/pdf_example.mat'; 
% pdf = load(fpdf).pdf_example; 
% Probably something like pdf(ista) = load(pdf_sta). 
pdf_file = load(fpdfs); 
pdfs = pdf_file.pdfs; 

% rough_scale = 5e-11; % How much to penalize roughness. 
rough_scale = 5e-8; % How much to penalize roughness. 


% Make a starting model. 
vgrid = vs_interp; % Probably a fine starting model. 
vgrid(isnan(vgrid)) = nanmean(nanmean(vgrid)); 
vgrid = vgrid + randn(size(vgrid))-0.5; % Give some error

% Plot the starting model. 
a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', vs(:,ivel),...
    'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', vgrid, ...
    'fignum', 3, 'title', 'V starting model'); 

%% Set up grid roughness calculations. Use 2nd derivative smoothing. 
dx2 = ((xgrid(1:end-2,:) - xgrid(3:end,:))/2).^2; % Squared, so sign doesn't matter. 
dy2 = ((ygrid(:,1:end-2) - ygrid(:,3:end))/2).^2; 
x_dx2 = xgrid(2:end-1,:); y_dx2 = ygrid(2:end-1,:); % x and y positions where we have dx. I think just for plotting. 
y_dy2 = ygrid(:,2:end-1); x_dy2 = xgrid(:,2:end-1); % y and x positions where we have dy
% Note for optimizing: Use a scalar (not matrix) for dx and dy. But as is, it's more versatile for map projections. 

% Example roughness to make sure the calculations are good. 
% This is the calculation to get roughness throughout the inversion
% (Basically). 
dvdx2 = (vgrid(1:end-2,:) - 2*vgrid(2:end-1,:) + vgrid(3:end,:))./dx2; 
dvdy2 = (vgrid(:,1:end-2) - 2*vgrid(1,2:end-1) + vgrid(:,3:end))./dy2; 
dvdx2 = dvdx2.^2; % Square, mostly to prevent negative values...
dvdy2 = dvdy2.^2; 
roughness = sum(dvdx2, 'all') + sum(dvdy2, 'all'); 
roughness = roughness * rough_scale; % How to get roughness penalty value. 

a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay,...
    'xgrid', x_dx2, 'ygrid', y_dx2, 'vgrid', dvdx2, ...
    'fignum', 4, 'title', 'X roughness'); colorbar(); 
a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay,...
    'xgrid', x_dy2, 'ygrid', y_dy2, 'vgrid', dvdy2, ...
    'fignum', 5, 'title', 'Y roughness'); colorbar(); 

% How to get "modelled" velocity throughout inversion. 
vsta = interpn(xgrid, ygrid, vgrid, stax, stay, 'linear'); % MODELLED velocity at each station. 

% Probability associated with a velocity. 
pv_mod = zeros(size(stax)); 
for ista=1:length(pv_mod); 
    pv_mod(ista) = interp1(pdfs(ista).mm, pdfs(ista).pm, vsta(ista), 'linear', 0 ); % Probability of a velocity at a specific station, from our velocity model surface. 
end

penalty = - sum(pv_mod); % Lower penalty is higher probability. 
penalty = penalty + roughness; 

fhand=@(vgrid)a3_1_penalty(vgrid,...
    pdfs, rough_scale, dx2, dy2, xgrid, ygrid, stax, stay); % Supply constants to the function handle. 
fhand(vgrid); 
% vgrid_start = vgrid; 

fhand2=@()fhand(vgrid); 

% % % %%
% % % for iii=1:100; 
% % %     fhand2(); 
% % % end

%% Efficient interpolation from station to pdf

% Go from different mm at each sta to common mm at each sta through
% interpolation. 
nmm = 300; 
nsta = length(pdfs); 
pdf_terp = zeros(nmm, nsta); 


all_mm = [pdfs(:).mm]; 
min_mm = min(min(all_mm)); 
max_mm = max(max(all_mm)); 
mm_terp = linspace(min_mm, max_mm, nmm)'; % Corresponds to pdf_terp dimension 1. 
dmm_di = mm_terp(2) - mm_terp(1); % How much does mm change per index. 

for ista = 1:nsta; 
    pdf_terp_i = interp1(pdfs(ista).mm, pdfs(ista).pm, mm_terp,...
        'cubic', 0 ); 
    pdf_terp(:, ista) = pdf_terp_i; 
end

% Hypothetical vsta for sesting
% vsta_mod = rand(nsta,1) .* .5 + mean(mm_terp); 
% vsta_mod(1) = min(mm_terp); 
% vsta_mod(2) = mean(mm_terp); 
% vsta_mod(3) = max(mm_terp); % For easy testing. 
vsta_mod = linspace(min(mm_terp), max(mm_terp), nsta)'
% Check best_ind_flce([1,3]). If it's 1, half length, then the whole length, then index selection is probably good. 


% % mm_dist = pdf_terp - vsta_mod';
% % mm_dist = (mm_terp - vsta_mod').^2; % Array. 
% mm_dist = (mm_terp - vsta_mod'); 
% % [min_dist, ]=min(mm_dist); 
% % fhandsort=@()sort(mm_dist); 
% geq0 = mm_dist>0; 
% find(geq0, 1); 
% timeit(fhandsort)


% % % best_ind = 1 + (vsta_mod - min(mm_terp))./dmm_di; 
% % % best_ind_fl = floor(best_ind); 
% % % best_ind_ce = ceil(best_ind); 
% % % % weight_fl = 1./(best_ind_fl - best_ind).^2; 
% % % ind_weight = 1./([best_ind-best_ind_fl, best_ind-best_ind_ce].^2); % In order of floor, then ceil. 
% % % ind_weight = ind_diff ./ sum(ind_diff,2); 

% More efficient. No sorting. 

best_ind = 1 + (vsta_mod - min(mm_terp))/dmm_di; % Start at one. Find (non integer) best index of array based on dindex/dmm
best_ind_flce = [floor(best_ind), ceil(best_ind)]; 
best_ind_flce(best_ind_flce<1)=1; % If lower than our bounds, just use our lowest part of pdf. That should be 0 anyway!
best_ind_flce(best_ind_flce>size(mm_terp,1)) = size(mm_terp,1); % If greater than our bounds, just use our last part of pdf. That should be 0 anyway!
ind_weight = 1./([best_ind - best_ind_flce].^2); 
ind_weight(or(isnan(ind_weight),ind_weight>99999)) = 99999; % In case of 1./0. 
ind_weight = ind_weight ./ sum(ind_weight,2); 
% % % % pdf_terp_flce = pdf_terp([[1:nmm];[1:nmm]],best_ind_flce')
% % % pdf_terp_flce = zeros(size(best_ind_flce,2), nsta); 
% % % pdf_terp_flce(1,:)=pdf_terp(1:nsta,best_ind_flce(:,1)); 
% % % pdf_terp_flce(2,:)=pdf_terp(2,best_ind_flce(:,1)); 

%%% Sort of complicated array accessing
% @(row, col, nrow, ncol)row + col*()

%%%
sta_hack = [[1:nsta];[1:nsta]]'; % Which station each model parameter corresponds, in best_ind_flce
best_ind_linear = sub2ind([nmm, nsta], best_ind_flce, sta_hack); % Linear (single dimension) coordinates to points of interest in pdf_terp

% access_pdf_terp = logical(zeros(size(pdf_terp)));
% access_pdf_terp( best_ind_linear ) = true; % 

% pdf_terp_flce(1)
% access_pdf_terp = logical(zeros(size(pdf_terp)));
% access_pdf_terp(best_ind_flce',[[1:nsta];[1:nsta]]) = true; 
pdf_mod = sum(pdf_terp(best_ind_linear).*ind_weight,2); % Interpolated pdf between two mm values
% Alteratnvei, probably slower way. 
% % % pdf_mod_matmult = zeros(length(ind_weight), nsta); 
% % % pdf_mod_matmult(best_ind_flce) = ind_weight; 

% sparse_mult = sparse(best_ind_flce', [[1:nsta];[1:nsta]], 1)
% sparse_mult = logical(sparse(best_ind_flce', [[1:nsta];[1:nsta]], 1))
% fhandtest=@()logical(sparse(best_ind_flce', [[1:nsta];[1:nsta]], 1)); 
% sparse_mult = (sparse(best_ind_flce', [[1:nsta];[1:nsta]], ind_weight)); 
% fhandtest=@()(sparse(best_ind_flce', [[1:nsta];[1:nsta]], ind_weight)); ; 
% sparse_mult = sparse([[1:nsta];[1:nsta]], best_ind_flce', ind_weight); 
% pdf_terp * sparse_mult; 
% timeit(fhandtest); 

% mm_terp(best_ind_fl) - vsta_mod % Test to make sure getting reasonable
% velocities. 

% fhandtest = @()floor(best_ind); 
% timeit(fhandtest)

% mmi = 

% Plot to check that we interpolated to a common mm correctly. 
figure(12); clf; hold on; 
tiledlayout(2,1,'TileSpacing','compact'); 
nexttile(); hold on; title('Interpolated') 
for ista = 1:nsta
    plot(ista  + pdf_terp(:,ista), mm_terp ); 
end

pdf_xplot = [1:nsta]' + pdf_mod; 
scatter(pdf_xplot, vsta_mod); 

nexttile(); hold on; title('Original'); 
for ista = 1:nsta; 
    plot(ista + pdfs(ista).pm, pdfs(ista).mm)
end
% scatter([1:length(pdf_mod)]' + pdf_mod, vsta_mod); 


%% Efficient interpolation from grid to stations. 
% Goal: Make a matrix which can multiply by (flattened array of) station
% valuaes, and immediately give back the interpolated station values. All
% distance calculations and such are done only once on the front end. 

% Determine which nodes are closest to each station and find their weights.
% gridi = zeros(size(xgrid)); 
ngrid = (size(xgrid,1)*size(xgrid,2)); 
gridi_r = [1:ngrid]'; 
vgrid_r = vgrid(gridi_r); 
gridi = reshape(gridi_r, size(xgrid,1), size(xgrid,2)); % Map between 2d grid and 1d. 
ngrid = length(gridi_r); 

nearesti = zeros(nsta, 4); % Don't think I need this? 

grid_terp = zeros(ngrid, nsta); % Maybe make sparse if needed

for ista = 1:nsta; 
    staxi = stax(ista); 
    stayi = stay(ista); 

    to_east  = xgrid >= staxi; 
    to_west  = xgrid <  staxi; 
    to_north = ygrid >= stayi; 
    to_south = ygrid <  stayi;

    xdist = xgrid - staxi;  
    ydist = ygrid - stayi; 
    tdist = sqrt(xdist.^2 + ydist.^2); 

    [easti,  ~] = find(to_east ); 
    [westi,  ~] = find(to_west ); 
    [~, northi] = find(to_north); 
    [~, southi] = find(to_south); 
    easti = min(easti); 
    westi = max(westi); 
    northi = min(northi); 
    southi = max(southi); 

    nei = gridi(easti, northi);
    nwi = gridi(westi, northi); 
    sei = gridi(easti, southi); 
    swi = gridi(westi, southi); 

    box_map = [nei, nwi, sei, swi]; 

    nearesti(ista, :) = box_map; 

    dist_to_crnr = tdist(box_map); 
    weight_to_crnr = 1./dist_to_crnr; % TODO think about if this is the best interpolation method for this. 
    weight_to_crnr = weight_to_crnr ./ sum(weight_to_crnr); 
    
    if any(weight_to_crnr) == inf; 
        error('Should handle station being right on a node'); 
    end

%     grid_terp(box_map,ista) = 1 ;
    grid_terp(box_map,ista) = weight_to_crnr; 

% % %     ne = to_north .* to_east; 
% % %     nw = to_north .* to_west; 
% % %     se = to_south .* to_east; 
% % %     sw = to_south .* to_west; 

%     xgrid(ne)
end
if any(isnan(grid_terp)); 
    error('There are nans in grid_terp. Figure out a solution')
end

% Now can matrix multiply vgrid_r' * weight_to_crnr to interpolate velocity at
% each station. Do not need to worry about distances again. 
% I thnk also gridi_r needs to be passed to flatten vgrid
vsta_terp_mat = (vgrid_r' * grid_terp)'; 
% fhand_matmult = @()(vgrid_r' * grid_terp)';
% fhand_matmult = @()(vgrid(gridi_r)' * grid_terp)';

%%% Now I just pass vgrid_r as argument and grid_terp as always present
%%% argument



%%
fhand=@(vgrid)a3_1_penalty_efficient(vgrid,...
    pdfs, rough_scale, dx2, dy2, xgrid, ygrid, stax, stay); % Supply constants to the function handle. 
fhand(vgrid); 
% vgrid_start = vgrid; 

fhand2=@()fhand(vgrid); 

% % % %%
% % % for iii=1:100; 
% % %     fhand2(); 
% % % end