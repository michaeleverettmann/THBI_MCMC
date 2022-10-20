clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
fpdfs = sprintf('%scompiled_pdfs_%s.mat',out_dir,STAMP); % File with pdfs. a2_1...m

mdls = load(fresults).mdls; 

%% parameters. 
rough_scale = 10e-9; % How much to penalize roughness.
max_inv_iterations = 15; % How many iterations to allow in inversion. 
version_surf = 3; 

to_invert = {"zsed", "zmoh"}; % Which model parameters to run. Those come first because they can influence later inversions.  
for inum = int16([5, 15, 20, 25, 30, 35, 40, ...
        45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 120, 140, 170, 210, 250, 300]/5); % Which depths/incidices to run. 
    to_invert{end+1} = inum; 
end

% Total hack to not worry as much about having depth slices of vs, but other things not being a function of depth. 
for iinv = to_invert; %!%! Add strings to the list to handle other parameters. Make sure they are always first in list. 
iinv = iinv{1}; 

%^%^ Handle whether doing depth or other model parameter inversion
v_at_depth = ~ strcmp(class(iinv), class("A string") ); % Use velocity from a depth, or one of the other parameters like moho depth. If a string is provided, we assume we are not using velocity at depth but another model parameter. %!%! Utilize v_at_depth
if v_at_depth 
    param = z_vs(iinv); %!%! Only do if v_at_depth. %!%! change variable depth. 
    this_inversion = sprintf('vs%1.0f',param); % String name affiliated with figures and files. %!%! change variable depth. 
else
    param = iinv; 
    this_inversion = sprintf('%s',param); 
end

mkdir(this_inversion); 

fprintf('Running inversion %s\n', this_inversion)

%% Reorganize structure. Could be more - or less - useful one way or other. 
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

% vs            = zeros(length(mdls.nwk),length(z_vs)); 
m_simple      = zeros(length(mdls.nwk), 1); 

nsta = length(mdls.nwk); 

for imd = 1:nsta; 

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
    if v_at_depth
% % %         dist_from_z = (abs(z_vs - mdl.Z)); %!%! These lines only revelvant if v_at_depth
% % %         [C,I ] = min(dist_from_z); %!%! These lines only revelvant if v_at_depth
% % %         vs(imd,:) = mdl.VSav(I'); % Velocity at depths. %!%! These lines only revelvant if v_at_depth
% % %         %!%! vs is derived from models. Just make it m_simple. 
        [~,iclosest_z] = min( abs(z_vs(iinv) - mdl.Z) ); 
        m_simple(imd) = mdl.VSav(iclosest_z);
    else
        m_simple(imd) = mdl.(iinv+"av"); % This is better than doing the median later, since median didn't work for some model parameters. 
    end
    
end


%% Setting up coordinate system. 
lonmin = min(mdls.lon-1); 
latmin = min(mdls.lat-1); 
lonmax = max(mdls.lon+1); 
latmax = max(mdls.lat+1); 

llminmax = [lonmin, lonmax, latmin, latmax]; % For easy passing to plotting functions

figure(1); clf; hold on; 
m_proj('lambert', 'long',[lonmin, lonmax],'lat',[latmin, latmax]);
[stax, stay] = m_ll2xy(mdls.lon, mdls.lat);  % Get stax and stay!!!

%% Plot stations positions. 
a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, ...
    'fignum', 1, 'title', 'station positions')

%% Setting up surface. 


if v_at_depth; %!%! nice
    ivel = find(z_vs == param); 
    iz = ivel; 
else
    iz = 'no different depths!'; 
end

nodes_per_degree = 4; 
nx = ceil((lonmax - lonmin) * nodes_per_degree); 
ny = ceil((latmax - latmin) * nodes_per_degree);  
edge_space = 0.2; % How far to go beyond max and min stax and y, in fraction. 


DX = max(stax) - min(stax); 
DY = max(stay) - min(stay); 
xline = linspace(min(stax) - edge_space * DX, max(stax) + edge_space * DX, nx)'; 
yline = linspace(min(stay) - edge_space * DY, max(stay) + edge_space * DY, ny)'; 
[xgrid, ygrid] = ndgrid(xline, yline);  
[longrid, latgrid] = m_xy2ll(xgrid, ygrid); 

m_interp = griddata(stax, stay, m_simple, xgrid, ygrid, 'cubic'); %!%! change variable vs_interp. Replace all m_simple. 
% m_interp(isnan(m_interp)) = nanmean(m_interp, 'all'); 

% % % %%% Handle points outside station bounds. 
% % % dist_handle = 0.03 * min([DX, DY]); 
% % % v_dft = mean(m_simple); 
% % % for ipt = 1:(nx*ny)
% % %     distpt = sqrt( (stax - xgrid(ipt)).^2 + (stay - ygrid(ipt)).^2 ); % distance of this point to each station
% % %     if min(distpt) > dist_handle;
% % %         sta_wt = 1./distpt; 
% % %         sta_wt = sta_wt / sum(sta_wt); 
% % %         dist_nearest = min(distpt); 
% % % %         sta_to_dft_wt = 1./((min(distpt)-dist_handle)/dist_handle); % Station weighting relative to default weighting. Station weighting is 1 at dist_handle from closest station. Station weighting should decrease moving away from stations. 
% % %         dft_to_sta_wt = ((dist_nearest-dist_handle) / dist_handle); % Station weighting relative to default weighting. Station weighting is 1 at dist_handle from closest station. Station weighting should decrease moving away from stations. 
% % %         disp([sta_to_dft_wt, v_new])
% % %         v_new = (sta_to_dft_wt * sum(sta_wt .* m_simple)) + ...
% % %             ((1-sta_to_dft_wt) * v_dft); 
% % %         m_interp(ipt) = v_new; 
% % %     end
% % % end

%%% Plot surface. 
a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', m_simple,... %!%! Replace vs(:,ivel)
    'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', m_interp, ...
    'fignum', 2, 'title', 'Simple v interpolation');  
% scatter(xgrid, ygrid, 5, 'k', 'filled')

exportgraphics(gcf, [this_inversion '/surface_simple_interpolation.pdf']); 




%% Make a starting model. 
mgrid = m_interp; % Probably a fine starting model.  %!%! Replace vgrid. Replace vs_interp. 
mgrid(isnan(mgrid)) = nanmean(nanmean(mgrid)); %!%! replace vgrid
% mgrid = mgrid + randn(size(mgrid))*0.5; % Give some error
% 
% 
%%% Handle points outside station bounds. 
for ipt = 1:(nx*ny)
    distpt = sqrt( (stax - xgrid(ipt)).^2 + (stay - ygrid(ipt)).^2 ); % distance of this point to each station
    distpt = distpt ./ (sqrt(DX.^2 + DY.^2)); 
    [distpts, distpti] = sort(distpt); 
%     diff()
%     exponent = cumprod(distpt(3)) * 8; 
    exponent = 4; 
    max_ind = size(m_simple, 1); 
%     max_ind = 30; 
    distpts = distpts(1:max_ind);
    vspt = m_simple( distpti(1:max_ind) ); %!%! Replace vs. Replace vspt. 
    sta_wt = (1./distpts).^exponent; 
    sta_wt = sta_wt / sum(sta_wt); 
    v_new = sum(sta_wt .* vspt); %!%! replace v_new. Replace vspt. 
    mgrid(ipt) = v_new; %!%! replace vgrid. Replace v_new. 
end

% Plot the starting model. 
a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', m_simple,... %!%! Replace vs(:,ivel)
    'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', mgrid, ...
    'fignum', 3, 'title', 'V starting model'); 
exportgraphics(gcf, [this_inversion '/surface_starting_model.pdf']); 


%% Set up grid roughness calculations. Use 2nd derivative smoothing. 
dx2 = ((xgrid(1:end-2,:) - xgrid(3:end,:))/2).^2; % Squared, so sign doesn't matter. 
dy2 = ((ygrid(:,1:end-2) - ygrid(:,3:end))/2).^2; 
x_dx2 = xgrid(2:end-1,:); y_dx2 = ygrid(2:end-1,:); % x and y positions where we have dx. I think just for plotting. 
y_dy2 = ygrid(:,2:end-1); x_dy2 = xgrid(:,2:end-1); % y and x positions where we have dy
% Note for optimizing: Use a scalar (not matrix) for dx and dy. But as is, it's more versatile for map projections. 

% Example roughness to make sure the calculations are good. 
% This is the calculation to get roughness throughout the inversion
% (Basically). 
dvdx2 = (mgrid(1:end-2,:) - 2*mgrid(2:end-1,:) + mgrid(3:end,:))./dx2; %!%! replace vgrid
dvdy2 = (mgrid(:,1:end-2) - 2*mgrid(1,2:end-1) + mgrid(:,3:end))./dy2; %!%! replace vgrid
dvdx2 = dvdx2.^2; % Square, mostly to prevent negative values...
dvdy2 = dvdy2.^2; 
roughness = sum(dvdx2, 'all') + sum(dvdy2, 'all'); % replace dvdx2 and dvdy2
roughness = roughness * rough_scale; % How to get roughness penalty value. 

a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay,...
    'xgrid', x_dx2, 'ygrid', y_dx2, 'vgrid', dvdx2, ...%Replace dvdx2 and dvdy2
    'fignum', 4, 'title', 'X roughness'); colorbar(); 
a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay,...
    'xgrid', x_dy2, 'ygrid', y_dy2, 'vgrid', dvdy2, ...%Replace dvdx2 and dvdy2
    'fignum', 5, 'title', 'Y roughness'); colorbar(); 

mgrid_start = mgrid; %!%! Replace vgrid_start and vgrid

%% Get PDF for stations. 
% Use an example pdf while developing surface inversion stuff. 
% fpdf = '~/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/pdf_example.mat'; 
% pdf = load(fpdf).pdf_example; 
% Probably something like pdf(ista) = load(pdf_sta). 
pdf_file = load(fpdfs); 
pdfs = pdf_file.pdfs_allparm; 


% If using different depths: 
% iz = 12; % Temporary. 
% param = pdf_file.pdfs_allparm(1).zatdep(iz); 
% fprintf('Temporarily using %1.0f depth for inversion. \n', param)

% pdfs = pdfs(:).vs{iz}; 

%%% Put this in function later. 
%!%! replace pdfs_vs. 
%!%! If not vs, access different part of pdfs. 
if v_at_depth
    pdfs_vs = pdfs(1).vs{1}; % Make a new structure (obnoxious). And have to start with the correct field names. Reason for new structure is that, I used a cell array for each different depth. Matlab doesn't actually access the nth stations ith cell array all in one call. 
    nsta = length(pdfs); 
    for ista = 1:nsta
        pdfs_vs(ista) = pdfs(ista).vs{iz}; %!%! Not access pdfs.vs
    end
    pdfs = pdfs_vs; %!%! replace pdfs_vs. 
else
    pdfs = [pdfs.(iinv)]; 
end
%%% Put this in function later. 

% this_inversion = sprintf('%s%s', 'vs', num2str(param) ); 


%% Prep for efficient inverison. Some constant variables. 

% Interpolate each pdf to common mm grid. 
% Important for computational efficient when calculating penalty. 
nsta = length(pdfs); 
nmm = 300; 
[pdf_terp, mm_terp, dmm_di] = p_prep_mm_to_pdf(pdfs, nmm); 

% Interpolate mm from the grid to stations.  
[fhand_vec, fhand_mat, grid_terp, nearesti, weighti...
    ] = p_prep_grid_to_sta_interp(...
    xgrid, ygrid, mgrid, stax, stay); %!%! replace vgrid

% The thing we want to minimize. Make a function handle. 
fhand_penalty=@(mgrid)a3_1_penalty_efficient(mgrid,...%!%! replace vgrid
    pdf_terp, rough_scale, dx2, dy2, xgrid, ygrid, stax, stay, ...
    nearesti, weighti, min(mm_terp), dmm_di, nmm, nsta); 

% Find the absolute max pdf sum that could be achieved
pdf_total_max = sum(max(pdf_terp,[],1)); 

%% Test the above efficient inversion interpolation stuff. 
msta_mod = linspace(min(mm_terp), max(mm_terp), nsta)'; % FOR TESTING Easy values, for testing. %!%! replace vsta_mod
pdf_mod = p_mm_to_pdf_dmdi(...
    msta_mod, pdf_terp, min(mm_terp), dmm_di, nmm, nsta); %!%! replace vsta_mod

% Plot to check that we interpolated to a common mm correctly. 
figure(12); clf; hold on; 
set(gcf, 'color', 'white')
tiledlayout(2,1,'TileSpacing','compact'); 
nexttile(); hold on; box on; set(gca, 'LineWidth', 1.5); 
title('Interpolated') 
for ista = 1:nsta
    plot(ista  + pdf_terp(:,ista), mm_terp ); 
end
scatter([1:nsta]' + pdf_mod, msta_mod); %!%! replace vsta_mod

nexttile(); hold on; box on; set(gca, 'LineWidth', 1.5); 
title('Original'); 
for ista = 1:nsta; 
    plot(ista + pdfs(ista).pm, pdfs(ista).mm)
end

%% Example of running the inversion. 
options = optimoptions("fminunc",Display="iter",...
    MaxFunctionEvaluations=inf,MaxIterations=max_inv_iterations,...
    Algorithm='quasi-newton');

opts = options;
opts.Algorithm = 'quasi-newton';
opts.HessianApproximation = 'lbfgs';
opts.SpecifyObjectiveGradient = false;
% overallTime = tic;

fprintf('Sum of max of pdf of each station: %1.1f.\nExcluding roughness, this is best possible penalty.\n', sum(max(pdf_terp)) )

% % % tic; 
% % % [vgrid_out,fval_out,flag_out,output_out] =...
% % %     fminunc(fhand_penalty, mgrid_start, opts);
% % % toc 


%% Iterate and see how model changes
tic; 
ii = 0; 
mgrid_temp = mgrid_start; %!%! Replace mgrid_temp and vgrid_start
opts_temp = opts; 
ii_vec = [0]; 
[pentot_ii,penpdf_ii,penprior_ii] = fhand_penalty(mgrid_temp); %!%! replace vgrid_temp
opts_temp.MaxIterations = 0; 
while ii < max_inv_iterations; 
%     opts_temp.MaxIterations = ceil((opts_temp.MaxIterations + 1)^(1.3)) ; 
%     possible_iters = ceil([opts_temp.MaxIterations, 10, opts_temp.MaxIterations - ii]); 
%     if ii < 5; possible_iters = [possible_iters, 1]; end
%     possible_iters(possible_iters < 1) = 1; 
%     opts_temp.MaxIterations = min(possible_iters); % Increase iteration gap as you go, but not above 10, and if we are at the end of the inversino, do enough iterations to just get to max iterations. 
    % How many iterations to do? Not efficient to start inversion too many times.     
    if ii < 5; 
        opts_temp.MaxIterations = 1; 
    elseif ii < 15; 
        opts_temp.MaxIterations = 2; 
    elseif ii < 50; 
        opts_temp.MaxIterations = 6; 
    elseif ii < 300; 
        opts_temp.MaxIterations = 15; 
    end
    if ii + opts_temp.MaxIterations > max_inv_iterations; 
        opts_temp.MaxIterations = min([1, max_inv_iterations - opts_temp.MaxIterations])
    end

    [mgrid_temp,fval_out,flag_out,output_out] =...%!%! replace vgrid_temp
        fminunc(fhand_penalty, mgrid_temp, opts_temp);%!%! replace vgrid_temp
    ii = ii + output_out.iterations; 
    ii_vec(end+1) = ii; 
    [pentot_ii(end+1),penpdf_ii(end+1),penprior_ii(end+1)] = fhand_penalty(mgrid_temp); %!%! replace vgrid_temp

end
mgrid_out = mgrid_temp; %!%! replace vgrid_temp, vgrid_out
toc

%% Plot how inversion progressed
% yscale_type = 'log'; 
for yscale_type = ["log", "linear"]; 


figure(13); clf; hold on; 
box on; set(gca, 'LineWidth', 1.5); 

title(sprintf('Surface inversion progress. Max possible: %1.1f', pdf_total_max),...
    'FontWeight','normal'); 
yyaxis('left'); 
ylabel('\Sigma P(m)')
set(gca, 'YScale', yscale_type, 'YDir', 'reverse')
grid on; 

plot(ii_vec, - penpdf_ii); 
scatter(ii_vec, - penpdf_ii, 'filled');  

% plot([min(ii_vec), max(ii_vec)], - pdf_total_max * ones(1,2) )

yyaxis('right'); 
ylabel('Roughness'); 
plot(ii_vec, penprior_ii); 
scatter(ii_vec, penprior_ii, 'filled'); 
set(gca, 'YScale', yscale_type)
xlabel('Iteration'); 
grid on; 


exportgraphics(gcf, sprintf('%s/surface_inversion_progress_%s.pdf',this_inversion, yscale_type)); 

end

%% Percent changes throughout inversion, per ii
fhand_prog = @(invval)(diff(invval) ./ diff(ii_vec))' ./ max(invval) * 100; 
[ii_vec(2:end)', fhand_prog(pentot_ii), fhand_prog(penpdf_ii), fhand_prog(penprior_ii)]



%%



% mgrid_out(isnan(m_interp)) = nan; 

% [mgrid_out,timetable.Fval("LBFGS_NoGrad"),timetable.Eflag("LBFGS_NoGrad"),output] =...
%     fminunc(fhand, mgrid_start, opts);
% 
% 
% timetable.Time("LBFGS_NoGrad") = toc(overallTime);
% timetable.Iters("LBFGS_NoGrad") = output.iterations;
% timetable.TimePerIter("LBFGS_NoGrad") =...
%     timetable.Time("LBFGS_NoGrad")/timetable.Iters("LBFGS_NoGrad");


% % % %% Example of running the inversion. Different inversion approach. No derivative. 
% % % options = optimset(Display="iter",...
% % %     MaxIter=100,PlotFcns=@optimplotfval, 'algorithm', );
% % % 
% % % fprintf('Sum of max of pdf of each station: %1.1f.\nExcluding roughness, this is best possible penalty.\n', sum(max(pdf_terp)) )
% % % 
% % % tic; 
% % % [mgrid_out,fval_out,flag_out,output_out] =...
% % %     fminsearch(fhand_penalty, mgrid_start, options);
% % % toc 
% % % 
% % % 
% % % 
% % % mgrid_out(isnan(m_interp)) = nan; 
% % % 
% % % % [mgrid_out,timetable.Fval("LBFGS_NoGrad"),timetable.Eflag("LBFGS_NoGrad"),output] =...
% % % %     fminunc(fhand, mgrid_start, opts);
% % % % 
% % % % 
% % % % timetable.Time("LBFGS_NoGrad") = toc(overallTime);
% % % % timetable.Iters("LBFGS_NoGrad") = output.iterations;
% % % % timetable.TimePerIter("LBFGS_NoGrad") =...
% % % %     timetable.Time("LBFGS_NoGrad")/timetable.Iters("LBFGS_NoGrad");

%% Plot inversion output. 
a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', m_simple,... %!%! replace vs(:,ivel)
    'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', mgrid_out, ... %!%! replace vgrid_out
    'fignum', 6, 'title', 'V output'); 
exportgraphics(gcf, sprintf('%s/surface_inversion_rough.pdf',this_inversion)); 

% Temporary. Mostly for testing. 
save('surface_out_example.mat', 'longrid', 'latgrid',...
    'xgrid', 'ygrid', 'llminmax'); % Longrid and stuff isn't going to change. 
save(sprintf('%s/surface_values_V%1.0f', this_inversion, version_surf), 'mgrid_out')


end