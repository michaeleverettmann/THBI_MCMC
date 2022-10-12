% % % clc; clear; 
% % % run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
% % % fpdfs    = sprintf('%scompiled_pdfs_%s.mat',out_dir,STAMP); % File with pdfs. a2_1...m
% % % 
% % % mdls = load(fresults).mdls; 
% % % 
% % % %% parameters. 
% % % rough_scale = 5e-9; % How much to penalize roughness.
% % % max_inv_iterations = 300; % How many iterations to allow in inversion. 
% % % z_vs = [0.5, 15, 50, 100, 150, 200]; % Depths where we want to solve for a velocity model. 
% % % ivel = 4; % Temporary. Expand to a loop. 
% % % 
% % % 
% % % 
% % % %% Reorganize structure. Could be more - or less - useful one way or other. 
% % % %%% BORROWED FROM A2
% % % models = struct(); 
% % % zmoh          = zeros(length(mdls.nwk),1); 
% % % zmohsig       = zeros(length(mdls.nwk),1);
% % % zsed          = zeros(length(mdls.nwk),1);
% % % zsedsig       = zeros(length(mdls.nwk),1);
% % % xicr          = zeros(length(mdls.nwk),1); 
% % % % xicrboundsmaybe= zeros(length(mdls.nwk),2);
% % % vpvs          = zeros(length(mdls.nwk),1); 
% % % % vpvssig       = zeros(length(mdls.nwk),1);
% % % 
% % % 
% % % vs            = zeros(length(mdls.nwk),length(z_vs)); 
% % % 
% % % for imd = 1:length(mdls.nwk); 
% % % 
% % %     % Make alternative format models variable. 
% % %     models(imd).nwk   = mdls.nwk  {imd}; 
% % %     models(imd).sta   = mdls.sta  {imd}; 
% % %     models(imd).lat   = mdls.lat  (imd); 
% % %     models(imd).lon   = mdls.lon  (imd); 
% % %     models(imd).model = mdls.model{imd}; 
% % %     models(imd).dir   = mdls.dir  {imd}; 
% % %     
% % %     % Extract map-view inversion results in cases where there's just one paramete for a station. 
% % %     mdl = mdls.model{imd}; 
% % %     zmoh              (imd) = mdl.zmohav; 
% % %     zmohsig           (imd) = mdl.zmohsig; 
% % %     zsed              (imd) = mdl.zsedav;     
% % %     zsedsig           (imd) = mdl.zsedsig; 
% % %     xicr              (imd) = mdl.xicrav;     
% % % %     xicrboundsmaybe   (imd,:) = mdl.xicrsig2'; 
% % %     vpvs              (imd) = mdl.vpvsav; 
% % % 
% % %     % Pseudo-tomography
% % %     dist_from_z = (abs(z_vs - mdl.Z)); 
% % %     [C,I ] = min(dist_from_z); 
% % %     vs(imd,:) = mdl.VSav(I'); % Velocity at depths. 
% % % end
% % % 
% % % 
% % % %% Setting up coordinate system. 
% % % lonmin = min(mdls.lon); 
% % % latmin = min(mdls.lat); 
% % % lonmax = max(mdls.lon); 
% % % latmax = max(mdls.lat); 
% % % 
% % % llminmax = [lonmin, lonmax, latmin, latmax]; % For easy passing to plotting functions
% % % 
% % % figure(1); clf; hold on; 
% % % m_proj('lambert', 'long',[lonmin, lonmax],'lat',[latmin, latmax]);
% % % [stax, stay] = m_ll2xy(mdls.lon, mdls.lat);  % Get stax and stay!!!
% % % 
% % % %% Plot stations positions. 
% % % a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, ...
% % %     'fignum', 1, 'title', 'station positions')
% % % 
% % % %% Setting up surface. 
% % % 
% % % nodes_per_degree = 4; 
% % % nx = ceil((lonmax - lonmin) * nodes_per_degree); 
% % % ny = ceil((latmax - latmin) * nodes_per_degree);  
% % % edge_space = 0.2; % How far to go beyond max and min stax and y, in fraction. 
% % % 
% % % 
% % % DX = max(stax) - min(stax); 
% % % DY = max(stay) - min(stay); 
% % % xline = linspace(min(stax) - edge_space * DX, max(stax) + edge_space * DX, nx)'; 
% % % yline = linspace(min(stay) - edge_space * DY, max(stay) + edge_space * DY, ny)'; 
% % % [xgrid, ygrid] = ndgrid(xline, yline);  
% % % [longrid, latgrid] = m_xy2ll(xgrid, ygrid); 
% % % 
% % % vs_interp = griddata(stax, stay, vs(:,ivel), xgrid, ygrid, 'cubic'); 
% % % % vs_interp(isnan(vs_interp)) = nanmean(vs_interp, 'all'); 
% % % 
% % % % Plot surface. 
% % % a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', vs(:,ivel),...
% % %     'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', vs_interp, ...
% % %     'fignum', 2, 'title', 'Simple v interpolation');  
% % % scatter(xgrid, ygrid, 5, 'k', 'filled')
% % % exportgraphics(gcf, 'surface_simple_interpolation.pdf'); 
% % % 
% % % 
% % % 
% % % %% Get PDF for stations. 
% % % % Use an example pdf while developing surface inversion stuff. 
% % % % fpdf = '~/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/pdf_example.mat'; 
% % % % pdf = load(fpdf).pdf_example; 
% % % % Probably something like pdf(ista) = load(pdf_sta). 
% % % pdf_file = load(fpdfs); 
% % % pdfs = pdf_file.pdfs; 
% % % 
% % % 
% % % % Make a starting model. 
% % % vgrid = vs_interp; % Probably a fine starting model. 
% % % vgrid(isnan(vgrid)) = nanmean(nanmean(vgrid)); 
% % % % vgrid = vgrid + randn(size(vgrid))*0.5; % Give some error
% % % 
% % % % Plot the starting model. 
% % % a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', vs(:,ivel),...
% % %     'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', vgrid, ...
% % %     'fignum', 3, 'title', 'V starting model'); 
% % % 
% % % %% Set up grid roughness calculations. Use 2nd derivative smoothing. 
% % % dx2 = ((xgrid(1:end-2,:) - xgrid(3:end,:))/2).^2; % Squared, so sign doesn't matter. 
% % % dy2 = ((ygrid(:,1:end-2) - ygrid(:,3:end))/2).^2; 
% % % x_dx2 = xgrid(2:end-1,:); y_dx2 = ygrid(2:end-1,:); % x and y positions where we have dx. I think just for plotting. 
% % % y_dy2 = ygrid(:,2:end-1); x_dy2 = xgrid(:,2:end-1); % y and x positions where we have dy
% % % % Note for optimizing: Use a scalar (not matrix) for dx and dy. But as is, it's more versatile for map projections. 
% % % 
% % % % Example roughness to make sure the calculations are good. 
% % % % This is the calculation to get roughness throughout the inversion
% % % % (Basically). 
% % % dvdx2 = (vgrid(1:end-2,:) - 2*vgrid(2:end-1,:) + vgrid(3:end,:))./dx2; 
% % % dvdy2 = (vgrid(:,1:end-2) - 2*vgrid(1,2:end-1) + vgrid(:,3:end))./dy2; 
% % % dvdx2 = dvdx2.^2; % Square, mostly to prevent negative values...
% % % dvdy2 = dvdy2.^2; 
% % % roughness = sum(dvdx2, 'all') + sum(dvdy2, 'all'); 
% % % roughness = roughness * rough_scale; % How to get roughness penalty value. 
% % % 
% % % a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay,...
% % %     'xgrid', x_dx2, 'ygrid', y_dx2, 'vgrid', dvdx2, ...
% % %     'fignum', 4, 'title', 'X roughness'); colorbar(); 
% % % a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay,...
% % %     'xgrid', x_dy2, 'ygrid', y_dy2, 'vgrid', dvdy2, ...
% % %     'fignum', 5, 'title', 'Y roughness'); colorbar(); 
% % % 
% % % % % % % % How to get "modelled" velocity throughout inversion. 
% % % % % % % vsta = interpn(xgrid, ygrid, vgrid, stax, stay, 'linear'); % MODELLED velocity at each station. 
% % % % % % % 
% % % % % % % % Probability associated with a velocity. 
% % % % % % % pv_mod = zeros(size(stax)); 
% % % % % % % for ista=1:length(pv_mod); 
% % % % % % %     pv_mod(ista) = interp1(pdfs(ista).mm, pdfs(ista).pm, vsta(ista), 'linear', 0 ); % Probability of a velocity at a specific station, from our velocity model surface. 
% % % % % % % end
% % % % % % % 
% % % % % % % penalty = - sum(pv_mod); % Lower penalty is higher probability. 
% % % % % % % penalty = penalty + roughness; 
% % % % % % % 
% % % % % % % fhand=@(vgrid)a3_1_penalty(vgrid,...
% % % % % % %     pdfs, rough_scale, dx2, dy2, xgrid, ygrid, stax, stay); % Supply constants to the function handle. 
% % % % % % % fhand(vgrid); 
% % % 
% % % vgrid_start = vgrid; 
% % % 
% % % %% Prep for efficient inverison. Some constant variables. 
% % % 
% % % % Interpolate each pdf to common mm grid. 
% % % % Important for computational efficient when calculating penalty. 
% % % nsta = length(pdfs); 
% % % nmm = 300; 
% % % [pdf_terp, mm_terp, dmm_di] = p_prep_mm_to_pdf(pdfs, nmm); 
% % % 
% % % % Interpolate mm from the grid to stations.  
% % % [fhand_vec, fhand_mat, grid_terp, nearesti, weighti...
% % %     ] = p_prep_grid_to_sta_interp(...
% % %     xgrid, ygrid, vgrid, stax, stay); 
% % % 
% % % % The thing we want to minimize. Make a function handle. 
% % % fhand_penalty=@(vgrid)a3_1_penalty_efficient(vgrid,...
% % %     pdf_terp, rough_scale, dx2, dy2, xgrid, ygrid, stax, stay, ...
% % %     nearesti, weighti, min(mm_terp), dmm_di, nmm, nsta); 
% % % 
% % % %% Test the above efficient inversion interpolation stuff. 
% % % vsta_mod = linspace(min(mm_terp), max(mm_terp), nsta)'; % FOR TESTING Easy values, for testing. 
% % % pdf_mod = p_mm_to_pdf_dmdi(...
% % %     vsta_mod, pdf_terp, min(mm_terp), dmm_di, nmm, nsta); 
% % % 
% % % % Plot to check that we interpolated to a common mm correctly. 
% % % figure(12); clf; hold on; 
% % % set(gcf, 'color', 'white')
% % % tiledlayout(2,1,'TileSpacing','compact'); 
% % % nexttile(); hold on; box on; set(gca, 'LineWidth', 1.5); 
% % % title('Interpolated') 
% % % for ista = 1:nsta
% % %     plot(ista  + pdf_terp(:,ista), mm_terp ); 
% % % end
% % % scatter([1:nsta]' + pdf_mod, vsta_mod); 
% % % 
% % % nexttile(); hold on; box on; set(gca, 'LineWidth', 1.5); 
% % % title('Original'); 
% % % for ista = 1:nsta; 
% % %     plot(ista + pdfs(ista).pm, pdfs(ista).mm)
% % % end
% % % 
% % % %% Example of running the inversion. 
% % % options = optimoptions("fminunc",Display="iter",...
% % %     MaxFunctionEvaluations=inf,MaxIterations=max_inv_iterations,...
% % %     Algorithm='quasi-newton');
% % % 
% % % opts = options;
% % % opts.Algorithm = 'quasi-newton';
% % % opts.HessianApproximation = 'lbfgs';
% % % opts.SpecifyObjectiveGradient = false;
% % % % overallTime = tic;
% % % 
% % % fprintf('Sum of max of pdf of each station: %1.1f.\nExcluding roughness, this is best possible penalty.\n', sum(max(pdf_terp)) )
% % % 
% % % tic; 
% % % [vgrid_out,fval_out,flag_out,output_out] =...
% % %     fminunc(fhand_penalty, vgrid_start, opts);
% % % toc 
% % % 
% % % 
% % % 
% % % vgrid_out(isnan(vs_interp)) = nan; 
% % % 
% % % % [vgrid_out,timetable.Fval("LBFGS_NoGrad"),timetable.Eflag("LBFGS_NoGrad"),output] =...
% % % %     fminunc(fhand, vgrid_start, opts);
% % % % 
% % % % 
% % % % timetable.Time("LBFGS_NoGrad") = toc(overallTime);
% % % % timetable.Iters("LBFGS_NoGrad") = output.iterations;
% % % % timetable.TimePerIter("LBFGS_NoGrad") =...
% % % %     timetable.Time("LBFGS_NoGrad")/timetable.Iters("LBFGS_NoGrad");
% % % 
% % % 
% % % % % % %% Example of running the inversion. Different inversion approach. No derivative. 
% % % % % % options = optimset(Display="iter",...
% % % % % %     MaxIter=100,PlotFcns=@optimplotfval, 'algorithm', );
% % % % % % 
% % % % % % fprintf('Sum of max of pdf of each station: %1.1f.\nExcluding roughness, this is best possible penalty.\n', sum(max(pdf_terp)) )
% % % % % % 
% % % % % % tic; 
% % % % % % [vgrid_out,fval_out,flag_out,output_out] =...
% % % % % %     fminsearch(fhand_penalty, vgrid_start, options);
% % % % % % toc 
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % vgrid_out(isnan(vs_interp)) = nan; 
% % % % % % 
% % % % % % % [vgrid_out,timetable.Fval("LBFGS_NoGrad"),timetable.Eflag("LBFGS_NoGrad"),output] =...
% % % % % % %     fminunc(fhand, vgrid_start, opts);
% % % % % % % 
% % % % % % % 
% % % % % % % timetable.Time("LBFGS_NoGrad") = toc(overallTime);
% % % % % % % timetable.Iters("LBFGS_NoGrad") = output.iterations;
% % % % % % % timetable.TimePerIter("LBFGS_NoGrad") =...
% % % % % % %     timetable.Time("LBFGS_NoGrad")/timetable.Iters("LBFGS_NoGrad");
% % % 
% % % %% Plot inversion output. 
% % % a3_2_plot_surface_simple(llminmax, 'stax', stax, 'stay', stay, 'stav', vs(:,ivel),...
% % %     'xgrid', xgrid, 'ygrid', ygrid, 'vgrid', vgrid_out, ...
% % %     'fignum', 6, 'title', 'V output'); 
% % % exportgraphics(gcf, './surface_inversion_rough.pdf'); 
% % % 
% % % % Temporary. Mostly for testing. 
% % % save('surface_out_example.mat', 'longrid', 'latgrid',...
% % %     'xgrid', 'ygrid', 'vgrid_out', 'llminmax'); 
% % % 
% % % 
