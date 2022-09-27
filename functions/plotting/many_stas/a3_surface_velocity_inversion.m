clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
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


%%
lonmin = min(mdls.lon); 
latmin = min(mdls.lat); 
lonmax = max(mdls.lon); 
latmax = max(mdls.lat); 

figure(1); clf; hold on; 

m_proj('lambert', 'long',[lonmin, lonmax],'lat',[latmin, latmax]);
m_coast('patch',[1 1 1]); % m_coast('patch',[1 .85 .7]);

[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end

m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[.3 .75 1]);

[stax, stay] = m_ll2xy(mdls.lon, mdls.lat);  

scatter(stax, stay); 

%%
ivel = 5; 

DX = max(stax) - min(stax); 
DY = max(stay) - min(stay); 
edge_space = 0.2; 
xline = linspace(min(stax) - edge_space * DX, max(stax) + edge_space * DX, 10)'; 
yline = linspace(min(stay) - edge_space * DY, max(stay) + edge_space * DY, 11)'; 
[xgrid, ygrid] = ndgrid(xline, yline);  

vs_interp = griddata(stax, stay, vs(:,ivel), xgrid, ygrid, 'cubic'); 

contourf(xgrid, ygrid, vs_interp, 15); 
scatter(stax, stay, 60, vs(:,ivel), 'filled' ); 

%% 

% Use an example pdf while developing surface inversion stuff. 
fpdf = '~/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/pdf_example.mat'; 
pdf = load(fpdf).pdf_example; 
rough_scale = 5e-11; 

vgrid = vs_interp; % Probably a fine starting model. 
vgrid(isnan(vgrid)) = nanmean(nanmean(vgrid)); 

vgrid = vgrid + randn(size(vgrid))-0.5; % Give some error

% dx = xgrid(2:end,:)-xgrid(1:end-1,:); 
% dy = ygrid(:,2:end)-ygrid(:,1:end-1); 
% dx2 = xgrid(1:end-2,:) - 2*xgrid(2:end-1,:) + xgrid(3:end,:);
% dy2 = xgrid(:,1:end-2) - 2*ygrid(1,2:end-1) + ygrid(:,3:end); 
dx2 = ((xgrid(1:end-2,:) - xgrid(3:end,:))/2).^2; 
dy2 = ((ygrid(:,1:end-2) - ygrid(:,3:end))/2).^2; 
x_dx2 = xgrid(2:end-1,:); y_dx2 = ygrid(2:end-1,:); 
y_dy2 = ygrid(:,2:end-1); x_dy2 = xgrid(:,2:end-1); 

% Options for optimizing. 
% Use one value for dx and dy. But as is, it's more versatile for map
% projections. 


dvdx2 = (vgrid(1:end-2,:) - 2*vgrid(2:end-1,:) + vgrid(3:end,:))./dx2; 
dvdy2 = (vgrid(:,1:end-2) - 2*vgrid(1,2:end-1) + vgrid(:,3:end))./dy2; 
% dvdx2(isnan(dvdx2)) = 2; 
% dvdy2(isnan(dvdy2)) = 2; 


% figure(2); clf; hold on; 
% contourf(x_dx2, y_dx2, dvdx2, 50); 
% colorbar(); 

roughness = sum(dvdx2 .^ 2, 'all') + sum(dvdy2 .^ 2, 'all'); 
roughness = roughness * rough_scale; 

vsta = interpn(xgrid, ygrid, vgrid, stax, stay, 'linear'); % Modelled velocity at each station. 

% ista = 1; 
pv_mod = zeros(size(stax)); 
for ista=1:length(pv_mod); 
    pv_mod(ista) = interp1(pdf.v, pdf.p, vsta(ista) ); % Probability of a velocity at a specific station, from our velocity model surface. 
end

penalty = - sum(pv_mod); % Lower penalty is higher probability. 
penalty = penalty + roughness; 

% pdf.p

% fun(vgrid)@(); 

global pdf rough_scale dx2 dy2 xgrid ygrid stax stay


a3_1_penalty(vgrid); 
%%
options = optimoptions("fminunc",Display="none",...
    MaxFunctionEvaluations=inf,MaxIterations=inf,...
    Display='iter');

opts = options;
opts.HessianApproximation = 'lbfgs';
opts.SpecifyObjectiveGradient = false;
overallTime = tic;

[~,timetable.Fval("LBFGS_NoGrad"),timetable.Eflag("LBFGS_NoGrad"),output] =...
    fminunc(@a3_1_penalty, vgrid, opts);


timetable.Time("LBFGS_NoGrad") = toc(overallTime);
timetable.Iters("LBFGS_NoGrad") = output.iterations;
timetable.TimePerIter("LBFGS_NoGrad") =...
    timetable.Time("LBFGS_NoGrad")/timetable.Iters("LBFGS_NoGrad");



