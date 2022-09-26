clc; clear; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
mdls = load(fresults).mdls; 

%% % Reorganize structure. Could be more - or less - useful one way or other. 
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

%%% Map-view style parameters
mapFigNum = 1; 
figure(mapFigNum); clf; hold on; set(gcf, 'color', 'white', 'pos', [-1152 297 1024 735]); 
tiledlayout(2,3,'TileSpacing','compact'); 
fun_mapbackground(mdls.lon, mdls.lat, zmoh, 'Moho depth (km)'); 
fun_mapbackground(mdls.lon, mdls.lat, zsed, 'Sed depth (km)'); 
fun_mapbackground(mdls.lon, mdls.lat, xicr, '\xi'); 
fun_mapbackground(mdls.lon, mdls.lat, zmohsig, 'Moho error (km)'); 
fun_mapbackground(mdls.lon, mdls.lat, zsedsig, 'Sed error (km)'); 
fun_mapbackground(mdls.lon, mdls.lat, vpvs, 'Vp/Vs'); 
exportgraphics(gcf, sprintf('%smap_view_mod.pdf', figPath))

%%% Pseudo-tomography
mapFigNum = 2; 
figure(mapFigNum); clf; hold on; set(gcf, 'color', 'white', 'pos', [-1152 297 1024 735]); 
tiledlayout(2,3,'TileSpacing','compact'); 
for ivs=1:length(z_vs); 
    fun_mapbackground(mdls.lon, mdls.lat, vs(:,ivs), sprintf('%1.1f km',z_vs(ivs))); 
end
exportgraphics(gcf, sprintf('%spseudo_tomography.pdf', figPath))
