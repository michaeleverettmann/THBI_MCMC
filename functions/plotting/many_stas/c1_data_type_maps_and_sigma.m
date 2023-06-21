% Make a map of where each data type is available. 
% Also get the average sigma for each data type. 

clc; clear; restoredefaultpath; 
run('a0_parameters_setup.m'); % !!! Set up all parameters and such in a0. Because there may be many scripts here dependent on those parameters. 
addpath('/Users/brennanbrunsvik/Documents/repositories/Base_code/colormaps/redblue'); 
tectpath1 = "/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/map_view/whitmeyer_karlstrom/layers/"; % Path to what is currently whitmeyer and karlstrom dataset. 
f_distance_pt_to_sta = './_distance_pt_to_sta.mat';
f_xsect_positions = './xsect_positions.mat'; 
compiled_results = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/compiled_results_standard.mat'; % SHould be generated from script: /Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/a1_collate_sta_models.m
pt_dist_nan = 100 /(6371 * 2 * pi / 360); % Don't plot if no station within this many km


% define colors. 
color_front    = [0,0,0]./255; % Grenville Front, Appalachian front.
color_thrust   = [125, 125, 125]./255;
color_mag      = [190, 000, 000]./255;
color_grav     = [000, 150, 000]./255;
color_ha       = [255, 000, 010]./255; % Harrisonburg Anomaly
color_rift     = [255, 111, 000]./255;
color_text     = 0.*[1,1,1];

% Define any colors. 
clr_tectfiles = struct("grv_frt", color_front, "MCR", color_rift, ...
    "Reelfoot", color_rift, "something_province", color_thrust); % Field names have to precicely correspond to file names from this tectonic dataset. 

%%% Cross section stuff, copy to cross-section code. 
version_surf = 7; 
ll_min_max_map = [-89  -72   32   46]; % Map view
xsect_positions = load(f_xsect_positions); 
lolim = xsect_positions.lolim; 
lalim = xsect_positions.lalim; 


n_contour = 30; 


depths = [25, 95, 145]; % Try loading these depths. Probably need to type manually for now, but could save as a .mat file in future. 
parms_other = ["zsed", "zmoh"]; 


%% Figure out which data types are available at each station. 
each_type = ["sig_SW_Ray_phV_eks",...
    "sig_SW_Ray_phV_dal",...
    "sig_SW_Ray_phV_lyneqhelm",...
    "sig_SW_Ray_phV_lynant",...
    "sig_SW_Lov_phV",...
    "sig_RF_Sp_ccp",...
    "sig_HKstack_P",...
    "sig_SW_HV"]; 


mdls = load(compiled_results).mdls; 

nsta = length(mdls.lon); 
ntypes = length(each_type); 


hastype    = logical(zeros(nsta, ntypes)); 
sigma_mu   = nan(nsta, ntypes); 
sigma_std  = nan(nsta, ntypes); 


for ista = 1:nsta;      
    md = mdls.model{ista}; 
    hyper = md.hyperparms;
    fnames = fieldnames(hyper); 
    types = string(fnames); 
    for itype = 1:ntypes; 
        hastype_i = any(each_type(itype) == types); 
        hastype(ista,itype) = hastype_i;  

        %%% Check for a datatype that wasn't manually added to each_type. Print to screen if something is found. 
        for itype2 = 1:length(types); 
            if ~ any( each_type == types(itype2) ); 
                fprintf('New type: %s\n', types(itype2))
            end
        end
        %%% End check. 
        

        %%% Add sigma value, if we have it. 
        if hastype_i; 
            sigma_mu (ista, itype) = hyper.(each_type(itype)).mu_log10 ; 
            sigma_std(ista, itype) = hyper.(each_type(itype)).std_log10; 
        end
        %%%
    end
end

perc_data_mis = 100-100*sum(hastype,'all')/(size(hastype,1)*size(hastype,2)); % Missing data percent
fprintf('Amount of data that was unavailable: %1.2f%%\n', perc_data_mis)

%% Check sigma inverted for each data type. 
sigma_mu_10 = 10 .^ sigma_mu; 
sigma_mean = nanmean(sigma_mu_10, 1); 
sigma_mean = log10(sigma_mean); 
n_stations_available = sum(~isnan(sigma_mu_10)); 
sigma_txt = "From script: /Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/c1_data_type_maps_and_sigma.m. \n"; 
for itype = 1:length(each_type); 
    sigma_txt = sigma_txt + sprintf("%s: sigma = %1.3f. N stations = %1.0f\n",...
        each_type(itype), sigma_mean(itype), n_stations_available(itype)); 
end
fid = fopen('check_each_data_type/average_sigma.txt', 'wt');
fprintf(fid, sigma_txt);
fclose(fid);
%%

% sfsmat = load('surface_out_example.mat'); xgrid = sfsmat.xgrid; ygrid = sfsmat.ygrid; llminmax = sfsmat.llminmax; latgrid = sfsmat.latgrid; longrid = sfsmat.longrid; 
% mdls = load(fresults).mdls; % For sta lon and lat so we know where to and not to plot. 


%% Any borders to plot
%brb2023.02.21 These were copied from the 2021 ENAM paper. They were made only roughly, so I need more accurate shape files. 
app_bord = [-73.9819443, -74.619146 , -75.520073 , -76.1132732, -76.7944389,-77.1789631, -77.6294266, -77.9920049, -78.2117582, -78.4383114,-79.1416076, -79.8888598, -80.2331436, -80.7605221, -81.2438831,-81.9140786, -82.4853569, -83.2434644, -84.3284171, -84.7700317,-85.1106908, -85.5944904, -86.0889287, -87.330563; 40.4970924,  40.7306085,  40.9467137,  41.1455697,  41.2695495,41.2365112,  41.0959121,  40.8719878,  40.7056279,  40.3967643,39.5548831,  38.4965935,  38.1259146,  37.7272803,  37.5445773,37.3439591,  37.1603165,  36.8268747,  36.1733569,  35.880149 ,35.4427709,  34.7777158,  34.1799976,  32.9349287]; 
app_bord2 = [-73.9819  -74.1099  -73.8461  -73.3406  -73.0769  -72.4395; 40.4971   40.9965   42.4397   43.8504   44.8247   45.6140];
gre_bord = [-82.8790832, -83.5164453, -83.6483134, -84.0878735, -84.3516096, -85.714246 , -87.3406184, -88.1098487; 43.9928145,  41.8040781,  40.8969058,  39.6733704,  38.9764925, 36.6155276,  34.9940038,  34.3978449]; 
app_bord = [flipud(app_bord')', app_bord2]; 

%%

ll_min_max_map = [-89  -68   26   46]; % Map view
figure(17); clf; hold on; set(gcf,'pos', [87 856 692 476]); 

m_proj('mercator', 'long',[ll_min_max_map(1), ll_min_max_map(2)],...
                   'lat',[ll_min_max_map(3), ll_min_max_map(4)]);

% State lines
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0.5*[1 1 1])
end

% Coast lines
% cst = m_coast('patch',[1 1 1], 'FaceAlpha', 0); 

% Cross section position, labels
xsect_letters = ["A", "B", "C", "D", "E", "F", "G", "H"]; 
for ixsect = 1:size(lolim,1); 
    Q1 = [lalim(ixsect, 1), lolim(ixsect, 1)];
    Q2 = [lalim(ixsect, 2), lolim(ixsect, 2)]; 
    [profd,profaz] = distance(Q1(1),Q1(2),Q2(1),Q2(2));
    gcarc = linspace(0, profd, 100)'; 
    d2km = 2 * pi * 6371 / 360; % Yes I know this is 111, but might as well be somewhat precise :) 
    dist_arc = d2km * gcarc; % km 
    [lat_surf_line, lon_surf_line] = reckon(Q1(1), Q1(2), gcarc, profaz); 
    m_plot(lon_surf_line, lat_surf_line, 'k', 'linewidth', 1); 

    line_shift = 0; 
    color_xsectlab = .8.*[1,1,1]; 
    m_text(lon_surf_line(1+line_shift )-0.4, lat_surf_line(1+line_shift  )+0.3, xsect_letters(ixsect), ...
        'color', color_text, 'fontweight', 'normal', 'units', 'data', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); 
    m_text(lon_surf_line(end-line_shift)+0.4, lat_surf_line(end-line_shift) -0.3 , xsect_letters(ixsect)+"'", ...
        'color', color_text, 'fontweight', 'normal', 'units', 'data', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); 

end

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

% % Label features. 
% m_text(-83, 39, 'GF', 'color', color_text, 'fontsize', 12, ...
%     'rotation', 75, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontweight', 'bold'); 
% m_text(-83, 36, 'AF', 'color', color_text, 'fontsize', 12, ...
%     'rotation', 45, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontweight', 'bold'); 

% % Label the plot
% m_text(-74, 30.5, '?????', 'Units','data', 'HorizontalAlignment', 'center',...
%     'VerticalAlignment', 'bottom'); 

% Grid, coastlines I think, etc. Background colors. 
m_grid('box','fancy','linestyle','none','gridcolor',.5 .*[1,1,1],...
    'backcolor','none');

% Now scatter stations. 
% dtype_colors = jet(ntypes); 
type_txt = ["Rayleigh Ekstrom", ...
    "Rayleigh Dalton", ...
    "Rayleigh Lynner Ambient", ...
    "Rayleigh Lynner Helmholtz", ...
    "Love Ekstrom", ...
    "Sp Hopper", ...
    "HK EARS", ...
    "HV Shen"]; % String for legend

dtype_colors = [1,0,0;0,1,0;0,0,1]; 

% cla; disp('clearing axis')
hold on; 
mkr_size = 15; 


% has_all_types = 
prct_stas_have = sum(hastype,1) / size(hastype, 1); 
b_miss_type = prct_stas_have~=1; 
type_plt_inds = find(b_miss_type); 

str_always_present = ("100%%: " + type_txt(~b_miss_type) + ",\n");
str_always_present = char(str_always_present.join('') );
str_always_present = str_always_present(1:end-3); % Get rid of last \n
str_always_present = ['\n' str_always_present '.']
str_always_present = sprintf(str_always_present); 
 
h_always_present = m_scatter(mdls.lon, mdls.lat, 13, 'k', 'filled', 'DisplayName', ...
    str_always_present); 

hnd_leg = []; 

for itype = type_plt_inds; 
    pltsta = hastype(:,itype); % Boolean, plot these stations. We have this data there. 
    slon = mdls.lon(pltsta); 
    slat = mdls.lat(pltsta); 
%     hsct = m_scatter(slon, slat, mkr_size); 

    % Marker orientations.
    ang = (find(itype==type_plt_inds)-1) / (length(type_plt_inds)) * 360; 
    dline = .4; 
    slond = slon + dline * cosd(ang); 
    slatd = slat + dline * sind(ang); 

    % Marker colors. 
    mkr_color = dtype_colors(find(itype==type_plt_inds),:); 

    % Legend text
    prc_avail = sum(hastype(:,itype))/size(hastype,1)*100; 
    txt_leg = sprintf('%3.0f%%, %s',prc_avail, type_txt(itype)); 


    hplt=m_plot([slon, slond]', [slat, slatd]', 'Color', mkr_color, ...
        'linewidth', 2, 'DisplayName', txt_leg); 

    hnd_leg(end+1) = hplt(1); % Only plot legend for one station
end

hnd_leg(end+1) = h_always_present; 

lgd = legend(hnd_leg, 'Location', 'southeast'); 


exportgraphics(gcf, sprintf('check_each_data_type/sta_map.pdf')); 