% This script mostly taken from the script /Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/collate_results/collate_results.m
recalculate = false; 
if ~recalculate; 
    fprintf('\n\nWarning: not recalculating PDF/Histogram through whole spaces.\n\n')
end



f_regions = 'choose_split_regions/tect_regions.mat'; % File with polygons
f_results = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_collate'; % Where are stations results stored? 
f_sta_list_all = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ENAM/batch/staList_all.txt'; 
f_stainf = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ENAM/INFO/stations.mat'; 

sta_dir = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_collate/%s_%s_dat1/standard/'; 
res_str = [sta_dir 'allmodels_perchain_orig.mat']; 
desired_chains = 12; 
desired_iter   = 16000; 

stainf = load(f_stainf); 
nt_sta = string(stainf.stainfo.nwk) + '_' + string(stainf.stainfo.stas); 
lon_sta = stainf.stainfo.slons; 
lat_sta = stainf.stainfo.slats; 

regions = load(f_regions);
regions_str = string(fieldnames(regions)); 
sta_list_all = string(table2cell(readtable(...
    f_sta_list_all, 'ReadVariableNames', 0))); 

sta_list_all = sta_list_all(:,1) + ' ' + sta_list_all(:,2); % Combine net and sta
sta_list_all = sta_list_all(1:end-1,:); % Remove non station final line


% have_results = nan(length(sta_list_all),1); 

% Medians of each velocity model. 
z_med = 0:1:300; % Interpolate to these velocities? 
v_med = nan(length(sta_list_all), length(z_med)); 

% For contouring 
vcont = 0.5:.025:5.5; 
mcont = zeros(length(z_med), length(vcont), length(sta_list_all)); 

% sta_list_all = sta_list_all(1:5); fprintf('Temp shorting stations\n')
if recalculate; 
    which_region = zeros(length(sta_list_all),1); 
    for istn = 1:length(sta_list_all); 
        fprintf('On station %1.0f\n', istn); 
        sta_str = sta_list_all(istn); 
        sta_str_splt = split(sta_str, ' '); 
        net = sta_str_splt(1); 
        sta = sta_str_splt(2); 
    
        file_res = sprintf(res_str, sta, net); 
        fpar = sprintf([sta_dir '%s'], sta, net, 'par.mat'); 
        fmod = sprintf([sta_dir '%s'], sta, net, 'final_model.mat'); 
        fmis = sprintf([sta_dir '%s'], sta, net, 'final_misfit.mat'); 
        fall = sprintf([sta_dir '%s'], sta, net, 'allmodels_perchain_orig.mat'); % These are big files! 
        fgd  = sprintf([sta_dir '%s'], sta, net, 'goodchains.mat'); 
    
        have_files = 0 ~= (exist(file_res, 'file') * ...
            exist(fpar, 'file') * exist(fmod, 'file') * exist(fmis, 'file') * ...
            exist(fall) * exist(fgd )); 
    
        if ~have_files; 
            have_results(istn) = false; 
            continue; 
        else
            have_results(istn) = true; 
        end
    
    
    
        % What is lat and lon of this sta? 
        ind_ll = find(strcmp(nt_sta, net + "_" + sta)); % Lat lon ind
        lat = lat_sta(ind_ll); 
        lon = lon_sta(ind_ll); 
    
        % Determine which polygon we are in. 
        which_poly = 0; % Not in any, by default
        for ipoly = 1:length(regions_str); 
            poly = regions.(regions_str(ipoly)); 
            is_in_poly = inpolygon(lon, lat, poly(:,1), poly(:,2)); 
            if is_in_poly; 
                which_poly = ipoly; 
            end
        end
        which_region(istn) = which_poly; 
    
        % Load and save velocity models. 
        model = load(fmod).final_model; 
        vs_terp = interp1(model.Z, model.VSav, z_med, 'linear'); 
        v_med(istn, :) = vs_terp; 
    
        % Get all the models. Contours. 
        suit = load(fall).allmodels_perchain; 
        gd   = load(fgd ).goodchains; 
        suit = suit(gd); % Only load good chains. 
    
        for ich = 1:length(suit); 
            mod = suit{ich}; 
            vs = mod.VS; 
            all_mods = zeros(length(z_med), length(mod)); 
            bad_mods = zeros(length(mod),1); 
    
            % Get chains that are good and interpolate them to common z
            for imod = 1:length(mod); 
                zimod = mod(imod).z; 
                vimod = mod(imod).VS; 
    
                if length(vimod) < 1; 
                    bad_mods(imod) = 1; 
                    continue
                end
                
                % Shift discontinuities. 
                z_dups = find(zimod(1:end-1) == zimod(2:end)); 
                if length(z_dups) > 0; 
                    zimod(z_dups) = zimod(z_dups) - 0.01; 
                    zimod(z_dups+1) = zimod(z_dups+1) + 0.01; 
                end
                all_mods(:,imod) = interp1(zimod, vimod, z_med ); 
            end
            all_mods = all_mods(:,bad_mods~=1);
            % Rough way to view result
    %         figure(1); clf; hold on; contourf(all_mods); colorbar; set(gca, 'ydir', 'reverse');   
    
            % Histogram of each chain 
            for iz = 1:length(z_med); 
                modz = all_mods(iz,:); 
                mcont(iz,:, istn) = mcont(iz,:, istn) + histc(modz, vcont); 
            end
    %         figure(2); clf; hold on; contourf(mcont(:,:,istn)); colorbar; set(gca, 'ydir', 'reverse');   
    
        end
    
    end
    
    have_results = logical(have_results); 
    
    save([pwd '/stas_collated.mat'], 'regions_str', 'z_med', 'which_region', 'v_med', 'mcont', 'vcont'); 
else
    %%
    clear; 
    load('./stas_collated.mat'); 
end

%% Calculate medians in each region
v_reg = zeros(length(regions_str), length(z_med)); 
for iregion = 1:length(regions_str); 
    v_reg_i = v_med(which_region==iregion,:); 
    v_reg(iregion, :) = median(v_reg_i, 1); % Median value at each depth. 
end

%% Plot mean velocities in the same plot
regions_titles = ["Appalachian", "Craton", "Grenville"]'; 
disp(regions_titles)

figure(2); clf; hold on; 
set(gcf,'pos', [1357 683 262 394]); 
title('Mean Vs profiles', 'fontweight', 'normal'); 
box on; grid on; set(gca,'LineWidth', 1.5); 
set(gca, 'YDir', 'reverse'); 
xlim([4, 4.8]); 
ylim([0, 250]); 
xlabel('Velocity (km/s)'); 
ylabel('Depth (km)'); 
        
for iregion = 1:length(regions_str); 
    v_i = v_reg(iregion, :); 
    plot(v_i, z_med, 'linewidth', 3, 'DisplayName', regions_titles(iregion));
end

legend(); 

% Now crust
figure(3); clf; hold on; 
box on; grid on; set(gca,'LineWidth', 1.5); 
set(gca, 'YDir', 'reverse'); 
xlim([1, 4.8]); 
ylim([0, 50]); 
xlabel('Velocity (km/s)'); 
ylabel('Depth (km)'); 
        
for iregion = 1:length(regions_str); 
    v_i = v_reg(iregion, :); 
    plot(v_i, z_med, 'linewidth', 3, 'DisplayName', regions_titles(iregion));
end

legend(); 

%% Plot all velocities
disp(regions_str)


plot_order = [3, 1, 2]; % Indecies corresponding to regions_titels

figure(1); clf; hold on; 
set(gcf, 'pos', [2570 777 560 420]); 
for irgn = 1:length(regions_str); % Plot in reverse order to match west to east 
    subplot(1,3,plot_order(irgn)); hold on; 
    box on; grid on; set(gca,'LineWidth', 1.5); 
    set(gca, 'YDir', 'reverse'); 
    xlim([4, 4.8]); 
    xlabel('Velocity (km/s)'); 
    ylabel('Depth (km)'); 
    title(regions_titles(irgn), 'FontWeight','normal'); 
    for istn = find(which_region==irgn); 
        v_i = v_med(istn,:); 
        plot(v_i, z_med, 'k', 'linewidth', 0.5); 
    end

    % Plot mean of the area
    plot(v_reg(irgn,:)', z_med', 'b', 'linewidth', 4); 
end

%% Calculate histograms at each depth in each region. 
plot_order = [3, 1, 2]; % Indecies corresponding to regions_titels

mcont_reg = zeros(length(z_med), length(vcont), length(regions_str)); 
for irgn = 1:length(regions_str); % Plot in reverse order to match west to east 
    for istn = find(which_region==irgn)'; % have to transpose for loop to work
        mcont_reg(:,:,irgn) = mcont_reg(:,:,irgn) + mcont(:,:,istn); 
    end
    mcont_reg(:,:,irgn) = mcont_reg(:,:,irgn) ./ sum(mcont_reg(:,:,irgn), 2) ./ mean(diff(vcont)); % Normalize to number of observations and to bin spacing. Make it an integral. 
end
%%
figure(4); clf; hold on; 
set(gcf, 'pos', [1070 777 560 420]); 
for irgn = 1:length(regions_str); % Plot in reverse order to match west to east 
    subplot(1,3,plot_order(irgn)); hold on; 
    box on; grid on; set(gca,'LineWidth', 1.5); 
    set(gca, 'YDir', 'reverse'); 
    xlim([4, 4.8]); 
    xlabel('Velocity (km/s)'); 
    ylabel('Depth (km)'); 
    title(regions_titles(irgn), 'FontWeight','normal'); 
    
    mcont_temp = mcont_reg(:,:,irgn); 
    mcont_temp(mcont_temp < 0.5)=nan; 
    contourf(vcont, z_med, mcont_temp, 20); 

end

%% Plot median and percentiles. 
figure(5); clf; hold on; 
set(gcf, 'pos', [500 815 233 382]); 

% Colors for each region. 
clr_reg = [[14, 207, 18]; ... % Appalachian
    [109, 64, 149]; ... % Craton
    [191, 94, 4]]./255; % Grennvile

make_lgd = []; % Hold legend handles. 

for irgn = 1:length(regions_str); % Plot in reverse order to match west to east 
    fprintf('Plotting %s\n',regions_titles(irgn))
    box on; grid on; set(gca,'LineWidth', 1.5); 
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca, 'YDir', 'reverse'); 
    xlim([4.3, 4.8]); 
    ylim([0, 250]); 
    xlabel('Velocity (km/s)'); 
    ylabel('Depth (km)'); 
    title('Regional Averages', 'FontWeight','normal'); 
    
%     medians = zeros(length(z_med), 1); 
    prct_vals = [0.1, 0.5, 0.9]; 
    percentiles = zeros(length(z_med), length(prct_vals)); 
    med_csum = nan(length(z_med)); 
    
    mcont_temp = mcont_reg(:,:,irgn); 
    csum = cumsum(mcont_temp, 2) .* mean(diff(vcont)); 
    csum = csum ./ max(csum,[],2); 

    uniqueify = linspace(0, 0.000001, length(vcont)); % Prevent non-unique error in interp
    for iprct = 1:length(prct_vals); 
        for iz = 1:length(z_med); 
            percentiles(iz, iprct) = interp1(csum(iz,:)+uniqueify, vcont+uniqueify, prct_vals(iprct)); 
        end
    end

    make_lgd(end+1) = plot(percentiles(:,prct_vals==0.5), z_med, ...
        'linewidth', 4, 'Color', clr_reg(irgn,:), 'DisplayName', regions_titles(irgn)); 
    plot(percentiles(:,1), z_med, ...
        'linewidth', 1, 'Color', clr_reg(irgn,:)); 
    plot(percentiles(:,end), z_med, ...
        'linewidth', 1, 'Color', clr_reg(irgn,:)); 

end
legend(make_lgd); 



%% Plot median and percentiles. FAST COPY PASTE CHANGE LIMITS FOR CRUST
figure(6); clf; hold on; 
set(gcf, 'pos', [500 815 233 382]); 

% Colors for each region. 
clr_reg = [[14, 207, 18]; ... % Appalachian
    [109, 64, 149]; ... % Craton
    [191, 94, 4]]./255; % Grennvile

make_lgd = []; % Hold legend handles. 

for irgn = 1:length(regions_str); % Plot in reverse order to match west to east 
    fprintf('Plotting %s\n',regions_titles(irgn))
    box on; grid on; set(gca,'LineWidth', 1.5); 
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(gca, 'YDir', 'reverse'); 
    xlim([3.0, 4.5]); 
    ylim([0, 50]); 
    xlabel('Velocity (km/s)'); 
    ylabel('Depth (km)'); 
    title('Regional Averages', 'FontWeight','normal'); 
    
%     medians = zeros(length(z_med), 1); 
    prct_vals = [0.1, 0.5, 0.9]; 
    percentiles = zeros(length(z_med), length(prct_vals)); 
    med_csum = nan(length(z_med)); 
    
    mcont_temp = mcont_reg(:,:,irgn); 
    csum = cumsum(mcont_temp, 2) .* mean(diff(vcont)); 
    csum = csum ./ max(csum,[],2); 

    uniqueify = linspace(0, 0.000001, length(vcont)); % Prevent non-unique error in interp
    for iprct = 1:length(prct_vals); 
        for iz = 1:length(z_med); 
            percentiles(iz, iprct) = interp1(csum(iz,:)+uniqueify, vcont+uniqueify, prct_vals(iprct)); 
        end
    end

    make_lgd(end+1) = plot(percentiles(:,prct_vals==0.5), z_med, ...
        'linewidth', 4, 'Color', clr_reg(irgn,:), 'DisplayName', regions_titles(irgn)); 
    plot(percentiles(:,1), z_med, ...
        'linewidth', 1, 'Color', clr_reg(irgn,:)); 
    plot(percentiles(:,end), z_med, ...
        'linewidth', 1, 'Color', clr_reg(irgn,:)); 

end
legend(make_lgd); 
%% Save all figs
exportgraphics(figure(1), 'vs_average_per_area.pdf'); 
exportgraphics(figure(2), 'vs_average_merged.pdf'); 
exportgraphics(figure(4), 'vs_heatmaps.pdf'); 
exportgraphics(figure(5), 'vs_percentiles.pdf'); 
exportgraphics(figure(6), 'vs_percentiles_crust.pdf'); 

