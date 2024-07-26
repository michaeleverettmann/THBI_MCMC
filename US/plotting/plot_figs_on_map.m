%% Startup. Path definitions
clc; clear; 
run('../../a0_STARTUP_BAYES.m');
addpath('~/MATLAB/borders');
%%
resdir_computer = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_cnsi'; % /Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_eri/cont_EProt-s1m1m2_testnwk_dat0/low_noise_sp
resdir_fig = '/Users/brennanbrunsvik/Documents/temp'; 
prior_path  = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ENAM/prior.mat' ; 
fstainfo = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/US/INFO/stations.mat';

stainfo = load(fstainfo); 
stainfo = stainfo.stainfo; 
stalats = stainfo.slats; 
stalons = stainfo.slons; 

netall = string(stainfo.nwk); 
staall = string(stainfo.stas); 

stamp = 'more_anis'; 

figure(1); clf; hold on; 
xlimits = [-130, -100]; 
ylimits = [30, 45]; 
xlim(xlimits); 
ylim(ylimits); 

[latbord, lonbord] = borders('states'); % add states map
for iline = 1:length(latbord); 
    line(lonbord{iline}, latbord{iline}, 'color', 'k'); 
end
% for iplace = 1:length(lonbord); 
%     m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
% end

% Loop through each folder
folds = dir(resdir_computer); 
for i = 1:length(folds)
    if folds(i).isdir && ~ismember(folds(i).name, {'.', '..'})
        
        % Get station name etc if we have a result here. 
        name = folds(i).name;
        fold = sprintf('%s/%s/%s',resdir_computer, name, stamp); 
        if ~ exist(sprintf('%s/final_model.mat',fold)); 
            continue
        end
        splt = split(name, '_'); 
        sta = splt{1}; 
        net = splt{2}; 
        
        % 
        iallsta = find((sta == staall) & (net == netall)); 
        if length(iallsta) == 0; 
            continue
        end 

        lat = stalats(iallsta); 
        lon = stalons(iallsta); 

        % scatter(lon, lat); 
        fig_model = sprintf('%s/heatmap_of_models.pdf',fold); % final_true_vs_pred_data_wavs.png
        % fig_model = sprintf('%s/final_true_vs_pred_data_wavs.png',fold);
        % img_model = fig_model; 
        img_model = PDFtoImg(fig_model); 

        if exist(img_model, 'file')
            % Insert the image at specified location
            img = imread(img_model);
            [img_height, img_width, ~] = size(img);
            
            % Define the position and size for the image
            img_pos_x = (lon-xlimits(1))/diff(xlimits); % lon;  % X position (longitude)
            img_pos_y = (lat-ylimits(1))/diff(ylimits); % lat;  % Y position (latitude)
            img_width_scaled = 0.05;  % Scaled width
            img_height_scaled = (img_height / img_width) * img_width_scaled;  % Scaled height based on aspect ratio
            
            % Insert the image
            axes('Position', [img_pos_x, img_pos_y, img_width_scaled, img_height_scaled]);
            imshow(img);
        end

        delete(img_model); 
        
        % break 

    end
end

exportgraphics(gcf, './test.pdf', 'Resolution',2400, 'ContentType','image'); 
