function []=plot_data_heatmap(datamap, trudata, predata, par); 

rf = datamap.RF_Sp_ccp(:); 
rf_obs = trudata.RF_Sp_ccp.PSV(:,1); 
z_obs = trudata.RF_Sp_ccp.zz; 

virt_lim_expect = max(abs(rf_obs)); 

nz=length(rf(1).rf); 
z = predata.RF_Sp_ccp.zz; % Temporary; 
rfgrid = zeros(nz,length(rf));
for irf = 1:length(rf); 
    rfgrid(:,irf) = rf(irf).rf; 
end

% Receiver function at iteration
fprintf('\nPlotting receiver function per iteration.\n'); 
figure(1); clf; hold on; box on; set(gca, 'lineWidth', 1.5); 
set(gcf, 'color', 'white', 'pos', [15 64 582 1255]); 
xlabel('Depth (km)'); 
ylabel('Iteration'); 
load('vik'); 
colormap(vik);
contourf(z,[1:size(rfgrid,2)],rfgrid', length(colormap()), 'edgecolor', 'none');
maxCol = 3/2 *  virt_lim_expect; %.5 * max(max(abs(rfgrid))); 
caxis([-maxCol,maxCol]); 
colorbar('northoutside'); 
exportgraphics(gcf, sprintf('%s/rf_at_iter_contour_%s.png',par.res.resdir,par.res.chainstr)); 



% PDF of receiver functions
nbin = 100; % Bins to apply in hist. 
rfvals = linspace(-3/2*max(abs(rf_obs)),3/2*max(abs(rf_obs)),nbin)'; % Values to give to hist. 
% rfvals = linspace(prctile(rfgrid,0.2,"all"), prctile(rfgrid,99.8,"all"), nbin)'; 


hrf = zeros(nz, nbin); % histogram of receiver function.
for iz = 1:nz; 
    hist_z = hist(rfgrid(iz, :)',rfvals); 
    hrf(iz,:) = hist_z'; % COunt rfs in each bin at each depth. 
    hrf(iz,:) = hrf(iz,:) ./ sum(hrf(iz,:)); % Normalize 
end
hrf(hrf==0)=nan; 

%%% PDFs of receiver functions
fprintf('\nPlotting receiver function PDF.\n'); 
figure(2); clf; hold on; set(gcf, 'color', 'white', 'pos', [252 210 920 500]); 
ncont = 100; 

med_rf = nanmedian(rfgrid,2); 
mea_rf = nanmean  (rfgrid,2); 

% subplot(2,1,1); cla; hold on; box on; set(gca, 'lineWidth', 1.5); 
% contourf(z, rfvals', hrf', ncont, 'edgecolor', 'none'); % pclr = pcolor(z, rfvals', hrf'); pclr.LineStyle = 'none'; 
% cbar = colorbar(); 
% cbar.Label.String = 'p(d)';
% colormap(viridis);
% xlabel('Depth'); 
% plot([min(z), max(z)], [0,0], 'k', 'linewidth', 0.25); 
% rf_obs_hand = plot(z_obs, rf_obs, 'k'               , 'linewidth', 2.5); 
% med_rf_hand = plot(z    , med_rf, 'red'             , 'linewidth', 1.5); 
% mea_rf_hand = plot(z    , mea_rf, ...
%     'color', [28, 221, 235]./256, 'linewidth', 1.5); 

subplot(1,1,1); cla; hold on; box on; set(gca, 'lineWidth', 1.5); 
hrf_log = log(hrf); 
hrf(isinf(hrf_log)) = -20; 
contourf(z, rfvals', hrf_log', ncont, 'edgecolor', 'none'); % pclr = pcolor(z, rfvals', hrf'); pclr.LineStyle = 'none'; 
cbar = colorbar(); 
cbar.Label.String = 'ln(p(d))';
colormap(viridis);
xlabel('Depth'); 
plot([min(z), max(z)], [0,0], 'k', 'linewidth', 0.25); 
rf_obs_hand = plot(z_obs, rf_obs, 'k'               , 'linewidth', 2.5); 
med_rf_hand = plot(z    , med_rf, 'red'             , 'linewidth', 1.5); 
mea_rf_hand = plot(z    , mea_rf, ...
    'color', [28, 221, 235]./256, 'linewidth', 1.5); 


ylim([ -3/2 * virt_lim_expect, 3/2 * virt_lim_expect ]); 

leg=legend([rf_obs_hand, med_rf_hand, mea_rf_hand],...
    'Observed', 'Median', 'Mean'); 
leg.Location = 'best'; 

exportgraphics(gcf, sprintf('%s/rf_pdf_%s.png',par.res.resdir,par.res.chainstr)); 

save(sprintf('%s/rf_pdf_info.mat',par.res.resdir), 'datamap', 'predata', 'par'); 


end