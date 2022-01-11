function plot_h_kappa_progress2(trudata, allmodelsOrig, resdir, chainNo, accept_info)

% Extract some info
hStack              = trudata.HKstack_P.H; 
kStack              = trudata.HKstack_P.K; 
eStack              = trudata.HKstack_P.Esum; 

ifaccept            = [accept_info.ifaccept            ]'; 
misfit              = [accept_info.misfit              ]'; 
temper              = [accept_info.temp                ]'; 

E2 = [misfit.E2]; 
E2 = [E2.HKstack_P]'; 

chi2 = [misfit.chi2]; 
chi2 = [chi2.HKstack_P]'; 

rms = [misfit.rms]; 
rms = [rms.HKstack_P]'; 

chi2sum = [misfit.chi2sum]; 
chi2sum = [chi2sum]'; 



% E2                  = [misfit.E2.HKstack_P             ]'; 
log_lik             = [accept_info.log_likelihood      ]'; 
lik                 = 10 .^ log_lik; 
iIter               = [accept_info.iter                ]'; 
sig                 = [accept_info.sig_hk]'; 

allmodels           = [accept_info(:).model            ]; 
hIter               = [allmodels.zmoh                  ]'; 
kIter               = [allmodels.vpvs                  ]'; 



% Start basic figure properties
figure(1); clf; hold on; 
set(gcf, 'pos', [1550, 410, 1500, 1300]); 

% First ax. 
ax1hk = subplot(2,2,1); cla; hold on; 
ax1   = start_hk_ax(ax1hk, kStack, hStack, eStack); 

plot(ax1, kIter(ifaccept), hIter(ifaccept), 'LineWidth', 1', 'Color', 'k'); 
scatter(ax1, kIter(ifaccept), hIter(ifaccept), 50, 'g'); 
scatter(ax1, kIter(~ifaccept), hIter(~ifaccept), 10, 'r'); 

% Second ax. 
ax2hk = subplot(2,2,2); cla; hold on; 
ax2   = start_hk_ax(ax2hk, kStack, hStack, eStack); 
plot(kIter(~ifaccept), hIter(~ifaccept), 'LineWidth', 1/4', 'Color', 'r'); 
plot(kIter(ifaccept), hIter(ifaccept), 'LineWidth', 1', 'Color', 'k'); 
scatter(ax2, kIter, hIter, 10, iIter); 


% Third ax. 
ax3hk = subplot(2,2,3); cla; hold on; 
ax3   = start_hk_ax(ax3hk, kStack, hStack, eStack); 

plot(kIter(~ifaccept), hIter(~ifaccept), 'LineWidth', 1/4', 'Color', 'r'); 
plot(kIter(ifaccept), hIter(ifaccept), 'LineWidth', 1', 'Color', 'k'); 
scatter(ax3, kIter, hIter, 40, chi2);
title('chi2'); 

exportgraphics(gcf, [resdir '/convergence_info_scatter_chain_' num2str(chainNo) '_hk_surface.png'], 'Resolution', 400)




%% Other plot. 
figure(2); clf; hold on; set(gcf, 'pos', [1345 1046 1985 785]);
nRow = 5; nCol = 2; 

subplot(nRow, nCol,1); hold on; grid on; box on;
plot(iIter, E2, 'k');
scatter(iIter(~ifaccept), E2(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), E2( ifaccept), 15, 'g', 'filled');
ylabel('E');
set(gca, 'yscale', 'log')


subplot(nRow, nCol,2); hold on; grid on; box on;
plot(iIter, chi2, 'k');
scatter(iIter(~ifaccept), chi2(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), chi2( ifaccept), 15, 'g', 'filled');
ylabel('chi2');
set(gca, 'yscale', 'log')


subplot(nRow, nCol,3); hold on; grid on; box on;
plot(iIter, rms, 'k');
scatter(iIter(~ifaccept), rms(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), rms( ifaccept), 15, 'g', 'filled');
ylabel('rms');
set(gca, 'yscale', 'log')


subplot(nRow, nCol,4); hold on; grid on; box on;
plot(iIter, chi2sum, 'k');
scatter(iIter(~ifaccept), chi2sum(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), chi2sum( ifaccept), 15, 'g', 'filled');
ylabel('chi2su');
set(gca, 'yscale', 'log')




subplot(nRow, nCol,5); hold on; grid on; box on;
plot(iIter, sig, 'k');
scatter(iIter(~ifaccept), sig(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), sig( ifaccept), 15, 'g', 'filled');
ylabel('sigma');
set(gca, 'yscale', 'log')


subplot(nRow, nCol,7); hold on; grid on; box on;
plot(iIter, log_lik, 'k');
scatter(iIter(~ifaccept), log_lik(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), log_lik( ifaccept), 15, 'g', 'filled');
ylabel('log\_lik');
set(gca, 'yscale', 'log')




subplot(nRow, nCol,6); cla; hold on; grid on; box on;
plot(iIter, sig, 'k');
scatter(iIter(~ifaccept), sig(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), sig( ifaccept), 15, 'g', 'filled');
ylabel('sigma');
currentXlim = xlim(); 
currentXlim(1) = floor(1/30 * currentXlim(2)); 
xlim(currentXlim); 
set(gca, 'yscale', 'log')

subplot(nRow, nCol,8); hold on; grid on; box on;
plot(iIter, log_lik, 'k');
scatter(iIter(~ifaccept), log_lik(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), log_lik( ifaccept), 15, 'g', 'filled');
ylabel('log\_lik');
xlim(currentXlim); 
set(gca, 'yscale', 'log')



subplot(nRow, nCol, 9); cla; hold on; grid on; box on; 
plot(iIter, temper); 
ylabel('Temperature'); 

% subplot(nRow, nCol, 10); cla; hold on; grid on; box on; 
% plot(iIter, temper); 
% ylabel('Temperature'); 
% set(gca, 'yscale', 'log')



exportgraphics(gcf, [resdir '/convergence_info_scatter_chain_' num2str(chainNo) '.png'], 'Resolution', 400)




end