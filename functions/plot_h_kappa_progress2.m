function plot_h_kappa_progress2(trudata, allmodelsOrig, resdir, chainNo, ...
    accept_info, par, HK_Esum)
%% brb2022.02.09
% Make several plots of h-kappa inversion progress. 
% This applies for just one chain. 


%% Extract some info
hStack              = trudata.HKstack_P.H; 
kStack              = trudata.HKstack_P.K; 
eStack              = trudata.HKstack_P.Esum; 
E_by_Esuper_max     = trudata.HKstack_P.E_by_Esuper_max; 

ifaccept            = [accept_info.ifaccept            ]'; 
ifaccept = logical(ifaccept); 
misfit              = [accept_info.misfit              ]'; 
temper              = [accept_info.temp                ]'; 
Pm_prior            = [accept_info.Pm_prior1           ]'; 
p_bd                = [accept_info.p_bd                ]'; % Prob associated with birth or death
ptbnorm             = [accept_info.ptbnorm             ]'; 
Emax_per_iter       = [accept_info.hk_Emax_per_iter    ]'; % Maximum energy in hk stack made at any iteration


E2 = [misfit.E2]; 
E2 = [E2.HKstack_P]'; 


misfSurf = hKappaError(par.datprocess.HKappa.scale_error,...
                       trudata.HKstack_P.Esum,...
                       par.datprocess.HKappa.min_error,...
                       'HK_Esum', HK_Esum); % Misfit value at all spots on E stack. 

chi2 = [misfit.chi2]; 
chi2 = [chi2.HKstack_P]'; 

rms = [misfit.rms]; 
rms = [rms.HKstack_P]'; 

chi2sum = [misfit.chi2sum]; 
chi2sum = [chi2sum]'; 

log_lik             = [accept_info.log_likelihood      ]'; 
lik                 = 10 .^ log_lik; 
iIter               = [accept_info.iter                ]'; 
sig                 = [accept_info.sig_hk]'; 

allmodels           = [accept_info(:).model            ]; 
hIter               = [allmodels.zmoh                  ]'; 
kIter               = [allmodels.vpvs                  ]'; 


%% Delayed rejection stuff
nak = [accept_info.non_acceptk]'; 
yShift = 0.05; % For accepted or rejected
yPos = nak - yShift; 
yPos(ifaccept) = yPos(ifaccept) + 2 * yShift; 

figure(1425); clf; hold on; set(gcf, 'pos', [2156 1099 851 184]); 
grid on; 
box on; 
title('Model acceptance and delayed rejection'); 
xlabel('Iteration'); 
ylabel('Number of perturbations'); 
set(gca, 'YTick', [1, 2]); 
set(gca, 'TickLength', [0,0]); 
ylim([1-3*yShift, 2+3*yShift]); 


% Two perturbations. Use only this for legend. 
thisCond = and( ifaccept, nak == 2); 
ac2_perc = sum(thisCond) / length(thisCond) * 100; 
ac2 = scatter(iIter(thisCond), yPos(thisCond), '|b', 'DisplayName', ...
    sprintf('Accepted 2P: %3.0f%%',ac2_perc)); 
thisCond = and(~ifaccept, nak == 2); 
rj2_perc = sum(thisCond) / length(thisCond)* 100; 
rj2 = scatter(iIter(thisCond), yPos(thisCond), '|r', 'DisplayName', ...
    sprintf('Rejected 2P: %3.0f%%',rj2_perc)); 


% One perturbation
thisCond = and( ifaccept, nak == 1); 
ac1_perc = sum(thisCond) / length(thisCond)* 100; 
ac1 = scatter(iIter(thisCond), yPos(thisCond), '|b', 'DisplayName', ...
    sprintf('Accepted 1P: %3.0f%%',ac1_perc)); 
thisCond = and(~ifaccept, nak == 1); 
rj1_perc = sum(thisCond) / length(thisCond)* 100; 
rj1 = scatter(iIter(thisCond), yPos(thisCond), '|r', 'DisplayName', ...
    sprintf('Rejected 1P: %3.0f%%',rj1_perc)); 

leg = legend([ac1, rj1, ac2, rj2], 'Location', 'east','NumColumns',2);
legTitle = title(leg, 'Percent iterations per category'); 

% % % scatter([accept_info.iter], [accept_info.non_acceptk], 40, [accept_info.ifaccept], 'filled'); 
% % % scatter([accept_info.iter], [accept_info.non_acceptk], 40, [accept_info.ifaccept], 'filled'); 
% % % cbar=colorbar(); 
% % % colormap('cool'); 
% % % cbar.Label.String = 'Accept?'; 
% % % set(gcf, 'color', 'white'); 
% % % grid on; 

perc2pert = sum([accept_info.non_acceptk]==2) ./ length([accept_info.non_acceptk]) .* 100; 
% title(sprintf('Perturbed twice %2.1f percent of time', perc2pert)); 
exportgraphics(gcf, [resdir, '/delayed_perturbation_' num2str(chainNo) '.pdf']);  

%% Plots with h-k as x and y in plots. 
figure(1); clf; hold on; 
set(gcf, 'pos', [1550, 410, 1500, 1300]); 

% First ax. 
ax1hk = subplot(2,2,1); cla; hold on; 
title('hk-stack - accepted and rejected models'); 
ax1   = start_hk_ax(ax1hk, kStack, hStack, eStack); hold on; 
colorbar(ax1hk, 'location', 'north'); 
plot(ax1, kIter(ifaccept), hIter(ifaccept), 'LineWidth', 1', 'Color', 'k'); 
scatter(ax1, kIter(~ifaccept), hIter(~ifaccept), 10 , 'r' ); 
scatter(ax1, kIter(ifaccept ), hIter(ifaccept ), 50 , 'g' ); 
scatter(ax1, kIter(1        ), hIter(1        ), 110, 'k'); 
scatter(ax1, kIter(end      ), hIter(end      ), 110, '*k'); 
xlabel('Vp/Vs'); 
ylabel('Hmoh'); 

% Second ax. 
ax2hk = subplot(2,2,2); cla; hold on; 
ax2   = start_hk_ax(ax2hk, kStack, hStack, eStack); 
plot(kIter(~ifaccept), hIter(~ifaccept), 'LineWidth', 1/4', 'Color', 'r'); 
plot(kIter(ifaccept), hIter(ifaccept), 'LineWidth', 1', 'Color', 'k'); 
scatter(ax2, kIter, hIter, 10, iIter); 
xlabel('Vp/Vs'); 
ylabel('Hmoh'); 

% Third ax. 
ax3hk = subplot(2,2,3); cla; hold on; 
ax3   = start_hk_ax(ax3hk, kStack, hStack, eStack); 
plot(kIter(~ifaccept), hIter(~ifaccept), 'LineWidth', 1/4', 'Color', 'r'); 
plot(kIter(ifaccept), hIter(ifaccept), 'LineWidth', 1', 'Color', 'k'); 
scatter(ax3, kIter, hIter, 40, chi2);
title('chi2'); 
xlabel('Vp/Vs'); 
ylabel('Hmoh'); 

% Four ax. 
ax4hk = subplot(2,2,4); cla; hold on;
colormap(ax4hk, flipud(colormap)); 
title('Error(h-k stack)'); 
xlabel('Vp/Vs'); 
ylabel('Hmoh'); 
ax4   = start_hk_ax(ax4hk, kStack, hStack, misfSurf); hold on; 
colorbar(ax4hk, 'location', 'north'); 
linkaxes([ax4, ax4hk]); 

exportgraphics(gcf, [resdir '/convergence_info_chain_' num2str(chainNo) '_v0.png'], 'Resolution', 400);

%% Plots of liklihood and other useful things, across iterations. 
figure(2); clf; hold on; set(gcf, 'pos', [1345 1046 1985 785]);
nRow = 6; nCol = 2; 

subplot(nRow, nCol,1); hold on; grid on; box on;
plot(iIter, E2, 'k');
scatter(iIter(~ifaccept), E2(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), E2( ifaccept), 15, 'g');
ylabel('E');
set(gca, 'yscale', 'log');

subplot(nRow, nCol,2); hold on; grid on; box on;
plot(iIter, chi2, 'k');
scatter(iIter(~ifaccept), chi2(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), chi2( ifaccept), 15, 'g');
ylabel('chi2');
set(gca, 'yscale', 'log');

subplot(nRow, nCol,3); hold on; grid on; box on;
plot(iIter, rms, 'k');
scatter(iIter(~ifaccept), rms(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), rms( ifaccept), 15, 'g');
ylabel('rms');
set(gca, 'yscale', 'log');

subplot(nRow, nCol,4); hold on; grid on; box on;
plot(iIter, chi2sum, 'k');
scatter(iIter(~ifaccept), chi2sum(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), chi2sum( ifaccept), 15, 'g');
ylabel('chi2su');
set(gca, 'yscale', 'log');

subplot(nRow, nCol,5); hold on; grid on; box on;
plot(iIter, sig, 'k');
scatter(iIter(~ifaccept), sig(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), sig( ifaccept), 15, 'g');
ylabel('sigma');
set(gca, 'yscale', 'log');

subplot(nRow, nCol,7); hold on; grid on; box on;
plot(iIter, log_lik, 'k');
scatter(iIter(~ifaccept), log_lik(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), log_lik( ifaccept), 15, 'g');
ylabel('log\_lik');
set(gca, 'yscale', 'log');

subplot(nRow, nCol,6); cla; hold on; grid on; box on;
plot(iIter, sig, 'k');
scatter(iIter(~ifaccept), sig(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), sig( ifaccept), 15, 'g');
ylabel('sigma');
currentXlim = xlim(); 
currentXlim(1) = floor(1/30 * currentXlim(2)); 
xlim(currentXlim); 
set(gca, 'yscale', 'log');

subplot(nRow, nCol,8); hold on; grid on; box on;
plot(iIter, log_lik, 'k');
scatter(iIter(~ifaccept), log_lik(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), log_lik( ifaccept), 15, 'g');
ylabel('log\_lik');
xlim(currentXlim); 
set(gca, 'yscale', 'log');

% brb2022.02.09 Don't remember wha tthe commented code is below. 
% % % subplot(nRow, nCol, 9); cla; hold on; grid on; box on; 
% % % dLogLik = log_lik(2:end) - log_lik(1:end-1); 
% % % ifacceptDLog = ifaccept(2:end); 
% % % iIterDLog = iIter(2:end); % at each iteration, show change from previous iteration
% % % % plot(iIterDLog, dLogLik, 'k', 'LineWidth', 0.5); 
% % % scatter(iIterDLog(~ifacceptDLog), dLogLik(~ifacceptDLog), 15, 'r', 'filled');
% % % scatter(iIterDLog( ifacceptDLog), dLogLik( ifacceptDLog), 15, 'g');
% % % ylabel('log_lik change');
% % % % set(gca, 'yscale', 'log'); 
% % % 
% % % 
% % % 
% % % clipVal = 5; 
% % % subplot(nRow, nCol, 10); cla; hold on; grid on; box on; 
% % % dLogLik = log_lik(2:end) - log_lik(1:end-1); % Recalculate for bug prevention
% % % ifacceptDLog = ifaccept(2:end); 
% % % dLogLik(dLogLik >  clipVal) =  clipVal; 
% % % dLogLik(dLogLik < -clipVal) = -clipVal; 
% % % iIterDLog = iIter(2:end); % at each iteration, show change from previous iteration
% % % % plot(iIterDLog, dLogLik, 'k', 'LineWidth', 0.5); 
% % % scatter(iIterDLog(~ifacceptDLog), dLogLik(~ifacceptDLog), 15, 'r', 'filled');
% % % scatter(iIterDLog( ifacceptDLog), dLogLik( ifacceptDLog), 15, 'g');
% % % ylabel('log_lik change');
% % % xlim(currentXlim); 
% % % 
% % % clipVal = 0.5; 
% % % subplot(nRow, nCol, 12); cla; hold on; grid on; box on; 
% % % dLogLik = log_lik(2:end) - log_lik(1:end-1); % Recalculate for bug prevention
% % % ifacceptDLog = ifaccept(2:end); 
% % % dLogLik(dLogLik >  clipVal) =  clipVal; 
% % % dLogLik(dLogLik < -clipVal) = -clipVal; 
% % % iIterDLog = iIter(2:end); % at each iteration, show change from previous iteration
% % % % plot(iIterDLog, dLogLik, 'k', 'LineWidth', 0.5); 
% % % scatter(iIterDLog(~ifacceptDLog), dLogLik(~ifacceptDLog), 15, 'r', 'filled');
% % % scatter(iIterDLog( ifacceptDLog), dLogLik( ifacceptDLog), 15, 'g');
% % % ylabel('log_lik change');
% % % xlim(currentXlim); 

subplot(nRow, nCol, 9); cla; hold on; grid on; box on; 
plot(iIter, Pm_prior, 'k', 'LineWidth', 1,'DisplayName', 'p(m) prior'); 
scatter(iIter, p_bd, 'blue', 'DisplayName', 'p b/d'); 
scatter(iIter, ptbnorm, 'magenta', 'DisplayName', '|dm/m|?'); 
legend('Location', 'north'); 

subplot(nRow, nCol, 10); cla; hold on; grid on; box on; 
plot(iIter, temper); 
ylabel('Temperature'); 

subplot(nRow, nCol, 11); cla; hold on; grid on; box on; 
plot(iIter, Emax_per_iter);
xlabel('Iter'); 
title(sprintf(...
    'Max hk energy at every iteration. Highest ever obtainable: %s',...
    E_by_Esuper_max)); 

set(gcf, 'color', 'white');
exportgraphics(gcf, [resdir '/convergence_info_chain_' num2str(chainNo) '_v1.png'], 'Resolution', 400);



%% Plots of liklihood change through iterations. 
figure(3); clf; hold on; set(gcf, 'pos', [45 1046 1985 785]);

subplot(nRow, nCol,1); hold on; grid on; box on;
plot(iIter, log_lik, 'k');
scatter(iIter(~ifaccept), log_lik(~ifaccept), 15, 'r', 'filled');
scatter(iIter( ifaccept), log_lik( ifaccept), 15, 'g');
ylabel('log\_lik');
xlim(currentXlim); 
set(gca, 'yscale', 'log');

currentXlim = xlim();

subplot(nRow, nCol, 3); cla; hold on; grid on; box on; 
dLogLik = log_lik(2:end) - log_lik(1:end-1); 
ifacceptDLog = ifaccept(2:end); 
iIterDLog = iIter(2:end); % at each iteration, show change from previous iteration
scatter(iIterDLog(~ifacceptDLog), dLogLik(~ifacceptDLog), 15, 'r', 'filled');
scatter(iIterDLog( ifacceptDLog), dLogLik( ifacceptDLog), 15, 'g');
ylabel('log_lik change');
xlim(currentXlim); 

clipVal = 5; 
subplot(nRow, nCol, 5); cla; hold on; grid on; box on; 
dLogLik = log_lik(2:end) - log_lik(1:end-1); % Recalculate for bug prevention
ifacceptDLog = ifaccept(2:end); 
dLogLik(dLogLik >  clipVal) =  clipVal; 
dLogLik(dLogLik < -clipVal) = -clipVal; 
iIterDLog = iIter(2:end); % at each iteration, show change from previous iteration
% plot(iIterDLog, dLogLik, 'k', 'LineWidth', 0.5); 
scatter(iIterDLog(~ifacceptDLog), dLogLik(~ifacceptDLog), 15, 'r', 'filled');
scatter(iIterDLog( ifacceptDLog), dLogLik( ifacceptDLog), 15, 'g');
ylabel('log_lik change');
xlim(currentXlim); 

clipVal = 0.5; 
subplot(nRow, nCol, 7); cla; hold on; grid on; box on; 
dLogLik = log_lik(2:end) - log_lik(1:end-1); % Recalculate for bug prevention
ifacceptDLog = ifaccept(2:end); 
dLogLik(dLogLik >  clipVal) =  clipVal; 
dLogLik(dLogLik < -clipVal) = -clipVal; 
iIterDLog = iIter(2:end); % at each iteration, show change from previous iteration
% plot(iIterDLog, dLogLik, 'k', 'LineWidth', 0.5); 
scatter(iIterDLog(~ifacceptDLog), dLogLik(~ifacceptDLog), 15, 'r', 'filled');
scatter(iIterDLog( ifacceptDLog), dLogLik( ifacceptDLog), 15, 'g');
ylabel('log_lik change');
xlim(currentXlim); 

% Positive change. Log scale. 
subplot(nRow, nCol, 2); cla; hold on; grid on; box on; 
dLogLik = log_lik(2:end) - log_lik(1:end-1); 
dLogLik(dLogLik <=  0.01) = nan; % Only show positive changes
ifacceptDLog = ifaccept(2:end); 
iIterDLog = iIter(2:end); % at each iteration, show change from previous iteration
scatter(iIterDLog(~ifacceptDLog), dLogLik(~ifacceptDLog), 15, 'r', 'filled');
scatter(iIterDLog( ifacceptDLog), dLogLik( ifacceptDLog), 15, 'g');
ylabel('log_lik change\newline > 0.01');
set(gca, 'yscale', 'log'); 
xlim(currentXlim); 

% Negative change. Log scale. 
subplot(nRow, nCol, 4); cla; hold on; grid on; box on; 
dLogLik = log_lik(2:end) - log_lik(1:end-1); 
dLogLik(dLogLik >= -0.01) = nan; % Only show positive changes
dLogLik = - dLogLik; 
ifacceptDLog = ifaccept(2:end); 
iIterDLog = iIter(2:end); % at each iteration, show change from previous iteration
scatter(iIterDLog(~ifacceptDLog), dLogLik(~ifacceptDLog), 15, 'r', 'filled');
scatter(iIterDLog( ifacceptDLog), dLogLik( ifacceptDLog), 15, 'g');
ylabel('log_lik change\newline < 0.01\newline negative');
set(gca, 'yscale', 'log'); 
set(gca, 'ydir', 'reverse'); 
xlim(currentXlim); 

% What percent is accepted within a window
window = 300; 
indsEval = [window:length(ifaccept)]'; 
perAcc = zeros(size(indsEval)); 
% dLogLikFull = dLogLik
for ii = indsEval'; 
    whichAcc = ifaccept(ii+1-window:ii); 
    whichAcc = int64(whichAcc); 
    perAcc(ii+1-window) = sum( whichAcc ); 
end
perAcc = perAcc / window; 

subplot(nRow, nCol, 8); cla; hold on; grid on; box on;
plot(indsEval - .5 * window, perAcc); 
ylabel('Percent accepted'); 

xlim(currentXlim); 

exportgraphics(gcf, [resdir '/convergence_info_chain_' num2str(chainNo) '_v2.png'], 'Resolution', 400);

end