%%% This suggests how to recover a pdf from MCMC. It's not perfect. 
% 1. Sort sampled velocity values v. 
% 2. Assume cumulative probability of getting the next higher v moves in
% equal steps. i.e. p_cumulative = interp1(v_sorted,linspace(0,1,nv) )
% 3. Get pdf from cumulative pdf using derivative. pdf = diff(pdf_cumulative) ./
% diff(v)
% 4. There will be some roughness due to imperfect sampling. We need to do
% smoothing. A big challenge here is the diff(v) in denominator... this
% makes things very instable if there was two very similar v values.
% I tried making v in two ways. 
% I tried using truly random v from normrnd. 
% This produces rough PDFs that need a lot of smoothing. 
% But, this is analogous to how MCMC really works. 
% I also tried making v with spacing determined by the pdf. Analytically,
% there was finer spacing where the pdf is higher. This returned almost
% perfect pdfs in the final plot. But, this isn't how MCMC really works...
% It does, however, demonstrate that this process is (at least partially)
% mathematically valid. 

figpath = '~/Documents/UCSB/ENAM/THBI_ENAM/figures/testing/hist_creation/'; 

nv = 2000; % Number of velocity in analytical pdf. 
n_samps = 2000; % Number of sampled velocities. An artificial MCMC result.
v = linspace(0.5, 5.5, nv)'; % V grid. 
sigma = 0.3; % Analytical error of v. 
mean_gaus = 3.4; % Analytical mean_gaus of v. 
smooth_window = .15; % How much to smooth final recovered pdf. In percent. 


pdf_an = exp(-(v-mean_gaus).^2./(2.*sigma.^2)); % Analytical pdf of velocity. 
pdf_an = pdf_an ./ trapz(v, pdf_an); % Normalize the pdf. 

pdf2 = true; 
pdf3 = true; 
if pdf2; 
    sigma2 = 0.1; 
    mean_gaus2 = 2; 
    pdf_an2 = exp(-(v-mean_gaus2).^2./(2.*sigma2.^2)); 
    pdf_an2 = pdf_an2 ./ trapz(v, pdf_an2); % Normalize the pdf. 
    pdf_an = pdf_an + pdf_an2; 
end 
if pdf3; 
    sigma3 = 0.5; 
    mean_gaus3 = 4.5; 
    pdf_an3 = exp(-(v-mean_gaus3).^2./(2.*sigma3.^2)); 
    pdf_an3 = pdf_an3 ./ trapz(v, pdf_an3); % Normalize the pdf. 
    pdf_an = pdf_an + pdf_an3; 
end 
pdf_an = pdf_an ./ trapz(v, pdf_an); % Normalize the pdf. 

pdf_cumulative = cumtrapz(v, pdf_an); % Cumulative pdf. 
% v_samps = interp1(pdf_cumulative, v, linspace(0,0.9999,n_samps)'); % NON random v that follows pdf. 
v_samps = normrnd(mean_gaus, sigma, [n_samps,1]); % Random v that follows pdf. Analogous to MCMC output. 
% v_samps = [min(v); v_samps; max(v)]; 
if pdf2; 
    v_samps = [v_samps; normrnd(mean_gaus2, sigma2, [n_samps,1])]; 
end
if pdf3; 
    v_samps = [v_samps; normrnd(mean_gaus3, sigma3, [n_samps,1])]; 
end
v_samps = sort(v_samps); % Sort to make comparable to a cumulative pdf. 
p_v_samps = interp1(v, pdf_an, v_samps); % P for each sampled v (only for the smoothly distributed v). 

fprintf('You ended up with %1.0f samples. ', length(v_samps))

figure(2); clf; hold on; 
box on; 
hist(v_samps, 100)
plot(v, pdf_an, 'DisplayName','PDF  (analytical)')
scatter(v_samps, p_v_samps, 15, 'filled', 'DisplayName', 'Sampled V'); 

%% One attempt to get pdf of v, using differences between v and all sampled v. 
vpen = zeros(size(v)); % A "penalty" approach. 
for iv = 1:length(v);
    vi = v(iv); 
    vdiff = abs(vi - v_samps); 
    vdiff = vdiff.^2; % Square
%     vdiff = abs(vdiff); 
    vpen(iv) = sum(vdiff); 
end

p_from_pen = 1./vpen; 
p_from_pen = p_from_pen /  (trapz(v, p_from_pen)); 
% p_from_pen = p_from_pen .^ (1.2); 
p_from_pen = p_from_pen /  (trapz(v, p_from_pen)); 


figure(3); clf; hold on; 
plot(v, 1./p_from_pen, 'DisplayName','1/penalty')
legend(); 


%% Cumulative probability style of recoverying pdf. 
% This almost works properly, and is very smooth, for a single Gaussian. It
% doesn't work for multiple gaussians sumed together though. 
prob_cum_vsamp = linspace(0, 1, length(v_samps))'; 
sort(v_samps); % SORT VSAMP!!! Is already done earlier, but left here for clarity. 
dv_p = diff(v_samps); % DV part of derivative. 
% newknt()
dv_p_n = length(dv_p); 
dv_p = [ones(size(dv_p)).*max(dv_p); dv_p; ones(size(dv_p)).*mean(dv_p(end-50:end))]; 
dv_p = smoothdata(dv_p, 'gaussian', n_samps .* smooth_window); 
dv_p = dv_p(dv_p_n+1:end-dv_p_n); 
dp_p = diff(prob_cum_vsamp); % Dp part of derivative. 
pv_recover = dp_p ./ dv_p; 
% pv_recover = smoothdata(pv_recover, 'gaussian', n_samps .* smooth_window); 
v_samp_recover = ( v_samps(1:end-1) + v_samps(2:end) ) ./ 2; % Because we loose one value doing derivative. 

%% Kernel density estimation, from Matlab. 
100; % nptsksd = int16((max(v_samps) - min(v_samps) ) ./ .01); % Number of points to use in estimating k density. I'm spacing it for now to have a point each 0.01 km/s
[pv_k, v_k] = ksdensity(v_samps, 'width', .04, 'NumPoints', nptsksd);%, 'Bandwidth',0.03); 
pv_k = pv_k / trapz(v_k, pv_k); 

%% Final plots. 
LW = 4; % Linewidths
figure(6); clf; hold on; 
set(gcf,'color', 'white'); 
histogram(v_samps, 100, 'Normalization', 'pdf', 'EdgeColor','auto'); 
plot(v_samp_recover, pv_recover, 'DisplayName','p(v) recovered','linewidth',LW); 
plot(v, pdf_an, 'DisplayName','p(v) analytical','linewidth',LW); 
plot(v, p_from_pen, 'DisplayName', 'p(v) from penalty approach.','linewidth',LW); 
plot(v_k, pv_k, 'DisplayName', 'p(v) kernel density estimation','linewidth',LW);      
this_title = sprintf('Compare analytical pdf, to pdf recovered from samples.\n Smooth window = %1.1f%%',...
    smooth_window*100); 
title(this_title, 'FontWeight','normal'); 
box on; 
legend(); 
xlabel('v'); 
ylabel('p(v)'); 
grid on; 

exportgraphics(gcf, [figpath 'hist_extimations.pdf']); 

pdf_example = struct('v', v_k, 'p', pv_k); 
% save('pdf_example.mat', 'pdf_example'); 