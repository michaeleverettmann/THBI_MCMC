
%% This cell for degrees of freedom. 
figure(1); clf; hold on; 
subplot(1,1,1); 
set(gcf, 'Position', [1457 1232 560 420]); 
set(gca, 'CameraPosition', [-133.4370 10.0526 125.2517]); 
% n=linspace(2,5,20)'; 
n = [2: 20]; 
phi = linspace(.01, 2, 1000); 
[n, phi] = ndgrid(n, phi); 
sigma = sqrt(phi./n); 
pdm = 1 ./ (sqrt(2 .* pi).*sigma).^n .* exp(-phi./(2.*sigma.^2)); 
pdm = log(pdm); 

% plot(n, pdm)
surf(n, phi, pdm, 'EdgeAlpha', 0.1); 
xlabel('N'); 
ylabel('phi'); 
zlabel('log10(p(d|m)'); 
% title(['phi = ' num2str(phi)]); 
title('p(d|m)')

colorbar(); 
grid(); 
exportgraphics(gcf, 'figures_general/pdm_phi_N.jpeg', 'Resolution', 500)



%% This cell looks at the error function. 
% load('/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ENAM/hk_temp_data.mat')
load('hk_temp_data.mat'); % Should be able to make this just by saving an h_kappa structure at some point during the inversion. 
hk = trudata.HKstack_P; 
k = hk.K; 
h = hk.H; 
Esum = hk.Esum / max(max(hk.Esum)); 
figure(2); clf; hold on; set(gcf, 'pos', [0 630 1909 1427]); 
CameraPosition = [-1.9281 186.9332 17.9083]; 
subplot(2,2,1); 
thing = surf(k', h, Esum', 'EdgeAlpha', 0)
xlabel('K'); 
ylabel('H'); 
set(gca, 'ydir', 'reverse')
title('Energy (normalized)'); 
% set(gca, 'CameraPosition', CameraPosition); 

pen1 = 1./Esum; 

clip1 = 6
pen1(pen1>clip1) = clip1; 
pen2 = (1-Esum) + 0.05; 
% pen2 = pen2 .^2
pen3 = pen2 .^ (1/4)

subplot(2,2,2);
thing2 = surf(k', h, pen1', 'EdgeAlpha', 0); 
xlabel('K'); 
ylabel('H'); 
set(gca, 'ydir', 'reverse')
title('1/E (clipped)'); 
% set(gca, 'CameraPosition', CameraPosition); 


subplot(2,2,3); 
thing3 = surf(k', h, pen2', 'EdgeAlpha', 0); 
xlabel('K'); 
ylabel('H'); 
set(gca, 'ydir', 'reverse'); 
title('1-Esum + small constant'); 
% set(gca, 'CameraPosition', CameraPosition); 



subplot(2,2,4); 
thing3 = surf(k', h, pen3', 'EdgeAlpha', 0); 
xlabel('K'); 
ylabel('H'); 
set(gca, 'ydir', 'reverse'); 
title('(1-Esum + small constant)^(1/2)'); 
% set(gca, 'CameraPosition', CameraPosition); 

exportgraphics(gcf, '../../figures_general/h_kappa_penalty_surfaces.jpeg', ...
    'Resolution', 300); 