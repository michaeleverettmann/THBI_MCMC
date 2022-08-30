figure(2120); clf; hold on; 
poslims = @()set(gca, 'ydir', 'reverse'); 

zz = prior.zatdep ; 

Xvs = [par.mod.sed.vsmin:0.04:par.mod.mantle.vsmax];
Xvp = [1.3*par.mod.sed.vsmin:0.04:par.mod.mantle.vsmax*2];
Nz  = length(zz);

hmps = zeros(Nz,length(Xvs));
hmpp = zeros(Nz,length(Xvp));
% nsaved = 0;

% loop over depths, adding 
for iz = 1:Nz
    hmps(iz,:) = hist(prior.VSmantle(:,iz),Xvs);
    hmpp(iz,:) = hist(prior.VSmantle(:,iz),Xvp);
end
% norm pdf to 1
hmps = hmps/size(prior.VSmantle,2);
hmpp = hmpp/size(prior.VSmantle,2);
% log pdf
hmps = log(hmps); hmps(isinf(hmps)) = -20;
hmpp = log(hmpp); hmpp(isinf(hmpp)) = -20;


%% Plot
fig_num = 25; 
figure(fig_num); clf; set(gcf,'pos',[-1011 247 841 721])
ax1 = axes(gcf,'pos',[0.13 0.11 0.35 0.815], 'linewidth', 1.5); hold on; box on; 
poslims(); 
ax2 = axes(gcf,'pos',[0.52 0.11 0.35 0.815], 'linewidth', 1.5); hold on; box on; 
poslims(); 
ax3 = axes(gcf,'pos',[0.91 0.11 0.02 0.815], 'linewidth', 1.5); hold on; box on; 


%% data
contourf(ax1,Xvs,zz,hmps,[-5:0.1:-0.1],'edgecolor','none');
contourf(ax2,Xvp,zz,hmpp,[-5:0.1:-0.1],'edgecolor','none');

%%
figure(2120); clf; hold on; 
set(gcf, 'color', 'white'); 
tiledlayout(2,2, 'TileSpacing', 'compact'); 

nexttile(1); 
hist(prior.VSsedtop, 20); 
title('VS sed top'); 

nexttile(3); 
hist(prior.VSsedbot, 20); 
title('VS sed bot'); 

nexttile(2); 
hist(prior.VScrusttop, 20); 
title('VS crust top'); 

nexttile(4); 
hist(prior.VScrustbot, 20); 
title('VS crust bot'); 