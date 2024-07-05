function [hmps,hmpp] = plot_HEATMAP_ALLMODELS(suite_of_models, final_model, par, ifsave, ofile)
% [hmps,hmpp] = plot_HEATMAP_ALLMODELS(suite_of_models,par,ifsave,ofile)
% 
% function to plot a heatmap pdf of the velocity models

if nargin<4 || isempty(ifsave)
    ifsave = false; % default is not to save
end
if nargin<5 || isempty(ofile)
    ofile = 'figs/heatmap_models_fig.pdf';
end

%% GET TARGET MODEL for comparison
sm = suite_of_models;

Xvs = [par.mod.sed.vsmin:0.04:par.mod.mantle.vsmax];
Xvp = [1.3*par.mod.sed.vsmin:0.04:par.mod.mantle.vsmax*2];
Xsa = 100*[(par.mod.crust.ximin-1):0.002:(par.mod.crust.ximax-1)]; 
Nz  = length(sm(1).Z);

hmps = zeros(Nz,length(Xvs));
hmpp = zeros(Nz,length(Xvp));
hmsa = zeros(Nz,length(Xsa)); 
nsaved = 0;

% loop over depths, adding 
for iz = 1:Nz
    hmps(iz,:) = hist(sm.VS(iz,:),Xvs);
    hmpp(iz,:) = hist(sm.VP(iz,:),Xvp);
    hmsa(iz,:) = hist(sm.Sanis(iz,:),Xsa); 
end
% norm pdf to 1
hmps = hmps/size(sm.VS,2);
hmpp = hmpp/size(sm.VS,2);
hmsa = hmsa/size(sm.Sanis,2); 
% log pdf
hmps = log(hmps); hmps(isinf(hmps)) = -20;
hmpp = log(hmpp); hmpp(isinf(hmpp)) = -20;
hmsa = log(hmsa); hmsa(isinf(hmsa)) = -20; 


%% Plot
fig_num = 25; 
figure(fig_num); clf; set(gcf,'pos',[-1011 247 841 721]); 
tiledlayout(1,3,"TileSpacing",'compact'); 
ax1 = nexttile(1); set(ax1, 'linewidth', 1.5); hold on; box on; 
ax2 = nexttile(2); set(ax2, 'linewidth', 1.5); hold on; box on; 
ax3 = nexttile(3); set(ax3, 'linewidth', 1.5); hold on; box on; 

%% data
contourf(ax1,Xvs,sm.Z,hmps,[-5:0.1:-0.1],'edgecolor','none');
contourf(ax2,Xvp,sm.Z,hmpp,[-5:0.1:-0.1],'edgecolor','none');
contourf(ax3,Xsa,sm.Z,hmsa,[-5:0.1:-0.1],'edgecolor','none');

%% mean model
fm_line_width = 2; 
fm_color      = 'red'; 
plot(ax1, final_model.VSav, final_model.Z, ...
    'LineWidth', fm_line_width, 'Color', fm_color); 
plot(ax2, final_model.VPav, final_model.Z, ...
    'LineWidth', fm_line_width, 'Color', fm_color); 
plot(ax3, final_model.Sanisav, final_model.Z, ...
    'LineWidth', fm_line_width, 'Color', fm_color); 

%% pretty
xlim(ax1,[3 4.9]);
xlim(ax2,[5.4 8.8]);
xlim(ax3,[-20, 20]);

set(ax1,'ydir','reverse','fontsize',15,'ytick',[0:25:max(sm.Z)],'color','none');
set(ax2,'ydir','reverse','fontsize',15,'yticklabel','','ytick',[0:25:max(sm.Z)],'color','none');
set(ax3,'ydir','reverse','fontsize',15,'yticklabel','','ytick',[0:25:max(sm.Z)],'color','none');

xlabel(ax1,'Vs (km/s)','fontsize',18)
ylabel(ax1,'Depth (km)','fontsize',18)
xlabel(ax2,'Vp (km/s)','fontsize',18)
xlabel(ax3, 'S anis (%)','fontsize',18)

% colours
caxis(ax1,[-5 0])
caxis(ax2,[-5 0])
caxis(ax3,[-5 0])
hcb = colorbar(ax3);
set(hcb,'pos',[0.151 0.13 0.02 0.215],'fontsize',14,'yaxislocation','right')
hcby = ylabel(hcb,'log_{10}(Probability)','fontweight','bold','verticalalignment','middle','HorizontalAlignment','center');
set(hcby,'rotation',270,'pos',[3.2 -2.5 0]);

%% Plot true model if synthetic. 
if par.inv.synthTest; 
    global TRUEmodel
    if ~isempty(TRUEmodel)
        Z = TRUEmodel.Z;
        vs = TRUEmodel.VS;
        vp = TRUEmodel.VP;
        sa = TRUEmodel.Sanis; 
        plot(ax1,vs,Z,'-b','Linewidth',2);
        plot(ax2,vp,Z,'-b','Linewidth',2);
        plot(ax3,vp,Z,'-b','Linewidth',2); 
    else
        warning('Why is TRUEmodel empty? Cant plot TRUEmodel'); 
    end
end

%% SAVE
if ifsave
    fprintf('   saving fig\n');
    save2pdf(fig_num,ofile,'/');
end

end