
%Script to grab a sta and make sensitivity kernels from the final model. 
clear
close all

ddir = dir('/Users/jon/Documents/Work/Codes/THBI/EAR/STAsinv/*_*_dat1');

for id = 1:length(ddir)%Loop through desired stations
    clearvars -except id ddir store
    dpath = [ddir(id).folder,'/',ddir(id).name,'/PHV_RF/'];%Get correct

    load([dpath,'par.mat']);
    load([dpath,'final_model.mat']);
    load([dpath,'final_predata.mat']);
    load([dpath,'allmodels_perchain.mat']);
    load([dpath,'posterior.mat']);

    %Get all models with anisotropy
    allmodels_collated = [];

    if par.inv.nchains==1
        allmodels_collated = allmodels_perchain;
        return
    end

    for iii = 1:length(allmodels_perchain)
        am = allmodels_perchain{iii};
        if isempty(am), continue; end
        am([am.iter]'<=par.inv.burnin) =  [];
        am = dealto(am,'chain',iii);

        allmodels_collated = [allmodels_collated,am];
    end

    %make the average model, but only take what's necessary.  
    Z = final_model.Z;
    Nz = length(Z);
    Nm = length(allmodels_collated);
    i50     = round(0.5*Nm);%50th percentile

    for ii = 1:Nm

        posterior.Sanis(:,ii) = linterp(allmodels_collated(ii).Z,allmodels_collated(ii).Sanis,Z);
        posterior.Panis(:,ii) = linterp(allmodels_collated(ii).Z,allmodels_collated(ii).Panis,Z);
    end

    % Get median S,P anis
    for iz = 1:Nz
        s_sort = sort(posterior.Sanis(iz,:));
        p_sort = sort(posterior.Panis(iz,:));
        Sanis(iz) = s_sort(i50);
        Panis(iz) = p_sort(i50);

    end

    model.Z = Z;
    model.Vs = final_model.Vsav;
    model.Vp = final_model.Vpav;
    model.rho = final_model.rhoav;
    model.Sanis = Sanis';
    model.Panis = Panis';
    par_mineos = struct('R_or_L','Ray','phV_or_grV','phV','ID','startA');

    %Input phv and periods from final model

    [phV,grV,eigfiles] = run_mineos(model,final_predata.SW_Ray_phV.periods,par_mineos,0,0,0);
    K = run_kernels(final_predata.SW_Ray_phV.periods,par_mineos,eigfiles,1,0,0);


    %% Plot

    f1 = figure(1);
    set(f1,'position',[613 82 1308 896]);
    clf

    ax1 = axes(f1,'position',[0.07 0.1 0.5 0.85]);hold on;
    box on;
%     ax2 = axes(f1,'position',[0.31 0.1 0.18 0.85]);hold on;
%     box on;
%     ax3 = axes(f1,'position',[0.55 0.1 0.18 0.85]);hold on;
%     box on;
%     ax4 = axes(f1,'position',[0.79 0.1 0.18 0.85]);hold on;
%    box on;
    for ik = 1:length(K)

        plot(ax1,K{ik, 1}.Vsv + K{ik, 1}.Vsh ,K{ik, 1}.Z./1e3,'linewidth',2);
        
        %plot(ax2,K{ik, 1}.Vsh,K{ik, 1}.Z./1e3,'linewidth',2);
        %plot(ax3,K{ik, 1}.Vpv,K{ik, 1}.Z./1e3,'linewidth',2);
        %plot(ax4,K{ik, 1}.Vph,K{ik, 1}.Z./1e3,'linewidth',2);
    end
    set(ax1,'ydir','reverse','fontsize',24,'linewidth',2);
    %set(ax2,'ydir','reverse','fontsize',24,'linewidth',2);
    %set(ax3,'ydir','reverse','fontsize',24,'linewidth',2);
    %set(ax4,'ydir','reverse','fontsize',24,'linewidth',2);
    title(ax1,'Vs Sensitivity, km/s','fontsize',22);
    ylabel(ax1,'Depth, km','fontsize',22);
    xlabel(ax1,'$\frac{\% \delta Phv}{\% \delta Vs}/$ m','fontsize',28,'interpreter','latex','fontweight','bold')
    %title(ax2,'Vsh Sensitivity, km/s','fontsize',22);
    %title(ax3,'Vpv Sensitivity, km/s','fontsize',22);
    %title(ax4,'Vph Sensitivity, km/s','fontsize',22);
    legend(ax1,{num2str(final_predata.SW_Ray_phV.periods)},'location','southeast','fontsize',24);
    ylim(ax1,[0 300]);
    %ylim(ax2,[0 300]);
    %ylim(ax3,[0 300]);
    %ylim(ax4,[0 300]);
    
save2jpg(f1,['Sensitivity_Kernel_',ddir(id).name(1:4)]);


store(id).K = K;
end 