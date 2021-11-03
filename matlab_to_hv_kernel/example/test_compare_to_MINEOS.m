clear 
% load TRUEmodel.mat
load examplemodel.mat
swperiods = [16,20,24,28,32,36,40,50,60,70,80,90]'; Np = length(swperiods);
TRUEmodel.Vs = TRUEmodel.VS; % It's not clear now if we use VS or vs... will fix this later. bb2021.10.04
TRUEmodel.Vp = TRUEmodel.VP; 
% HVtool
[HVr,HVK_1,phV_1,grV_1] = run_HVkernel(TRUEmodel,swperiods,'test1',1);
% MINEOS

%%% bb2021.11.03 This doesn't seem to be working below right now. 
% [phV_2,grV_2] = run_mineos(TRUEmodel,swperiods,'R','test1',0);
par_mineos = struct('R_or_L', 'R', 'ID', 'test1'); 
[phV_2,grV_2] = run_mineos(TRUEmodel,swperiods,par_mineos);
% [phVkernels_2] = run_kernels(swperiods,'R','ph','test1',1,0);
[phVkernels_2] = run_kernels(swperiods,par_mineos);


% compare

%% Compare phV, grV
figure(11), clf,
hold on
plot(swperiods,[phV_1(:),phV_2(:)],'linewidth',2)
plot(swperiods,[grV_1(:),grV_2(:)],'linewidth',2)

%% Report the differences 
disp('****************')
disp('Comparing phase and group velocities from Mineos and HV code. ')
disp('Using periods:') 
swperiods
disp('(phv_HVkernel - phv_mineos)/phv_HVkernel, at each period, =:')
(phV_1 - phV_2)./phV_1
disp('(grv_HVkernel - grv_mineos)/phv_HVkernel, at each period, =:')
(grV_1 - grV_2)./phV_1
disp('If these are sufficiently similar, you are maybe good to go?')
disp('****************')

%% Compare kernels
figure(88), clf; set(gcf,'pos',[680 385 999 713]);
ax1 = subplot(1,2,1); cla, hold on;
ax2 = subplot(1,2,2); cla, hold on;   
h = {}; % zeros(Np,1);
for ip = 1:Np
        plot(ax1,HVK_1{ip}.Kph_Vs/1e3,0.5*(HVK_1{ip}.Z1+HVK_1{ip}.Z2),'--','linewidth',2);
        plot(ax2,HVK_1{ip}.Kph_Vp/1e3,0.5*(HVK_1{ip}.Z1+HVK_1{ip}.Z2),'--','linewidth',2)        

        h{ip}=plot(ax1,phVkernels_2{ip}.Vsv+phVkernels_2{ip}.Vsh,phVkernels_2{ip}.Z/1e3,'linewidth',2);
        plot(ax2,phVkernels_2{ip}.Vpv+phVkernels_2{ip}.Vph,phVkernels_2{ip}.Z/1e3,'linewidth',2)        
        labstrs = ''; labstrp = '';   
end
%% axes properties
set(ax1,'ydir','reverse','ylim',[0 400],'fontsize',16)
set(ax2,'ydir','reverse','ylim',[0 400],'fontsize',16)
xlabel(ax1,'Sensitivity to $V_S$',...
    'interpreter','latex','fontsize',22)
title(ax1,labstrs);
xlabel(ax2,'Sensitivity to $V_P$',...
    'interpreter','latex','fontsize',22)
title(ax2,labstrp);
ylabel(ax1,'Depth (km)','interpreter','latex','fontsize',22)
ylabel(ax2,'Depth (km)','interpreter','latex','fontsize',22)
%% legend
% hl = legend(ax1,h,num2str(round(swperiods)),'location','southeast');
% set(get(hl,'title'),'string','Period (s)','fontsize',18,'fontweight','bold')
% set(hl,'fontsize',16,'Linewidth',2)
