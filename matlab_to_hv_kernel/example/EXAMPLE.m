% Example of running HV kernels. You need the Fortran codes compiled. 
% must be in this directory to run
clear 
addpath('../'); 
run('../../a0_STARTUP_BAYES.m'); 
load examplemodel.mat
TRUEmodel.Vs = TRUEmodel.VS; % It's not clear now if we use VS or vs... will fix this later. bb2021.10.04
TRUEmodel.Vp = TRUEmodel.VP; 

% Temporary. Re-create model from Tanimoto et al. 2009. 
Z1 = TRUEmodel.Z < 15; 
Z3 = TRUEmodel.Z > 30; 
Z2 = and(~Z1 , ~Z3); 
TRUEmodel.VS(Z1) = 3.0; 
TRUEmodel.VP(Z1) = 5.5; 
TRUEmodel.rho(Z1) = 2.5; 
TRUEmodel.VS(Z2) = 3.6; 
TRUEmodel.VP(Z2) = 6.5; 
TRUEmodel.rho(Z2) = 2.85; 
TRUEmodel.VS(Z3) = 4.5; 
TRUEmodel.VP(Z3) = 7.8; 
TRUEmodel.rho(Z3) = 3.3; 


swperiods = [16,20,24,28,32,36,40,50,60,70,80,90]'; 
swperiods = [10]'; 
Np = length(swperiods);
% HVtool
[HVr,HVK,phV,grV] = run_HVkernel(TRUEmodel,swperiods,'test1',1,1,1);

figure(11);clf;
hp = plot(swperiods,[phV,grV],'o-','Linewidth',2);
set(gca,'fontsize',16)
hl=legend(hp,'Phase Velocity','Group Velocity','location','southeast'); set(hl,'fontsize',15)
xlabel('Period (s)','fontsize',18,'fontweight','bold');
ylabel('Velocity (km/s)','fontsize',18,'fontweight','bold');

figure(12);clf;
hp = plot(swperiods,HVr,'o-','Linewidth',2);
set(gca,'fontsize',16)
xlabel('Period (s)','fontsize',18,'fontweight','bold');
ylabel('H/V ratio','fontsize',18,'fontweight','bold');

exportgraphics(figure(3), 'sens_kers_HV.pdf')

