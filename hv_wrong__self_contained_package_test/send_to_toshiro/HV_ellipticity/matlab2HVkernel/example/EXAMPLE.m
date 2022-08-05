% must be in this directory to run
clear 
load examplemodel.mat
swperiods = [16,20,24,28,32,36,40,50,60,70,80,90]'; Np = length(swperiods);
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
