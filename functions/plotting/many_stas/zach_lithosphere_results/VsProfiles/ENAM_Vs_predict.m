%% Script to predict VS(T,z) to match the clustered Vs profiles
%% from Brennan's three grouped areas (craton, Appalachians, Grenville)
% clear all
% close all

%% parameters of atten/V
f = 0.1;
d = 1e-3;

%grav
g = 9.81;

%% parameters of thermal structure
dz = 1000;       % depth steps (m)
z = [0:dz:300*1e3]';% depth (m)

%% paths
run('~/Dropbox/MATLAB/vbr-stable/vbr_init.m')
addpath('~/Dropbox/MATLAB/seis_tools/ABERSHACKER16/')
addpath('/Users/zeilon/Dropbox/Work/Cratons') % for geotherm calculator
nasdatapath = '/Volumes/data/';


% for heat flow, see
% https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1002/ggge.20271
% Davies, J. H. (2013), Global map of solid Earth surface heat flow,
% Geochem. Geophys. Geosyst., 14, 4608– 4622, doi:10.1002/ggge.20271.

%% loop over structures - aim to step from E to W (young to old)
run_names = {'pied','gren','crat'};
z_mohs   = [ 35       ,43        ,41   ]*1e3; % moho depth (km)
q_0s     = [ 48       ,45        ,40   ]*1e-3; % surface heat flux (W/m^2)
T_ps     = [ 1350     ,1350      ,1350 ]; % Potential temperatures (C)
A0_ucs   = [ 1.3      ,1.4       ,1.17  ]*1e-6; % volumetric heat production, upper crust (W/m^3)

z_mohs   = [ 25       ,40        ,50   ]*1e3; % moho depth (km)
q_0s     = [ 48       ,45        ,40   ]*1e-3; % surface heat flux (W/m^2)
T_ps     = [ 1350     ,1350      ,1350 ]; % Potential temperatures (C)
A0_ucs   = [ 1.3      ,1.4       ,1.2  ]*1e-6; % volumetric heat production, upper crust (W/m^3)

lithmantle_mineralogy = 'lherzolite';

%% defaults (some will change)
  par = struct(...
    'z_ucb',17e3,...   % upper/lower crustal boundary(m)
    'z_moh',30e3,...   % moho depth (m)
    'z_lab',nan,... % lab depth (m) < will come out of top-down inversion
    'q_0',84e-3,...  % Surface heat flux (W/m^2)
    'q_m',nan,...    % Moho heat flux (W/m^2), often 30e-3; < will come out of top-down inversion
    'T_0',10,...       % Z=0 temperature (C)
    'k0',4.55,...     % ref thermal conductivity in mantle(W/m/K)
    'kav_m',nan,...  % average mantle thermal conductivity in mantle(W/m/K), often 3.35 
    'k0_uc',2.5,...     % thermal conductivity in crust(W/m/K)
    'k0_lc',2.5,...     % thermal conductivity in crust(W/m/K)
    'A0_uc',2.3e-6,... % volumetric heat production, upper crust (W/m^3)
    'A0_lc',0.8e-6,... % volumetric heat production, lower crust (W/m^3)
    'dA_lc',10e3,... % e-folding length of A in lower crust (m)
    'A_m',0.01e-6,... 0.01e-6,... % volumetric heat production, mantle (W/m^3)
    'q_lab',10e-3,...  % lab heat flux (W/m^2)
    'dTdz_ad',0.35e-3,...     % adiabatic gradient (K/m)
    'T_p',1350,...      % potential temperature (C)
    'Rho_oma',3.25); % density of mantle at reference temperature (0 C)



%% paths
addpath('~/Dropbox/MATLAB/seis_tools/anelasticity/ANEL_MODELS/')

for imod = 1:length(run_names)
    % parameters for this run
    par.run_name = run_names{imod};
    par.z_moh = z_mohs(imod);
    par.q_0 = q_0s(imod);
    par.T_p = T_ps(imod);
    par.A0_uc = A0_ucs(imod);

    %% get temperature
    [T,q,z_lab] = geotherm_continent_fn(z,par,'top_down',1);
   
    if isempty(z_lab)||isnan(z_lab),z_lab=max(z); end
    
    figure(33),clf
    subplot(121),plot(T,z/1e3,'linewidth',2),set(gca,'ydir','rev');xlabel('T (C)');ylabel('Depth (km)')
    subplot(122),plot(q*1e3,z/1e3,'linewidth',2),set(gca,'ydir','rev');xlabel('q (mW/m^2)');

    %% Density and pressure
    iuc = find(z <= par.z_ucb);
    ilc = find(par.z_ucb <= z & z <= par.z_moh);
    iml = find(par.z_moh <= z & z <= z_lab);  % mantle lith
    ima = find(z_lab <= z);                   % mantle asth
    % prelim density and pressure
    rho = nan(size(T));
    rho(iuc) = 2650;
    rho(ilc) = 2750;
    rho(iml) = 3300*(1 - 4e-5*T(iml));
    rho(ima) = 3300*(1 - 4e-5*T(ima));
    P = cumtrapz(z,rho*9.81);
    
    % Load mineral database
    [minpropar, compar]=ah16_loaddb('AbersHackerMacroJan2016.txt');
    
    % nth pass - (first ignore P,T sensitivity, then use it
    for kk = 1:4
        fprintf('Loop %u/5 for pressure...',kk)
        clear modu
        % upper crust - from New Theory of the Earth
        uc_mins ={'qz' ,'or','lAb','an','phl','mu','hb','en','fs','fo'}; 
        uc_modes=[ 24  , 12 , 20  , 19 , 4   , 1 , 5   , 5  , 6  , 3];       % These are volume fraction modes
        for k = 1:length(iuc)
            modu(iuc(k),1)=ah16_rockvel(T(iuc(k)),P(iuc(k))/1e9, minpropar, uc_mins,uc_modes);
        end
        % lower crust - granulite
        lc_mins ={'hAb' ,'an','en','fs','di','hed','tr','fact','ts','parg','mt'}; 
        lc_modes=[ 21   , 21 , 9  , 9  , 5  , 9   , 2  , 7    , 9  , 7    , 1];       % These are volume fraction modes
        for k = 1:length(ilc)
            modu(ilc(k),1)=ah16_rockvel(T(ilc(k)),P(ilc(k))/1e9, minpropar, lc_mins,lc_modes);
        end
        % lithospheric mantle
        switch lithmantle_mineralogy
            case 'lherzolite' % mantle = lherzolite
                ma_mins ={'fo' ,'fa','en' ,'fs','di','sp'}; 
                ma_modes=[45.51,5.07,22.48,2.5,20.03,4.5];       % These are volume fraction modes
                for k = 1:length(iml)
                    modu(iml(k),1)=ah16_rockvel(T(iml(k)),P(iml(k))/1e9, minpropar, ma_mins,ma_modes);
                end
            case 'pyrolite' % mantle = lherzolite
                ma_mins ={'alm' ,'gr','py' ,'fo','fa','en','fs','di','hed'}; 
                ma_modes=[2.1   ,1.1 ,10.8 ,54.5,6.7 ,14.7,1.5 ,7.4 ,1.2];       % These are volume fraction modes
                for k = 1:length(iml)
                    modu(iml(k),1)=ah16_rockvel(T(iml(k)),P(iml(k))/1e9, minpropar, ma_mins,ma_modes);
                end
            case 'harzburgite' % mantle = harzburgite
                ma_mins ={'fo', 'fa', 'en','fs'}; 
                ma_modes=[72.48,7.52,18.24,1.76];       % These are volume fraction modes
                for k = 1:length(iml)
                    modu(iml(k),1)=ah16_rockvel(T(iml(k)),P(iml(k))/1e9, minpropar, ma_mins,ma_modes);
                end
        end
        % asthenosphereic mantle
        %(pyrolite)
%         ma_mins ={'alm' ,'gr','py' ,'fo','fa','en','fs','di','hed'}; 
%         ma_modes=[2.1   ,1.1 ,10.8 ,54.5,6.7 ,14.7,1.5 ,7.4 ,1.2];       % These are volume fraction modes
        %(lherzolite)
        ma_mins ={'fo' ,'fa','en' ,'fs','di','sp'}; 
        ma_modes=[45.51,5.07,22.48,2.5,20.03,4.5];       % These are volume fraction modes
        for k = 1:length(ima)
            modu(ima(k),1)=ah16_rockvel(T(ima(k)),P(ima(k))/1e9, minpropar, ma_mins,ma_modes);
        end
        
        rho = [modu.r]';
        P_ = cumtrapz(z,rho*9.81);
        fprintf('. average update %.1f MPa\n',rms(P_-P)/1e6);
        P = P_;
    end
    % grab rock densities, velocities and unrelaxed moduli
    Vp = [modu.vp]';
    Vs_anh(:,imod) = [modu.vs]';
    rho = [modu.r]';
    Gu = [modu.g]'*1e9;
    P_base(imod,1) = P(end);
    
    T_plot(:,imod) = T;
    
    %% predict velocities (and Q)
    % anharmonic (unrelaxed)
%     par.muU0       = 72.45;      % GPa unrelaxed compliance at 0 Celsius, 0 Pa
%     par.dmuU_dT    = -1.094e-2;   % GPa/K compliance variation with T
%     par.dmuU_dP    = 1.987;      % dim'less compliance variation with P
%     Gu = par.muU0 + par.dmuU_dT*T + par.dmuU_dP*(P/1e9);

    % grain size shenanigans
    d = d.*ones(length(z),1);
% %     d(iml) = 1e-4;

    % melt shenanigans
    vfac = ones(length(z),1);
%     phis = zeros(length(z),1);
%     phis(z<100e3) = 0.1;
%     vfac =  melt_visc_reduce_TH09(phis)/8;
    

    % takei
    for iz = 1:length(z)
        [ J1(iz,1),J2(iz,1),tauM(iz,1),qinv_takei(iz,1),G_takei(iz,1)] = takei2017 ( T(iz)+273,d(iz),P(iz)/1e9,2*pi*f,1,Gu(iz)/1e9 );
    end
    VS_takei(:,imod) = sqrt(G_takei./rho);
    Q_takei(:,imod) = 1./qinv_takei; Q_takei(Q_takei>1e3) = 1e3;

    % JF10
%     for iz = 1:length(z)
%         [ J1(iz,1),J2(iz,1),tauM(iz,1),qinv_jf10(iz,1),G_jf10(iz,1)] = JF2010 ( T(iz)+273,d(iz),P(iz)/1e9,2*pi*f,vfac(iz),Gu(iz) );
%     end
%     VS_jf10(:,imod) = sqrt(G_jf10./rho);
%     Q_jf10(:,imod) = 1./qinv_jf10; Q_jf10(Q_jf10>1e3) = 1e3;

    % PM13
    for iz = 1:length(z)
        [ J1(iz,1),J2(iz,1),tauM(iz,1),qinv_pm13(iz,1),G_pm13(iz,1)] = PM2013 ( T(iz)+273,d(iz),P(iz)/1e9,2*pi*f,vfac(iz),Gu(iz)/1e9 );
    end
    VS_pm13(:,imod) = sqrt(G_pm13./rho);
    Q_pm13(:,imod) = 1./qinv_pm13; Q_pm13(Q_pm13>1e3) = 1e3;

    % fj05
    for iz = 1:length(z)
        [ J1(iz,1),J2(iz,1),tauM(iz,1),qinv_fj05(iz,1),G_fj05(iz,1)] = FJ2005 ( T(iz)+273,d(iz),P(iz)/1e9,2*pi*f,vfac(iz),Gu(iz) );
    end
    VS_fj05(:,imod) = sqrt(G_fj05./rho);
    Q_fj05(:,imod) = 1./qinv_fj05; Q_fj05(Q_fj05>1e3) = 1e3;

end


%% PLOTTING
figure(55);clf;set(gcf,'pos',[244 259 756 540]);
ax1 = axes('pos',[0.12 0.11 0.34 0.84]); hold on
ax2 = axes('pos',[0.53 0.11 0.39 0.84]); hold on

% Velocity
for imod = 1:length(run_names)
    if imod > size(Vs_anh,2), continue; end
    if isempty(Vs_anh(1,imod)), continue; end
%     if imod==1, clr = 'r'; elseif imod==2, clr ='b'; end
    clr = colour_get(imod,length(run_names)+1,0,parula);
    plot(ax2,Vs_anh(:,imod)*1e3,z,'linestyle','-','linewidth',1.5,'color','k'); 
    plot(ax2,VS_takei(:,imod),z,'linestyle','--','linewidth',1.5,'color',clr); 
%     plot(ax2,VS_jf10(:,imod),z,'linestyle','-','linewidth',1.5,'color',clr); 
%     plot(ax2,VS_pm13(:,imod),z,'linestyle','-.','linewidth',1.5,'color',clr); 
%     plot(ax2,VS_fj05(:,imod),z,'linestyle',':','linewidth',1.5,'color',clr); 
end

% Quality
% plot(ax3,Q_takei,z,'linestyle','--','linewidth',1.5); 
% plot(ax3,Q_jf10,z,'linestyle','--','linewidth',1.5); 
% plot(ax3,Q_pm13,z,'linestyle','-.','linewidth',1.5); 
% plot(ax3,Q_fj05,z,'linestyle',':','linewidth',1.5); 
% set(ax3,'ydir','reverse','xlim',[0 400],'ylim',[0 max(z)])

% Temp
for imod = 1:length(run_names)
    if imod > size(Vs_anh,2), continue; end
    if isempty(Vs_anh(1,imod)), continue; end
%     if imod==1, clr = 'r'; elseif imod==2, clr ='b'; end
    clr = colour_get(imod,length(run_names)+1,0,parula);
    plot(ax1,T_plot(:,imod),z,'linewidth',1.5,'color',clr); 
end

% pretty
set(ax1,'ydir','reverse','ylim',[0 max(z)],...
    'linewidth',2,'layer','top','box','on',...
    'ytick',[0:50:300]*1e3,'yticklabels',[0:50:300],'fontsize',17)
set(ax2,'ydir','reverse','xlim',[3.6 4.9]*1e3,'ylim',[0 max(z)],...
    'linewidth',2,'layer','top','box','on',...
    'yticklabels',[],'fontsize',17)


xlabel(ax1,'\textbf{Temp.} ($^{\circ}$C)','fontsize',20,'interpreter','latex')
xlabel(ax2,'$\mathbf{V_S}$ (km/s)','fontsize',20,'interpreter','latex')
ylabel(ax1,'\textbf{Depth} (km)','fontsize',20,'interpreter','latex')

% legend
% hl = legend(ax2,{'Takei17','JF10','PM13','FJ05'},'location','southwest');

%% plot comparison with observed clustered profiles
load('enam_mcmc_average_vel_profiles.mat');
h_inv = plot(ax2,profiles_export*1e3,z_export*1e3,'-k','linewidth',3);
eus_cls = {[0.3 0.9 0.4],[0.7 0.1 0.6],[0.7 0.6 0.1]};
eus_rgns = {'margin','craton','gren'};
ylim([ax2,ax1],[0 270]*1e3)

%% retrieve SEMUCB global comparison models
addpath([nasdatapath,'models_seismic/SEMum2_avg_VS/']);
a = SEMum2_avgprofiles(0,[nasdatapath,'models_seismic/SEMum2_avg_VS/']);
% h_sem_yoce = plot(ax2,a.ocean_0_25*1e3,a.Z*1e3,':c','linewidth',2);
h_sem_phan = plot(ax2,a.cont_Phan*1e3,a.Z*1e3,':','linewidth',2,'color',[188,143,143]./255);
h_sem_arch = plot(ax2,a.cont_Arch*1e3,a.Z*1e3,'--','linewidth',2,'color',[188,143,143]./255);

%% predicted differential topography
fprintf('predicted differential topography (negative means in the air)\n')
(P_base - mean(P_base))/rho(1)/1000

%% predicted LABs - second derivative messing around
figure(99); clf,hold on
for ip = 1:3
    fprintf('--------------\n\n%s profile\n',eus_rgns{ip})
    zmoh(ip) = z_export(find(profiles_export(:,ip)>4.2,1,'first'));

    [a,b] = LAB_finder(profiles_export(:,1),z_export,zmoh(ip),1)

    smthz = moving_average(profiles_export(:,ip),10,1);
    d2vdz1 = gradient(smthz);
    d2vdz2 = gradient(gradient(smthz));

    subplot(311), hold on, 
    plot(z_export,d2vdz2,'linewidth',2,'color',eus_cls{ip})
    xlim([50 270]),ylim([-0.0002 0.0002]), yline(0)
    xlabel('depth, km'),ylabel('second derivative')

    subplot(312), hold on
    plot(z_export,d2vdz1,'linewidth',2,'color',eus_cls{ip})
    xlim([50 270]),ylim([-0.002 0.002]), yline(0)
    xlabel('depth, km'),ylabel('first derivative')

    subplot(313), hold on
    plot(z_export,smthz,'linewidth',2,'color',eus_cls{ip})
    xlim([50 270])
    xlabel('depth, km'),ylabel('velocity (km/s)')

    % method 1: stationary points of first derivative
    invg = crossing(d2vdz2); % find maximum gradient points (turning)
    invg(d2vdz1(invg)>0) = []; % have to be negative velocity gradients
    invg(z_export(invg)>250) = []; % no nvg deeper than 250
    invg(z_export(invg)<zmoh(ip)) = []; % no nvg deeper than 250
    [~,imax] = max(smthz(1:200));

    if any(invg<imax) 
        zmld(ip) = z_export(invg(find(invg<imax,1,'first')));
    else
        zmld(ip) = nan;
    end    

    invg(invg<imax) = [];% has to be deeper than the Vs maximum
    invg = invg(mindex(d2vdz1(invg))); % choose the one associated with steepest nvg
    zlab(ip) = z_export(invg(find(invg>imax,1,'first')));

    
    % method 2: velocity minimum
    [itrn,ztrn] = crossing(d2vdz1,z_export); % find maximum gradient points (turning)
    ztrn(z_export(itrn)<zlab(ip)) = []; % must be deeper than method 1 lab
    itrn(z_export(itrn)<zlab(ip)) = []; % must be deeper than method 1 lab
    zlab2(ip) = ztrn(find(itrn>imax,1,'first'));

    % method 3: halfway between max and min points
    itrn = crossing(d2vdz1); % find maximum gradient points (turning)
    itrn(z_export(itrn)<zmoh(ip)) = []; % must be in mantle
    itrn(itrn < (imax-1)) = [];
    imin = itrn(end);
    if length(itrn)>2, itrn = itrn(1:2);end
    zlab3(ip) = 0.5*diff(z_export(itrn)) + z_export(imax);
    
    % method 4: 90% between max and min points
    [itrn0,ztrn] = crossing(d2vdz1,z_export); % find maximum gradient points (turning)
    itrn = itrn0;
    itrn(z_export(itrn)<zmoh(ip)) = []; % must be in mantle
    itrn(itrn < (imax-1)) = [];
    V_mami = profiles_export(itrn,ip);
    [~,zzo] = crossing(profiles_export(:,ip),z_export,V_mami(2)+0.1*(max(V_mami)-min(V_mami)));
    zzo(zzo<=z_export(imax)) = [];
    zzo(zzo>=z_export(imin)) = [];
    zlab4(ip) = zzo(1);

    fprintf('    Z_lab = %.0f km, Z_mld = %.0f km\n',zlab(ip) ,zmld(ip))
    fprintf('OR  Z_lab = %.0f km\n',zlab2(ip))
    fprintf('OR  Z_lab = %.0f km\n',zlab3(ip))
    fprintf('OR  Z_lab = %.0f km\n',zlab4(ip))
end
legend(eus_rgns)
% xlim([50 250]),ylim([-0.0002 0.0002])
return
save


% %% plot stations
% for imod = 1:length(run_names)
%     if imod > size(Vs_anh,2), continue; end
%     if isempty(Vs_anh(1,imod)), continue; end
%     if Vs_anh(1,imod)==0, continue; end
% 
%     switch run_names{imod} % {'young_phaneroz','archean1','K22A-WY','WVOR'}
%         case 'archean1'
%             stanw = 'SCIA_US';
%         case 'K22A-WY'
%             stanw = 'K22A_TA';
%         case 'WVOR'
%             stanw = 'WVOR_US';
%     end
% 
%     inv_mod = load([nasdatapath,'THBI/NWUS/STASinv/',stanw,'_dat30/NWUS_CCP_HK/final_model.mat']);
%     inv_mod = inv_mod.final_model;
%     h_inv = plot(ax2,inv_mod.VSav*1e3,inv_mod.Z*1e3,'-k','linewidth',3);
%     h_inv1 = plot(ax2,inv_mod.VSsig1*1e3,inv_mod.Z*1e3,'--','linewidth',2,'color',[0.5 0.5 0.5]);
% 
% end
% return
% 
% % eymn_mod = load([nasdatapath,'BayesianJointInv/NWUS_STASinv/EYMN_US_dat30/NWUS_CCP_HK/final_model.mat');
% % eymn_mod = eymn_mod.final_model;
% % h_eymn = plot(ax2,eymn_mod.VSav*1e3,eymn_mod.Z*1e3,'-r','linewidth',3);
% 
% k22a_mod = load([nasdatapath,'BayesianJointInv/NWUS_STASinv/K22A_TA_dat30/NWUS_CCP_HK/final_model.mat']);
% k22a_mod = k22a_mod.final_model;
% h_k22a = plot(ax2,k22a_mod.VSav*1e3,k22a_mod.Z*1e3,'-m','linewidth',3);
% h_k22a1 = plot(ax2,k22a_mod.VSsig1*1e3,k22a_mod.Z*1e3,'--m','linewidth',2);
% 
% boz_mod = load([nasdatapath,'BayesianJointInv/NWUS_STASinv/BMO_US_dat30/NWUS_CCP_HK/final_model.mat']);
% boz_mod = boz_mod.final_model;
% h_boz = plot(ax2,boz_mod.VSav*1e3,boz_mod.Z*1e3,'-g','linewidth',3);
% 
% 
% boz_mod = load([nasdatapath,'BayesianJointInv/NWUS_STASinv/SCIA_US_dat30/NWUS_CCP_HK/final_model.mat']);
% boz_mod = boz_mod.final_model;
% h_boz = plot(ax2,boz_mod.VSav*1e3,boz_mod.Z*1e3,'-g','linewidth',3);
% h_boz = plot(ax2,boz_mod.VSsig1*1e3,boz_mod.Z*1e3,'--g','linewidth',2);
