%% Script to predict VS(T,z) to match the clustered Vs profiles
%% from Brennan's three grouped areas (craton, Appalachians, Grenville)
clear all
close all

%% parameters of atten/V
f = 0.1;
d = 1e-3;
ifsave = true;
mantle_mineralogy = 'lherzolite';

ifmelt = false; % choice to use dry melting + poroelastic effect...
% note, assume melt has ONLY the alpha effect on anelasticity; no viscosity
% jump

anelmodel = 'takei';

%grav
g = 9.81;

%% parameters of P,T structure
dP = 0.1;           % pressure steps (GPa)
P = [1.2:dP:10.5]';% pressure range (GPa)
T = [300:20:1700]';% temperature range (C)

Nt = length(T);
Np = length(P);

[PP,TT] = meshgrid(P,T);

%% paths
run('~/Dropbox/MATLAB/vbr-stable/vbr_init.m')
addpath('~/Dropbox/MATLAB/seis_tools/ABERSHACKER16/')
addpath('/Users/zeilon/Dropbox/Work/Cratons') % for geotherm calculator
nasdatapath = '/Volumes/data/';


% Load mineral database
[minpropar, compar]=ah16_loaddb('AbersHackerMacroJan2016.txt');
 
% lithospheric mantle
switch mantle_mineralogy
    case 'lherzolite' % mantle = lherzolite
        ma_mins ={'fo' ,'fa','en' ,'fs','di','sp'}; 
        ma_modes=[45.51,5.07,22.48,2.5,20.03,4.5];       % These are volume fraction modes
        
    case 'pyrolite' % mantle = lherzolite
        ma_mins ={'alm' ,'gr','py' ,'fo','fa','en','fs','di','hed'}; 
        ma_modes=[2.1   ,1.1 ,10.8 ,54.5,6.7 ,14.7,1.5 ,7.4 ,1.2];       % These are volume fraction modes
        
    case 'harzburgite' % mantle = harzburgite
        ma_mins ={'fo', 'fa', 'en','fs'}; 
        ma_modes=[72.48,7.52,18.24,1.76];       % These are volume fraction modes    
end

%% calculate anharmonic velocities
Vs_anh = nan(Nt,Np);
Vp_anh = nan(Nt,Np);
Gu     = nan(Nt,Np);
rho    = nan(Nt,Np);
fprintf('Getting anharmonic velocities... ')
for it = 1:Nt
for ip = 1:Np
    modu=ah16_rockvel(T(it),P(ip), minpropar, ma_mins,ma_modes);
    % grab rock densities, velocities and unrelaxed moduli
    Vs_anh(it,ip) = modu.vs*1e3;
    Vp_anh(it,ip) = modu.vp*1e3;
    Gu(it,ip)     = modu.g*1e9;
    rho(it,ip)    = modu.r;
end
end
fprintf(' done.\n')

%% calculate effects of melt
phis = zeros(Nt,Np);
Tm_K= nan(Nt,Np);
C_h2o_solid = 50;% ppm
C_co2_solid = 0;
max_phi_wtpct = 5; % maximum melt fraction, in wt%. Assume that any more goes away
if ifmelt
    fprintf('Calculating melt fractions and solidi, poroelastic effect... ')
    Vs_anh_m = nan(Nt,Np);
    Vs_anh_nom = Vs_anh;
    Gu_nom = Gu;

    meltstr = '_wmelt';
    for ip = 1:Np
    for it = 1:Nt
        melt_solution = solve_for_stable_melt_fraction_ze(T(it)+273,P(ip),C_h2o_solid,C_co2_solid);
        Tm_K(it,ip) = melt_solution.Tsolidus_K;
        phis(it,ip) = melt_solution.phi;
        if phis(it,ip) >max_phi_wtpct/100, phis(it,ip) = max_phi_wtpct/100; end
        [ ~, Vs_anh_m(it,ip) ] = melt_poroelastic_Takei( Vp_anh(it,ip),Vs_anh(it,ip),phis(it,ip) );
        if ~isreal(Vs_anh_m(it,ip)),keyboard; end
    end
    end
    % calculate new Gu, with poroelastic effect
    Gu = Vs_anh_m.^2 .* rho;
    fprintf(' done.\n')

else
    meltstr = '_nomelt';
end
alpha_melt = 26; %[Mei et al., 2002]
vfac = exp(-alpha_melt*phis); %cf. Holtzman 2016. Assume no c1 term



%% calculate anelastic velocities
Vs_ane = nan(Nt,Np);
tau_M  = nan(Nt,Np);
fprintf('Calculating anelastic velocities...')
for it = 1:Nt
for ip = 1:Np
switch anelmodel
    case 'takei'
    % takei
    [ ~,~,tau_M(it,ip),~,G_r] = takei2017( T(it)+273,d,P(ip),2*pi*f,vfac(it,ip),Gu(it,ip)/1e9,0);
    Vs_ane(it,ip) = sqrt(G_r./rho(it,ip));
    case 'JF10'
    % JF10
    [ ~,~,tau_M(it,ip),~,G_r] = JF2010( T(it)+273,d,P(ip),2*pi*f,vfac(it,ip),Gu(it,ip)/1e9 );
    Vs_ane(it,ip) = sqrt(G_r./rho(it,ip));
    case 'PM2013'
    % PM13
    [ ~,~,tau_M(it,ip),~,G_r] = PM2013( T(it)+273,d,P(ip),2*pi*f,vfac(it,ip),Gu(it,ip)/1e9 );
    Vs_ane(it,ip) = sqrt(G_r./rho(it,ip));
    case 'FJ05'
    % fj05
    [ ~,~,tau_M(it,ip),~,G_r] = FJ2005( T(it)+273,d,P(ip),2*pi*f,vfac(it,ip),Gu(it,ip)/1e9 );
    Vs_ane(it,ip) = sqrt(G_r./rho(it,ip));
end % switch
end % P loop
end % T loop
fprintf(' done.\n')


% calculate eta
G = Vs_ane.^2.*rho;
eta_log10 = log10(G.*tau_M.*vfac);

%% little plotting
if ifmelt
figure(672),clf,set(gcf,'Position',[23 164 1359 402])
%anharmonic velocity, no melt effect
subplot(131)
contourf(TT,PP,Vs_anh_nom,linspace(3600,5000,35))
colormap(gca,flipud(turbo)),caxis([3800 5000])
hcb = colorbar;
xlabel('T (C)')
ylabel('P (GPa')
ylabel(hcb,'Vs-anh-nom (m/s)')
%anharmonic velocity, WITH melt effect
subplot(132)
contourf(TT,PP,Vs_anh_m,linspace(3600,5000,35))
colormap(gca,flipud(turbo)),caxis([3800 5000])
hcb = colorbar;
xlabel('T (C)')
ylabel('P (GPa')
ylabel(hcb,'Vs-anh-wm (m/s)')
%anelastic velocity
subplot(133)
contourf(TT,PP,Vs_ane,linspace(3600,5000,35))
colormap(gca,flipud(turbo)),caxis([3800 5000])
hcb = colorbar;
xlabel('T (C)')
ylabel('P (GPa')
ylabel(hcb,'Vs-ane (m/s)')
end


figure(673),clf,set(gcf,'Position',[23 464 1059 402])
%velocity
subplot(121)
contourf(TT,PP,Vs_ane,linspace(3600,5000,35))
colormap(gca,flipud(turbo)),caxis([3800 5000])
hcb = colorbar;
xlabel('T (C)')
ylabel('P (GPa')
ylabel(hcb,'Vs (m/s)')
% viscosity
subplot(122)
contourf(TT,PP,eta_log10,linspace(17,23,35))
colormap(gca,flipud(parula)),caxis([18 23])
hcb = colorbar;
xlabel('T (C)')
ylabel('P (GPa')
ylabel(hcb,'log_{10}\eta (Pa s)')

if ifsave
    save2pdf(673,['Vs_T_P_',anelmodel,meltstr]);
end

%% save
if ifsave
    ofile = ['VsTP_lookup_',anelmodel,meltstr];
    save(ofile,'T','P','Vs_ane',"tau_M",'d','f','mantle_mineralogy','anelmodel','rho','eta_log10','phis');
end


return