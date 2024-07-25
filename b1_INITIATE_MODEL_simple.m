function model = b1_INITIATE_MODEL(par,selev,ifplot,hkStack)
% model = b1_INITIATE_MODEL(par)
% 
% Function to make the nstas x 1-D model as per the parameters defined in
% the input structure par. If relevant, give station elevation in km.
% 
% % model parameters:
%   1. hsed
%   2. vs_sed_top
%   3. vs_sed_bottom
%   4. hcrust
%   5. vpvs_crust (may be fixed)
%   ++ V_knots_crust (there are N_knots_crust of these)
%   ++ V_knots_mantle (there are N_knots_mantle of these)
%
%  ==> 5 + Nkc + Nkm parms 
% 
%       (c.f. 13 for Shen and Ritz, who fixed vpvs and Nkc=4, Nkm=5

%%% From z0_...
if nargin <2 || isempty(selev)
    selev = zeros(par.mod.nstas,1);
end

% figure(1); clf; hold on; 
% set(gca, 'YDir', 'reverse'); 
% ylabel('Depth'); 
p = iasp91(); 
% in_mod = p.depth < (par.mod.maxz+200); 
% names = {'depth', 'vp', 'vs', 'rho'};% , 'qk', 'qu'}; 
% for iname = 1:length(names); 
%     name = names{iname}; 
%     vals = p.(name); 
%     vals = vals(in_mod); 
%     p.(name) = vals; 
% end
% plot(p.vs , p.depth, 'DisplayName', 'Vs' ); 
% plot(p.vp , p.depth, 'DisplayName', 'Vp' ); 
% plot(p.rho, p.depth, 'DisplayName', 'rho'); 

% Crustal thickness. Just pick in the middle of priors. Was 42.5, which is close to average US crustal thickness. brb2024.07.24
h_crust = mean([par.mod.crust.hmax, par.mod.crust.hmin]); 

% Make Vs sediment. Determine the starting values based on par allowed values, and increasing. 
h_sed = 0.5; 
vs_sed_allowance = par.mod.sed.vsmax - par.mod.sed.vsmin; 
vs_sed_mid = mean([par.mod.sed.vsmax, par.mod.sed.vsmin]); 
vs_sed = [vs_sed_mid - vs_sed_allowance/4, vs_sed_mid+vs_sed_allowance/4]; 

% Crust knots. 
kvs_crust = [3.36 3.68 4.0]'; % 3.36 from lowest iasp91 crustal velocity (no ocean). ~4.0 km/s at ~40 km depth in Brunsvik et al 2024. 
cknots = linspace(h_sed, h_sed+h_crust,length(kvs_crust)-1)'; % Not sure why there is 1 less knot than velocity value. 
fcknots = (cknots-h_sed)/h_crust;
k_crust = length(cknots);

% Mantle knots. 
% Remove layers from our starting mantle model. Allows us to use interp1. 
p_d_terp = p.depth; 
p_v_terp = p.vs; 
is_lay = p_d_terp(2:end) == p_d_terp(1:end-1); 
repl_with = (p_v_terp(find(is_lay)) + p_v_terp(find(is_lay)+1))./2; 
p_v_terp(find(is_lay)+1) = repl_with; 
p_v_terp = p_v_terp(~is_lay); 
p_d_terp = p_d_terp(~is_lay); 

% Make depths and interp velocity at ddepths. 
z_terp = linspace(h_sed+h_crust, par.mod.maxz, par.mod.mantle.kmin+par.mod.mantle.kdstd)'; % Depth of moho in iasp91, to max depth. Use an amount of knots such that we might reasonably end up adding or removing knots. 
kvs_mantle = interp1(p_d_terp, p_v_terp, z_terp); 
mknots = linspace(h_sed+h_crust, par.mod.maxz, par.mod.mantle.kmin+par.mod.mantle.kdstd-1)'; 
fmknots = (mknots-h_sed-h_crust)/(par.mod.maxz-h_sed-h_crust);
k_mantle = length(mknots);

% xi and vpvs. 
vpvs_crust = 1.75; 
xi_mantle = ones(size(par.mod.mantle.xidepths)); % 1.0; 
n_xi_mantle = length(xi_mantle); 
xi_crust = 1; 

% DEPTHS
cminz = h_sed;
cmaxz = h_sed+h_crust;
zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';
mminz = h_sed+h_crust;
mmaxz = par.mod.maxz + selev;
zm = unique([mminz:par.mod.dz:mmaxz,mmaxz])';

% Make spline basis. 
[ cspbasis ] = make_splines( cknots,par, zc);
[ mspbasis ] = make_splines( mknots,par, zm);

%% MAKE ALL PARAMETER STRUCTURES
sed = struct('h',h_sed,'VS',vs_sed);
crust = struct('h',h_crust,'Nsp',k_crust+1,'Nkn',k_crust,'VS_sp',kvs_crust,'vpvs',vpvs_crust,'xi',xi_crust,'splines',cspbasis,'knots',cknots,'fknots',fcknots,'z_sp',zc);
mantle = struct('Nkn',k_mantle,'Nsp',k_mantle+1,'VS_sp',kvs_mantle,'xi',xi_mantle,'splines',mspbasis,'knots',mknots,'fknots',fmknots,'z_sp',zm);
data = ([]);

%% DATA PARMS
hparm = struct([]);
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id};
    [ ~,sigma_prior ] = get_sigma_min_prior( dtype,par );
    hparm(1).(['sig_',dtype]) = sigma_prior;
end

%% OVERALL PARMS
mod = par.mod; 
M = 2 + 1 + k_crust + k_mantle ...
      + double(mod.sed.hmin~=mod.sed.hmax) ...
      + double(mod.crust.vpvsmin~=mod.crust.vpvsmax) ...
      + double(mod.crust.ximin~=mod.crust.ximax) ... 
      + double(mod.mantle.ximin~=mod.mantle.ximax) * n_xi_mantle ;

% %% MODEL WITH ALL PARMS
model = struct('sedmparm',sed,'crustmparm',crust,'mantmparm',mantle,'datahparm',hparm,...
               'M',M,'selev',selev);
% 
% %% TURN PARMS INTO REAL MODEL
model = make_mod_from_parms(model,par);

%% plot
if ifplot
    figure(31), clf
    plot(model.VS,model.z,'b',model.VP,model.z,'r',model.rho,model.z,'k')
    set(gca,'ydir','reverse')
end