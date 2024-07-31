function [ model ] = make_mod_from_parms( model,par, options )
    arguments
        model
        par
        options.use_splines = true; 
        options.vs = struct('crust',[],'mantle',[]); 
        options.z  = struct('crust',[],'mantle',[]); 
    end
% [ model ] = make_mod_from_parms( model )
% 
% make the 1D model vectors and important model values from the three
% structures sedparm, crustparm, and mantparm that actually house the
% important values that go towards building the model.

mps = model.sedmparm;
mpc = model.crustmparm;
mpm = model.mantmparm;

%% SEDIMENTS
zs0     = [0,mps.h]'; % This is all that's initially parameterised with sediment. 
vs_sed0 = mps.VS(:);

zs = unique([0:par.mod.dz:mps.h,mps.h])'; % Now up-interpolate to steps of dz
vs_sed = linterp(zs0, vs_sed0, zs); % Now up-interpolate to steps of dz

vp_sed = sed_vs2vp(vs_sed);
rho_sed = sed_vs2rho(vs_sed);
xi_sed = ones(size(zs));

if diff(zs0)==0; zs=[]; vs_sed=[]; vp_sed=[]; rho_sed=[]; xi_sed=[]; end

%% CRUST
if options.use_splines; 
    cminz = mps.h;
    cmaxz = mps.h+mpc.h;
    zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';
    vs_crust = sum(mpc.splines*diag(mpc.VS_sp),2);
else
    vs_crust = options.vs.crust; 
    zc = options.z.crust; 
end
vp_crust = vs_crust*mpc.vpvs; %!%!%! brb2022.06.03 TODO here can make the od vpvs jump at moho work better by tapering the transition of  vpvs crust to vpvs mantle
rho_crust = sed_vs2rho(vs_crust); % use same expression as for sed (assume not ultramafic or other poorly-fit rocks)
xi_crust = mpc.xi*ones(size(zc));


%% MANTLE
% try
if options.use_splines; 
    mminz = cmaxz;
    mmaxz = par.mod.maxz + model.selev; 
    zm = unique([mminz:par.mod.dz:mmaxz,mmaxz])';
    vs_mantle = sum(mpm.splines*diag(mpm.VS_sp),2);
else
    vs_mantle = options.vs.mantle; 
    zm = options.z.mantle; 
end
vp_mantle = mantle_vs2vp(vs_mantle,zm );
rho_mantle = mantle_vs2rho(vs_mantle,zm );

%%% Handle xi through mantle. 
% This commented version used Gaussian interpolation. It seems worse for our purposes than pchip interpolation. 
% % % % Handle xi. Something like splines, but not adjustable in depth. Gaussian interpolation. 
% % % zxi = [par.mod.crust.hmax; 100; 150; 200; max(zm)]; % Depths where we know xi
% % % xi_at_z = [mpm.xi(1); mpm.xi(1); mpm.xi(2); 1; 1]; % Xi values, corresponding to zxi
% % % 
% % % sigma = (zxi(end-1)-zxi(2))/length(zxi(2:end-1)) * 4/5; % Standard deviation for Gaussian weighting. Adjust by a constant at the end for... taste. 
% % % xi_interp = zeros(size(zm)); % Initialize the interpolated xi array
% % % for i = 1:length(zm) % Loop over each zm value to interpolate
% % %     dist = zxi - zm(i); 
% % %     weights = exp(-((dist).^2) / (2 * sigma^2)); % Compute Gaussian weights based on distance to each zxi
% % %     weights = weights / sum(weights); % Normalize the weights
% % %     xi_interp(i) = sum(weights .* xi_at_z); % Compute the weighted average of xi_at_z
% % % end
% % % xi_mantle = xi_interp; 

%Handle xi through mantle using pchip interpolation. 
xidepths = par.mod.mantle.xidepths; 
mpmxi = mpm.xi; 
xidepths(1) = max([xidepths(1), model.crustmparm.h]); % If our shallowest mantle xi value was in the crust, move the depth to the Moho. 

zxi     = [model.crustmparm.h; xidepths+0.00001; 150; max(zm)]; % Depths where we know xi
xi_at_z = [mpmxi(1);           mpmxi           ;      1; 1      ]; % Xi values, corresponding to zxi

% Interpolate to zm
xi_mantle = interp1(zxi, xi_at_z, zm, 'pchip');

% Plotting results
plot_interp = false;
if plot_interp; 
    figure(4001); clf; hold on;
    plot(xi_mantle, zm, 'k'); % Interpolated xi values
    scatter(xi_at_z, zxi, 'filled'); % Known xi values at zxi
    set(gca, 'YDir', 'reverse');
    xlabel('Xi'); ylabel('Depth');
    title('Xi Interpolation');
end
%%% Done with xi. 

%% COLLATE
zz = [zs;zc;zm];
zz0 = zz-model.selev; % true depth, from sea level/ref. ellipsoid
Nz = length(zz);
zsed = mps.h;
zmoh = mps.h+mpc.h;
vs = [vs_sed;vs_crust;vs_mantle];
vp = [vp_sed;vp_crust;vp_mantle];
rho = [rho_sed;rho_crust;rho_mantle];
xi = [xi_sed;xi_crust;xi_mantle];

%% OUTPUT
model.z    = zz;
model.z0   = zz0;
model.VS   = vs;
model.VP   = vp;
model.rho  = rho;
model.Nz   = Nz;
model.zsed = zsed;
model.zmoh = zmoh; % Model.zmoh is absolute moho depth. mps.h + mpc.h
model.vpvs = mpc.vpvs;
model.cxi  = mpc.xi;
model.mxi  = mpm.xi;
model.fdVSsed = 100*diff(vs(zz==zsed))./mean(vs(zz==zsed)); if zsed==0, model.fdVSsed = nan; end
model.fdVSmoh = 100*diff(vs(zz==zmoh))./mean(vs(zz==zmoh));
model.Sanis = 100*(xi-1); % in percentage 
model.Panis = model.Sanis * par.mod.misc.psi_to_xi; % brb20240705 xi psi will need to be related through some constant parameter. 

% re-order fields of model structure 
forder ={'z','z0','VS','VP','rho','Nz','zsed','zmoh','vpvs','cxi','mxi','selev',...
         'fdVSsed','fdVSmoh','Sanis','Panis',...
         'sedmparm','crustmparm','mantmparm','datahparm','M'};
     
model = orderfields(model,forder);

ifplot = false; 
if ifplot; %brb20240704 Testing new mantle anisotropy 
    figure(2038); clf; hold on; 
    plot(model.Sanis, model.z, 'DisplayName','S anis');
    plot(model.Panis, model.z, 'DisplayName','P anis'); 
    set(gca, 'YDir', 'reverse'); 
    legend(); 
end

