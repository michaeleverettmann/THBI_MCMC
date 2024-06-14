function [HK_A] = HK_stack_anis_wrapper_no_grid(par,model,RF,tt,rayp,options)

    arguments
        par % Eventually contain phase weights
        model % Model structure
        RF,
        tt,
        rayp, 
        options.phase_wts = [.5, .3, .2]; % brb2022.02.08 1/3 for each gave h-kappa stacks that don't look like others... They had way to tight of constraint on moho depth (for synthetics). I'm making a compromise. 
        options.posterior = [] % This shouldn't be needed
        options.plotVelAng = false; 
    end 

%%% Heavily based on HKstack_anis_wrapper.m

phase_wts = options.phase_wts; 

% Get vsv, vpv. If there is anisotropy, this influences HK calculation.
isCrust = model.z < model.zmoh; 
rho = mean(model.rho(isCrust)); % Medium density of crust. Not sure how best to average this yet.
eta = 1; % IMPORTANT if eta is ever not 1, we need to modify the code. 
vsAv = (1./mean(1./model.VS(isCrust)) + mean(model.VS(isCrust)) ) / 2; % mean VS. Decide which type of mean to take later. 
if isempty(options.posterior); % Have to calculate averages from 1-D profile. 
    xi = model.Sanis(isCrust)/100+1; % Radial S anisotropy. add one because for some reason here anisotropy is + or -
    xi = (1/mean(1./xi) + mean(xi)) / 2; % Voigt-Reuss-Hill average
    phi = model.Panis(isCrust)/100+1; % Radial P anisotropy
    phi = (1/mean(1./phi) + mean(phi)) / 2; 
else % If we have the posterior, get average values from there. 
    posterior = options.posterior; 
    xi = median(posterior.xicrust); 
    phi = median(posterior.phicrust); 
end

% Assume ray path change is insignificant when adding anisotropy
incAngVs = rayp2inc(rayp/111.1949, vsAv             ); 
incAngVp = rayp2inc(rayp/111.1949, vsAv * model.vpvs); % PN 20.6 was the incidence angle here. Seems reasonable. Same as for the real rays for US.CEH

% Christof calculations. There a loops in this function for the different rays. 
[velVs,~] = christof_radial_anis(vsAv, vsAv * model.vpvs,...
    xi, phi, eta, rho, incAngVs); % Velocities for ray with vs incidence angle
[velVp,~] = christof_radial_anis(vsAv, vsAv * model.vpvs,...
    xi, phi, eta, rho, incAngVp); % Velocities for ray with vp incidence angle
vsvAn = velVs(2,:); % Vsv for s ray 
vpAn  = velVp(1,:); % Vp for p ray

vpvsEffective = vpAn./vsvAn; 

% HKstack_no_grid is mostly vectorized. 
[HK_A] = HKstack_no_grid(RF, tt, ...
    rayp'/111.1949, phase_wts, ...
    vsAv, model.zmoh, model.vpvs,...
    'vpvsANIS_vpvs',vpvsEffective./model.vpvs,'Vsv_av',vsvAn); 

% Plot just to make sure above is working
if options.plotVelAng
    figure(2001); clf; hold on; 
    set(gcf, 'color', 'white', 'pos', [339 481 798 255]);
    sgtitle('Anisotropy velocities for P and Ps phases'); 
    subplot(1,2,1); hold on; 
    scatter(incAngVs, velVs(3,:), 'DisplayName', 'Vsh (s ray)'); 
    scatter(incAngVs, velVs(2,:), 'DisplayName', 'Vsv (s ray)'); 
    box on; grid on; xlabel('Incidence angle'); ylabel('Velocity'); 
    leg = legend(); leg.Location = 'best'; 
    subplot(1,2,2); hold on; box on; grid on; 
    scatter(incAngVp, velVp(1,:), 'DisplayName', 'Vp (p ray)'); 
    box on; grid on; xlabel('Incidence angle'); ylabel('Velocity'); 
    leg = legend(); leg.Location = 'best'; 
end

end