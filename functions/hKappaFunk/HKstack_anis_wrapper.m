function [HK_A, HK_H, HK_K, t_pred] = HKstack_anis_wrapper(par, model, ...
                                    RF, tt, rayp, options)  
    arguments
        par % Eventually contain phase weights
        model % Model structure
        RF,
        tt,
        rayp, 
        options.ifplot = false; 
        options.phase_wts = [.5, .3, .2]; % brb2022.02.08 1/3 for each gave h-kappa stacks that don't look like others... They had way to tight of constraint on moho depth (for synthetics). I'm making a compromise. 
        options.hBounds = [10, 70]; 
        options.kBounds = [1.5, 2.1]; 
        options.kNum = 201; % kNum ~= hNum to help catch errors. 
        options.hNum = 200; 
        options.posterior = []
    end

% MCMC/THBI wrapper to get relevant model/data and use it for radially 
% anisotropic HK stack. 

% % % % % % Earth flattening transformation
a = 6371; 
zflt = @(z)-a .* log((a-z)./a); 
vflt = @(z,v)a./(a-z) .* v; 
% r = a - model.z; 
% zzf = -a .* log(r./a); 
% vsf = a./r .* vs; 
% vpf = a./r .* vp; 
% zz = zzf; 
% vs = vsf; 
% vp = vpf; 
model.z_orig = model.z; 
model.zmoh_orig = model.zmoh; 
model.z    = zflt(model.z   ); 
model.zmoh = zflt(model.zmoh); 
model.VS = vflt(model.z_orig, model.VS); 
model.VP = vflt(model.z_orig, model.VP); 

% Choose weighting for different phases. There are a few valid choices.
phase_wts = options.phase_wts; 

% Get average crust velocity. There are again a few choices. 
% Get vsv, vpv. If there is anisotropy, this influences HK calculation.
isCrust = model.z <= model.zmoh; 
isCrustF = find(isCrust); 
if model.z(isCrustF(end-1)) == model.z(isCrustF(end)); % Handle the duplicated moho depth for crust and mantle
    isCrustF(end) = []; % Remove mantle value
    isCrust(:) = false; 
    isCrust(isCrustF) = true; 
end
dz = zeros(size(model.z)); % For weighted mean, integral style. 
dz(1:end-1,1) = dz(1:end-1,1) + 0.5 .* diff(model.z);
dz(2:end  ,1) = dz(2:end  ,1) + 0.5 .* diff(model.z); 
dzC = dz(isCrust); 

% Function to calculate mean accounting for variable dz spacing. Simple integral mean. 
meanInt = @(m,dz)sum(dz.*m)/sum(dz); 
% vrAv = @(m,dz)( (1./meanInt(1./m,dz) + meanInt(m,dz) ) / 2 ); % Voigt Reuss Hill average (or whatever it's called). 
vrAv = @(m,dz)( (1./meanInt(1./m,dz) ) ); % This averaging approach is slightly more consistent with propmat than the voigt reuss. Makes sense... it's all about travel times. 

rho = meanInt(model.rho(isCrust), dzC); 
eta = 1; % IMPORTANT if eta is ever not 1, we need to modify the code. 
vsAv = vrAv(model.VS(isCrust), dzC); 
vpAv = vrAv(model.VP(isCrust), dzC); 

 % TODO need to incocporate sediment vpvs and xi into the averages properly. Tough to decide how to balance them being model parameters that we plot and look at, while HK sampling should be based on different parameters.  brb 2022.08.31, I cant really tell the difference between hk stacks or the positions of plotted points when changing between these two values. 
vpvs_av_iso = vpAv / vsAv; % including in sediment. 

if isempty(options.posterior); % Have to calculate averages from 1-D profile. 
    xi = model.Sanis(isCrust)/100+1; % Radial S anisotropy. add one because for some reason here anisotropy is + or -
    xi = vrAv(xi, dzC); 
    phi = model.Panis(isCrust)/100+1; % Radial P anisotropy
    phi = vrAv(phi, dzC); 
else % If we have the posterior, get average values from there. 
    posterior = options.posterior; 
    xi = median(posterior.xicrust); 
    phi = median(posterior.phicrust); 
end
    

% % %%brb2022.11.16 Old version of HK xi correction, where I only solved Christoffel equations for the vs and vp value in our model. Quicker. 
% % % % Assume ray path change is insignificant when adding anisotropy
% % % incAngVs = rayp2inc(rayp/111.1949, vsAv              ); 
% % % incAngVp = rayp2inc(rayp/111.1949, vsAv * vpvs_av_iso);
% % % 
% % % [velVs,~] = christof_radial_anis(vsAv, vsAv * vpvs_av_iso,...
% % %     xi, phi, eta, rho, incAngVs); % Velocities for ray with vs incidence angle
% % % [velVp,~] = christof_radial_anis(vsAv, vsAv * vpvs_av_iso,...
% % %     xi, phi, eta, rho, incAngVp); % Velocities for ray with vp incidence angle
% % % vsvAn = velVs(2); % Vsv for s ray 
% % % vpAn  = velVp(1); % Vp for p ray
% % % 
% % % % fprintf('\nvsvAn = %1.5f : vpAn = %1.5f\n', vsvAn, vpAn)
% % % 
% % % vpvsEffective = vpAn/vsvAn; 
% % % 
% % % [HK_A, HK_H, HK_K, t_pred] = HKstack(RF, tt, ...
% % %     rayp/111.1949, phase_wts, ...
% % %     vsAv, ... % isotropic
% % %     linspace(options.hBounds(1), options.hBounds(2), options.hNum)',...
% % %     linspace(options.kBounds(1), options.kBounds(2), options.kNum), ...
% % %     'vpvsANIS_vpvs',vpvsEffective/vpvs_av_iso,'Vsv_av',vsvAn); 
% % %%brb2022.11.16 END 
[HK_A, HK_H, HK_K, t_pred] = hk_anis(...
    RF, tt, rayp, vsAv, rho, xi, phi, eta); %brb2022.11.16 Solves Christoffel for each k, because each k values produces a different vp. Shouldn't matter for MCMC because we only sample HK stack at the h and k value in our model anyway! 
fhand=@()hk_anis(...
    RF, tt, rayp, vsAv, rho, xi, phi, eta); 
timeit(fhand)

% The rest of the script expects transposed hk stack. For IRIS ears data loaded from Zach's functions, K is first dimension, H is second bb2021.12.31
HK_A = HK_A'; 
HK_H = HK_H'; 
HK_K = HK_K'; 

% options.ifplot = true; if options.ifplot; warning('Setting options.ifplot = true'); end; 
if options.ifplot; 
    figure(198); clf; hold on; set(gcf,'color','white');
    subplot(1,1,1); hold on; 
    xlabel('kappa'); ylabel('H'); title('Synthetic H-kappa stack (black dot ignores sediment kappa)'); 
    set(gca, 'ydir', 'reverse');         
    sf = pcolor(HK_K, HK_H, HK_A'); %Can't use surf. It's 3d. It always covers other plotting objects. 
    sf.EdgeAlpha = 0; 
    colorbar(); 
    xlim([min(HK_K), max(HK_K)]); 
    ylim([min(HK_H), max(HK_H)]); 


    % An unfinished attempt to interpolate to find the exact point of maximum h-k energy. 
%         F = griddedInterpolant({HK_K, HK_H},HK_A,'spline');
%         HK_K2 = linspace(min(HK_K)+.1, max(HK_K)-.1, 1000); 
%         HK_H2 = linspace(min(HK_H)+.1, max(HK_H)-.1, 1001); 
%         [HK_K2, HK_H2] = ndgrid(HK_K2, HK_H2); 
%         HK_A2 = F(HK_K2, HK_H2); 
%         Emax = max(max(HK_A2));
%         [xmax,ymax] = find(Emax==HK_A2); 
%         scatter(gca, HK_K2(xmax), HK_H2(ymax), 50, 'red')

    % Simple wasy to estimate maximum h-k energy point. 
    Emax = max(max(HK_A));
    [xmax,ymax] = find(Emax==HK_A); 
    kBest = HK_K(xmax); 
    hBest = HK_H(ymax); 

    % Plot position of max energy. 
    scatter(gca, kBest, hBest, 50, 'red')
    text(kBest, hBest, sprintf('H = %2.2f, k = %1.3f',...
        HK_H(ymax),HK_K(xmax)), 'verticalalignment', 'bottom' )

    % Plot true position of max energy. TODO add to legend. 
    try
        scatter(model.vpvs, model.zmoh, 50, 'k'); 
    catch
        warning('Could not plot model.vpvs and/or model.zmoh'); 
    end

    % Plot the receiver function and the predicted times of different phases. This is a function in HKstack.m 
    HKstack(RF, tt, ...
        rayp/111.1949, phase_wts, vsAv, ... 
        linspace(10, 70, 2000)',linspace(1.5, 2.1, 2001), ...
            'plotPicks',true,'kChoice',kBest,'hChoice',hBest, ...
            'vpvsANIS_vpvs',vpvsEffective/model.vpvs,'Vsv_av',vsvAn); 

end

end