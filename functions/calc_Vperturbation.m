function [ modptb ] = calc_Vperturbation( model0,model1,ifplot )
%[ modptb ] = calc_Vperturbation( model0,model1,ifplot )
% 
% Function to calculate the fractional perturbation of model1 from the
% reference model0, as a function of depth, for vsv,vsh,vpv,vph,rho
% 
% perturbations are returned for all depths in corresponding card file
% (from surface to Earth centre)

if nargin <3 || isempty(ifplot)
    ifplot=false;
end

% ifplot = true; warning('if plot = true')
plot_upscale_interp = ifplot; 

%% brb2022.06.30. Interpolate models onto upsampled basis, then take m1 back to m0.z. 
% If a discontinuity depth changes, we would have problems. 
% Z indicies of m0 and m1 might not correspond to similar depths, and dv at some depth
% might not have the desired meaning anymore in our kernels K.
% Then using K to find dc stops working well. 
% Solution: Sample m1 onto m0.z, but do so by taking the most
% representative average of m1 in the viscinity of m0.z. 
% For a coarse m0.z, discontinuities can't be captured properly through
% simple interpolation. 
if ~ isequal(model0.z, model1.z); 
    if plot_upscale_interp; 
        model1_plot = model1; 
        model0_plot = model0; 
    end
   
    each_param = {'Panis', 'Sanis', 'rho', 'VP', 'VS', 'z0', 'z'}; % List all parameters in model that have nz values. 
    dz0 = diff(model0.z); 
    dz1 = diff(model1.z); 
    disc0 = find(dz0==0); 
    disc1 = find(dz1==0); 
    model0_unique = model0; 
    model1_unique = model1; 
    model0_unique.z(disc0  ) = model0_unique.z(disc0  ) - 0.0001; % Remove duplicate z values so we can use interp1. 
    model0_unique.z(disc0+1) = model0_unique.z(disc0+1) + 0.0001; 
    model1_unique.z(disc1  ) = model1_unique.z(disc1  ) - 0.0001; 
    model1_unique.z(disc1+1) = model1_unique.z(disc1+1) + 0.0001; 

    % Establish z points for upscaling, and spacing between points for weighting. 
    z_upsc = unique([...
        0; 0.01; 0.02; 0.03; 0.04; 0.05; 0.06; 0.08; 
        logspace(log10(0.1), log10(10),30)';...
        linspace(10, 100, 200)'; ...
        linspace(100, 300, 100)'; ...
        model0_unique.z]); % Will upscale/upsample to this basis, then back down to model0 basis. 
    mid_z = (z_upsc(1:end-1) + z_upsc(2:end) ) / 2; 
    dz_mid = zeros(size(z_upsc)); 
    dz_mid(2:end  ) =                   z_upsc(2:end  ) - mid_z; % Distance from z_upsc to midpoint above it. 
    dz_mid(1:end-1) = dz_mid(1:end-1) + mid_z - z_upsc(1:end-1); % Distance from z_upsc to midpoint below it. To check: does sum(dz_mid) = max(z_upsc)? Yes - brb2022.07.01. 
    old_z_midpts = [min(model0.z)-1;...
        (model0.z(2:end) + model0.z(1:end-1))./2;...
        max(model0.z)+1 ]; % For deciding how to associate z_upsc model values to which old z indecies. 

    % Upscale model1.  
    for i_param = 1:length(each_param); 
        param = each_param{i_param}; 
        model1_unique.(param) = interp1(...
            model1_unique.z, model1_unique.(param), z_upsc, 'linear'); % Linear avoids instabilities around discontinuities that we might find with cubic, splines, etc. brb2022.07.01
    end
    model1_unique.z = z_upsc; 

    nz0 = length(model0.z); 
    for i_param = 1:length(each_param); 
        param         = each_param{i_param}; 
        array_upsamp  = model1_unique.(param); % Upsampled resolution
        array_newbase = nan(size(model0.z)); % Lower resolution - main result of this block of code. 

        for idep = 1:nz0;
            this_dep = and( (old_z_midpts(idep  ) <= model1_unique.z) ,...
                            (old_z_midpts(idep+1) >=  model1_unique.z) ); % These upsampled incicies will correspond to the previous indecies. 
            
            if sum(this_dep) == 0; % If no upsampled depths in this range, just pick nearest depth. 
                [dist_best, i_dist] = min(abs(...
                    old_z_midpts(idep) - model1_unique.z)); 
                this_dep = i_dist; 
            end
                
            array_newbase(idep) = sum(...
                array_upsamp(this_dep) .* dz_mid(this_dep)) ...
                ./ sum(dz_mid(this_dep)); % Simple integral average 
        end
        if any(isnan(array_newbase)); % Should not occur very often
            fprintf('\nNan values during calc_Vperturbation. brb2022.07.29 There should already be code in place that fixes this. '); 
        end
        model1.(param) = array_newbase; 
    end    
    model1.z = model0.z; 

    if plot_upscale_interp; 
        figure(1005); clf; hold on; 
        h = tiledlayout(1,2,'TileSpacing', 'compact'); 
        nexttile(); cla; hold on; set(gca, 'ydir', 'reverse'); 
        box on; set(gca, 'LineWidth', 1.5); ylim([-5, 300]); 
        hmo  = plot(model0_plot.VS  , model0_plot.z); % Handle model 0. 
        hm1  = plot(model1_plot.VS  , model1_plot.z); 
        hm1i = plot(model1     .VS  , model1     .z); % Handle model 1 interpolated/averaged onto model 0 basis. 
        leg = legend([hmo hm1 hm1i], {'m0', 'm1 orig', 'm1 averaged onto m0 basis'}, 'location', 'best');    
        title('VS used for making Kbase'); 

        nexttile(); cla; hold on; set(gca, 'ydir', 'reverse'); 
        box on; set(gca, 'LineWidth', 1.5); ylim([-5, 60]); 
        hmo  = plot(model0_plot.VS       , model0_plot.z); % Handle model 0. 
        hm1  = plot(model1_plot.VS - 0.4 , model1_plot.z); 
        hm1i = plot(model1     .VS + 0.4 , model1     .z); % Handle model 1 interpolated/averaged onto model 0 basis. 
        title('Highlighting possible problems\nShifted for clarity'); 
    end
end

%% Radial anisotropy
xi0  = 1 + model0.Sanis/100;  % assumes Sanis is a percentage of anis about zero
phi0 = 1 + model0.Panis/100;  % assumes Panis is a percentage of anis about zero
[ vsv0,vsh0 ] = VsvVsh_from_VsXi( model0.VS,xi0 );
[ vpv0,vph0 ] = VpvVph_from_VpPhi( model0.VP,phi0 );

xi1  = 1 + model1.Sanis/100;  % assumes Sanis is a percentage of anis about zero
phi1 = 1 + model1.Panis/100;  % assumes Panis is a percentage of anis about zero
[ vsv1,vsh1 ] = VsvVsh_from_VsXi( model1.VS,xi1 );
[ vpv1,vph1 ] = VpvVph_from_VpPhi( model1.VP,phi1 );


%% get in card format inc. all depths in PREM
card0 = write_cardfile([],model0.z,vpv0,vsv0,model0.rho,[],[],vph0,vsh0);
card1 = write_cardfile([],model1.z,vpv1,vsv1,model1.rho,[],[],vph1,vsh1);

%% parse depths and nodes
zz0 = card0.depth;
zz1 = card1.depth;
N = length(zz0); % this many kernel rows

% find coincident nodes
[zzbotha,zbothi0a,zbothi1a] = intersect(zz0,zz1,'stable'); % gets first indices of discs
[zzbothb,zbothi0b,zbothi1b] = intersect(zz0,zz1,'legacy'); % gets second indices of discs
zbothi0 = union(zbothi0a,zbothi0b);
zbothi1 = union(zbothi1a,zbothi1b);

zdisc0 = card0.depth(diff(card0.depth)==0);
zdisc1 = card1.depth(diff(card1.depth)==0);

% account for 1's discontinuities matching single nodes in 0.
if length(zbothi0)~=length(zbothi1)
    d1 = setdiff(zdisc1,zdisc0);
    if any([d1])
        for ii = 1:length(d1)
            zbothi0 = sort([zbothi0;find(zz0==d1(ii))]);
        end
    end
end
if ~isequal(zz0(zbothi0),zz1(zbothi1)) % Because estimated dV for multiplying with surface waves kernels will not make sense if there are different depth vector lengths? I dont know... brb2022.05.19
    error('something wrong with indices'); 
end

%find non-coincident points
zdiffi0 = setdiff([1:N]',zbothi0);

%% voigt average velocities
[card0.vsvgt,card0.vpvgt] = voigtav(card0.vsh,card0.vsv,card0.vph,card0.vpv,card0.eta);
[card1.vsvgt,card1.vpvgt] = voigtav(card1.vsh,card1.vsv,card1.vph,card1.vpv,card1.eta);


%% calc dvals
% initialise dvals
dvsv = nan(N,1);
dvsh = nan(N,1);
dvpv = nan(N,1);
dvph = nan(N,1);
drho = nan(N,1);
dvsav = nan(N,1);
dvpav = nan(N,1);

% % insert dvals for coincident points
dvsv(zbothi0) = card1.vsv(zbothi1)./card0.vsv(zbothi0) - 1;
dvsh(zbothi0) = card1.vsh(zbothi1)./card0.vsh(zbothi0) - 1;
dvpv(zbothi0) = card1.vpv(zbothi1)./card0.vpv(zbothi0) - 1;
dvph(zbothi0) = card1.vph(zbothi1)./card0.vph(zbothi0) - 1;
drho(zbothi0) = card1.rho(zbothi1)./card0.rho(zbothi0) - 1;
dvsav(zbothi0) = card1.vsvgt(zbothi1)./card0.vsvgt(zbothi0) - 1;
dvpav(zbothi0) = card1.vpvgt(zbothi1)./card0.vpvgt(zbothi0) - 1;

if ~isempty(zdiffi0); 
    warning('m0.z ~= m1.z. This is odd: m1 should have been projected onto m0.z basis. brb2022.07.11'); 
    % insert dvals for non-coincident points (resolve onto card0 basis)
    dvsv(zdiffi0) = linterp(zz1,card1.vsv,zz0(zdiffi0))./card0.vsv(zdiffi0) - 1;
    dvsh(zdiffi0) = linterp(zz1,card1.vsh,zz0(zdiffi0))./card0.vsh(zdiffi0) - 1;
    dvpv(zdiffi0) = linterp(zz1,card1.vpv,zz0(zdiffi0))./card0.vpv(zdiffi0) - 1;
    dvph(zdiffi0) = linterp(zz1,card1.vph,zz0(zdiffi0))./card0.vph(zdiffi0) - 1;
    drho(zdiffi0) = linterp(zz1,card1.rho,zz0(zdiffi0))./card0.rho(zdiffi0) - 1;
    dvsav(zdiffi0) = linterp(zz1,card1.vsvgt,zz0(zdiffi0))./card0.vsvgt(zdiffi0) - 1;
    dvpav(zdiffi0) = linterp(zz1,card1.vpvgt,zz0(zdiffi0))./card0.vpvgt(zdiffi0) - 1;
end

% fix nans
dvsv(card0.vsv==0) = 0; if any(isnan(dvsv)), error('nans in dvsv'); end
dvsh(card0.vsh==0) = 0; if any(isnan(dvsh)), error('nans in dvsh'); end
dvpv(card0.vpv==0) = 0; if any(isnan(dvpv)), error('nans in dvpv'); end
dvph(card0.vph==0) = 0; if any(isnan(dvph)), error('nans in dvph'); end
drho(card0.rho==0) = 0; if any(isnan(drho)), error('nans in drho'); end
dvsav(card0.vsvgt==0) = 0; if any(isnan(dvsav)), error('nans in dvsvgt'); end
dvpav(card0.vpvgt==0) = 0; if any(isnan(dvpav)), error('nans in dvpvgt'); end


%% output
modptb =struct('Z',zz0,'dvsv',dvsv,'dvsh',dvsh,'dvpv',dvpv,'dvph',dvph,'drho',drho,'dvsav',dvsav,'dvpav',dvpav);

%% plot
if ifplot
    figure(44); clf; set(gcf,'pos',[331 384 1537 614])
    ax1 = subplot(1,5,1); 
    ax2 = subplot(1,5,2); 
    ax3 = subplot(1,5,3); 
    ax4 = subplot(1,5,4); 
    ax5 = subplot(1,5,5); 
    
    plot(ax1,dvsv,zz0,'linewidth',2);
    plot(ax2,dvsh,zz0,'linewidth',2);
    plot(ax3,dvpv,zz0,'linewidth',2);
    plot(ax4,dvph,zz0,'linewidth',2);
    plot(ax5,drho,zz0,'linewidth',2);
    
    set(ax1,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax2,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax3,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax4,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax5,'ydir','reverse','fontsize',15,'ylim',[0 500])
    
    title(ax1,'Vsv','fontsize',20);xlabel(ax1,'$\mathbf{dV_{SV}/V_{SV}}$','interpreter','latex','fontsize',19)
    title(ax2,'Vsh','fontsize',20);xlabel(ax2,'$\mathbf{dV_{SH}/V_{SH}}$','interpreter','latex','fontsize',19)
    title(ax3,'Vpv','fontsize',20);xlabel(ax3,'$\mathbf{dV_{PV}/V_{PV}}$','interpreter','latex','fontsize',19)
    title(ax4,'Vph','fontsize',20);xlabel(ax4,'$\mathbf{dV_{PH}/V_{PH}}$','interpreter','latex','fontsize',19)
    title(ax5,'rho','fontsize',20);xlabel(ax5,'$\mathbf{d\rho/\rho}$','interpreter','latex','fontsize',19)
    ylabel(ax1,'\textbf{Depth (km)}','interpreter','latex','fontsize',19)
    sgtitle('Comparisons from cardfile'); 

    figure(45); clf; set(gcf,'pos',[331 384+800 1537 614])
    ax1 = subplot(1,5,1); 
    ax2 = subplot(1,5,2); 
    ax3 = subplot(1,5,3); 
    ax4 = subplot(1,5,4); 
    ax5 = subplot(1,5,5); 
    
    plot(ax1,model0.VS,model0.z,'linewidth',2, 'color', 'k');
    plot(ax1,model1.VS,model1.z,'linewidth',2, 'color', 'r');
    plot(ax2,(model1.VS-model0.VS)./model0.VS,model0.z,'linewidth',2, 'color', 'k');
    plot(ax3,model0.VP,model0.z,'linewidth',2, 'color', 'k');
    plot(ax3,model1.VP,model1.z,'linewidth',2, 'color', 'r');
    plot(ax4,(model1.VP-model0.VP)./model0.VP,model0.z,'linewidth',2, 'color', 'k');
    plot(ax5,(model1.rho-model0.rho)./model0.rho,model0.z,'linewidth',2, 'color', 'k');

    
    set(ax1,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax2,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax3,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax4,'ydir','reverse','fontsize',15,'ylim',[0 500])
    set(ax5,'ydir','reverse','fontsize',15,'ylim',[0 500])
    
    title(ax1,'Vs','fontsize',20);%xlabel(ax1,'$\mathbf{dV_{SV}/V_{SV}}$','interpreter','latex','fontsize',19)
    title(ax2,'Dvs','fontsize',20);%xlabel(ax2,'$\mathbf{dV_{SH}/V_{SH}}$','interpreter','latex','fontsize',19)
    title(ax3,'Vp','fontsize',20);%xlabel(ax3,'$\mathbf{dV_{PV}/V_{PV}}$','interpreter','latex','fontsize',19)
    title(ax4,'Dvp','fontsize',20);%xlabel(ax4,'$\mathbf{dV_{PH}/V_{PH}}$','interpreter','latex','fontsize',19)
    title(ax5,'Drho','fontsize',20);%xlabel(ax5,'$\mathbf{d\rho/\rho}$','interpreter','latex','fontsize',19)
    ylabel(ax1,'\textbf{Depth (km)}','interpreter','latex','fontsize',19)
    sgtitle('(m1-m0)/m0: Comparison directly between m0 and m1 variables. '); 
    
end

end

