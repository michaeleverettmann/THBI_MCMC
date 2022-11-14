function [model,laymodel,par] = z0_SYNTH_MODEL_splinemoduse(par,ifplot,options)
    arguments
        par
        ifplot = false 
        options.xi_crust = 1.05; 
        options.h_crust = 45; 
        options.vpvs_crust = 1.75; 
    end

xi_crust = options.xi_crust; 



global TRUEmodel TLM

% ver = 'orig'; % Index for model parameters to choose from. 
versta = [par.stadeets.sta]'; % Treat the station as determining what model to load. 

if strcmp(versta, 'teststa'); versta = 'orig'; end; % Backward compatibility. 

use_splines = false; 

%% CHOOSE CUSTOM KEY PARAMETERS
% For important parameters I modified, I'll try to do a %<
% Or maybe put the parameters I want to modify from base at the top
if strcmp(versta, 'orig'); % Original(ish) version. 
    selev = 0;
    h_sed = 0;
    h_crust = 45;

    vs_sed = [3.3 3.3];

    kvs_crust = [3.3 3.33 3.74 3.81]';
    cknots = linspace(h_sed, h_sed+h_crust,3)';
    fcknots = (cknots-h_sed)/h_crust;
    k_crust = length(cknots);

    kvs_mantle = [4.1 4.3  4.4 4.0 4.1 4.43]';
    mknots = [h_sed+h_crust, 60, 110, 140, par.mod.maxz]';
    fmknots = (mknots-h_sed-h_crust)/(par.mod.maxz-h_sed-h_crust);
    k_mantle = length(mknots);

%     vpvs_crust = 1.8;
    vpvs_crust = options.vpvs_crust; 
%     xi_crust = 1.05;
    xi_mantle = 1.0; 
    
    use_splines = true; 
elseif strcmp(versta, 'sed_deep'); % Original(ish) version. 
    selev = 0; 
    h_sed = 4.1; %<
    h_crust = 45;

    vs_sed = [1.5 2.5]; %<

    kvs_crust = [3.3 3.33 3.74 3.81]';
    cknots = linspace(h_sed, h_sed+h_crust,3)';
    fcknots = (cknots-h_sed)/h_crust;
    k_crust = length(cknots);

    kvs_mantle = [4.1 4.3  4.4 4.0 4.1 4.43]';
    mknots = [h_sed+h_crust, 60, 110, 140, par.mod.maxz]';
    fmknots = (mknots-h_sed-h_crust)/(par.mod.maxz-h_sed-h_crust);
    k_mantle = length(mknots);

%     vpvs_crust = 1.8;
    vpvs_crust = options.vpvs_crust; 

%     xi_crust = 1.05;
    xi_mantle = 1.0; 
    
    use_splines = true; 
elseif strcmp(versta, 'crat_2mld'); % Original(ish) version. 
    selev = 0; 
    h_sed = 1; 
    h_crust = 45;

    vs_sed = [1.5 2.5]; 

    kvs_crust = [3.3 3.33 3.74 3.81]';
    cknots = linspace(h_sed, h_sed+h_crust,3)';
    fcknots = (cknots-h_sed)/h_crust;
    k_crust = length(cknots);

    kvs_mantle = [          4.1    4.3,  3.8,  4.4,  4.55, 4.65, 4.3, 4.8, 4.8,  4.825  4.7,  4.6, 4.5, 4.5  ]'; %<
    mknots =     [h_sed+h_crust,      60,   80,   100,  120,  140,  160, 180,  200,   220, 240, 260,    par.mod.maxz   ]'; %<
    
    
    fmknots = (mknots-h_sed-h_crust)/(par.mod.maxz-h_sed-h_crust);
    k_mantle = length(mknots);

%     vpvs_crust = 1.8;
    vpvs_crust = options.vpvs_crust; 


%     xi_crust = 1.05;
    xi_mantle = 1.0; 
    
    use_splines = true; 
end



if use_splines; 
       
    %% DERIVATIVE PARMS
    % DEPTHS
    cminz = h_sed;
    cmaxz = h_sed+h_crust;
    zc = unique([cminz:par.mod.dz:cmaxz,cmaxz])';
    mminz = h_sed+h_crust;
    mmaxz = par.mod.maxz + selev;
    zm = unique([mminz:par.mod.dz:mmaxz,mmaxz])';
    % CRUST splines
    % dzsp = (cmaxz-cminz)/(k_crust-2);
    % knots = [repmat(cminz,1,3),cminz:dzsp:cmaxz,repmat(cmaxz,1,3)]';
    % sp = fastBSpline.lsqspline(knots,2,linspace(cminz,cmaxz,k_crust)',kvs_crust); % dummy velocities as placeholder
    % cspbasis = sp.getBasis(zc); cspbasis = cspbasis(:,2:end-1);
    [ cspbasis ] = make_splines( cknots,par, zc);
    % MANTLE splines
    % dzsp = (mmaxz-mminz)/(k_mantle-2);
    % knots = [repmat(mminz,1,3),mminz:dzsp:mmaxz,repmat(mmaxz,1,3)]';
    % sp = fastBSpline.lsqspline(knots,2,linspace(mminz,mmaxz,k_mantle)',kvs_mantle); % dummy velocities as placeholder
    % mspbasis = sp.getBasis(zm); mspbasis = mspbasis(:,2:end-1);
    [ mspbasis ] = make_splines( mknots,par, zm);
    % OVERALL
    M = 1 + 2 + 1 + k_crust + k_mantle + 1;

    %% MAKE ALL PARAMETER STRUCTURES
    sed = struct('h',h_sed,'VS',vs_sed);
    crust = struct('h',h_crust,'Nsp',k_crust+1,'Nkn',k_crust,'VS_sp',kvs_crust,'vpvs',vpvs_crust,'xi',xi_crust,'splines',cspbasis,'knots',cknots,'fknots',fcknots,'z_sp',zc);
    mantle = struct('Nkn',k_mantle,'Nsp',k_mantle+1,'VS_sp',kvs_mantle,'xi',xi_mantle,'splines',mspbasis,'knots',mknots,'fknots',fmknots,'z_sp',zm);
    % data = struct('sigmaPsRF',par.mod.data.prior_sigmaPsRF,...
    %               'sigmaSpRF',par.mod.data.prior_sigmaSpRF,...
    %               'sigmaSW',par.mod.data.prior_sigmaSW);
    data = ([]);

    %% MODEL WITH ALL PARMS
    model = struct('sedmparm',sed,'crustmparm',crust,'mantmparm',mantle,...
                   'datahparm',data,'M',M,'selev',selev);

    %% TURN PARMS INTO REAL TARGET MODEL
    TRUEmodel = make_mod_from_parms(model,par);
    % % % plot_PARAMETERISATION( TRUEmodel )
    % % % plot_quickmodel(par,TRUEmodel,TRUEmodel); warning('Remove this line'); 
    % same format...
end


if ~use_splines; 
    % Use station name to determine which velocity values go in. 
    % First word in station name is which base velocity model to use from SEM
    % Remaining words go in at any order. Examples: sed0 means 0 depth sediment.
    % See below for which key words will do what to your models. 
    
    % Test: ver = 'cont_EProt-sed4-mld1-mld2-zmoh25'
    
    % v_mod: tells which SEM model to use. 
    % v_alt: alterations we will manually do. v_alt explanation: 
    % s: sediment. 
    % m: mid-lithospheric velocity discontinuity. 
    % h: Moho depth. 
    versplit = split(versta, '-'); 
    v_mod = versplit{1}; 
    if length(versplit) == 2;
        v_alt = versplit{2}; 
    elseif length(versplit) == 1; 
        v_alt = ''; 
    else
        error('You put in a strange format for synthetic station name'); 
    end; 

    
    % Base crust velocities. 
    % % % is_crust = find(TRUEmodel.z==45); 
    % % % is_crust = is_crust(1); 
    % % % vs_c = TRUEmodel.VS(1:is_crust); 
    % % % z_c = TRUEmodel.z(1:is_crust); 
    vs_c = [3.3000    3.3067    3.3162    3.3284    3.3435    3.3612    3.3818    3.4051    3.4312 3.4600    3.4916    3.5260    3.5617    3.5955    3.6272    3.6567    3.6840    3.7093 3.7324    3.7534    3.7722    3.7889    3.8035    3.8100]'; 
    z_c = [0     2     4     6     8    10    12    14    16    18    20    22    24    26    28    30    32    34 36    38    40    42    44    45]'; 
    
    % Increase velocity at depth a little
    vs_c = vs_c + 0.1 * z_c./max(z_c); 
    fprintf('Synth model: Increasing lower sed velocity a little (0.1 km/s)\n'); 
    
    % Load mantle model. 
    if strcmp(v_mod, 'layermantle') 
        z0 = [0:2:300]; 
        vs0 = ones(size(z0)); 
        vs0(z0<80) = 3.9; 
        vs0(z0>=80) = 4.2; 
        vs0(z0>=150) = 4.5; 
        fprintf('Synth model: Making simple layered mantle.\n'); 
    else 
        SEMum2_avg = SEMum2_avgprofiles(); 
        z0  = SEMum2_avg.Z'; 
        vs0 = SEMum2_avg.(v_mod);
        fprintf('Synth model: Getting basic mantle model from SEMum2.\n'); 
    end
    
    if contains(v_alt, 'h25'); 
        zmoh = 25; 
        killc = z_c >= zmoh; 
        vs_base = interp1(z_c, vs_c, zmoh); 
        vs_c = [vs_c(~killc); vs_base]; 
        z_c  = [z_c(~killc) ; zmoh   ];
        fprintf('Synth model: Making shallow moho at %1.0f km.\n', zmoh); 
    end
        
        

    % "continuous" model. 
    z_m = [0:par.mod.dz:300]'; 
    z_m = unique(sort([z_m; z_c(end)])); % Make sure Moho depth is in mantle model. 
    z_m = z_m(z_m>=z_c(end)); % Just keep mantle depths and vals from SEM model. 
    vs_m = interp1(z0, vs0, z_m, 'spline', 'extrap');
    
    
    dvs_moh = vs_m(1) - vs_c(end); 
    dvs_moh_max = 0.5; 
    if dvs_moh > dvs_moh_max;
        fprintf('Synth model: Lowering uppermost mantle velocity to prevent large zmoh vs jump. New v jump of %1.3f, fractional perturbation of %1.3f\n',dvs_moh_max, dvs_moh_max./vs_c(end));  
        dvs_shift = dvs_moh - dvs_moh_max; % How much we need to decrease vs_m at zmoh
        vs_m_mod = interp1([min(z_m), 100, 105, max(z_m)], ... % Make the vs change stop by 100 km depth more or less. 
            [dvs_shift,0, 0, 0],...
            z_m, 'makima'); % Basically a linear interpolation, but makima makes this a little smoth
        vs_m = vs_m - vs_m_mod; 
%         figure(1); clf; hold on; plot(vs_m_mod, z_m); set(gca, 'ydir', 'reverse')
    end
        

    % % % Quick plot to make sure interpolation and extrapolation aren't messing anything up too much. 
    % figure(1); clf; hold on; set(gca, 'ydir', 'reverse'); plot(vs, z); plot(vs0, z0); 
    
    if contains(v_alt, 'm1'); 
        fprintf('Synth model: Adding low velocity at 80 km.\n'); 
        ftr_amp   = -0.2; 
        ftr_pos   = 80; 
        ftr_width = 10; 
        ftr = ftr_amp * exp( -(z_m-ftr_pos).^2 / (2*ftr_width^2) ); % Feature to add
        vs_m = vs_m + ftr; 
    end
    
    if contains(v_alt, 'm2')
        fprintf('Synth model: Adding low velocity at 150 km.\n'); 
        ftr_amp   = -0.2; 
        ftr_pos   = 150; 
        ftr_width = 10; 
        ftr = ftr_amp * exp( -(z_m-ftr_pos).^2 / (2*ftr_width^2) ); % Feature to add
        vs_m = vs_m + ftr;     
    end
    
    

    if contains(v_alt, 's1'); 
        sed = struct('h',1,'VS',[2.0, .95*vs_c(1)]);
%         sed = struct('h',1,'VS',[1.3, 2.3]);
    elseif contains(v_alt, 's4'); 
        sed = struct('h',4,'VS',[2.0, .95*vs_c(1)]);
    elseif contains(v_alt, 's4sharp'); % See what happens if the sediment/basement vs jump is sharp. I think this produces strong changes to HV, the inversion cannot accurately recover the sharpness of the feature, so HV isn't fit well, and everything goes to shit. Testing... brb2022.08.30
        sed = struct('h',4,'VS',[2.0, 2.5]);
    elseif contains(v_alt, 's2'); 
        sed = struct('h',2,'VS',[2.0, .95*vs_c(1)]);
    elseif contains(v_alt, 's0'); 
        sed = struct('h',0,'VS',[3.3 3.3]);
    else
        sed = struct('h',0,'VS',[3.3 3.3]);
    end
    % Don't worry about making this a discontinuity. 
    % It's useful anyway to test a parameterisation different than what we use for forward modelling. 
    killc = z_c < sed.h; 
    vs_c  = vs_c(~killc); 
    z_c   = z_c (~killc); 

    vpvs_crust = options.vpvs_crust; 
    
    %% MAKE ALL PARAMETER STRUCTURES
    crust = struct('h',max(z_c)-sed.h,'vpvs',vpvs_crust,'xi',xi_crust); % warning('Trying to crustal anisotropy'); 
    mantle = struct('xi',1);
    model = struct('sedmparm',sed,'crustmparm',crust,'mantmparm',mantle,...
                   'M', nan, 'datahparm', nan, 'selev',0);

               
    %% TURN PARMS INTO REAL TARGET MODEL
    TRUEmodel = make_mod_from_parms(model,par,'use_splines',false,...
        'vs', struct('crust', vs_c, 'mantle', vs_m),...
        'z' , struct('crust', z_c, 'mantle', z_m));
    
    if ifplot 
        figure(1); clf; hold on; box on; grid on; set(gcf, 'color', 'white'); 
        xlabel('Vs (km/s'); ylabel('Depth (km)'); title('Synthetic model'); 
        plot(TRUEmodel.VS, TRUEmodel.z);
        set(gca, 'ydir', 'reverse');           
    end

    h_crust = model.sedmparm.h + model.crustmparm.h; 
end



TRUEmodel.Z = TRUEmodel.z;
TRUEmodel.vs = TRUEmodel.VS;
TRUEmodel.vp = TRUEmodel.VP;


% % RHO is effed up = we don't agree on rho scaling, but it matters a lot for
% % the elastic moduli that propmat actually uses. 
% rho = linterp(trymod.z,trymod.rho,Z);
% % vs = linterp(trymod.z,trymod.VS,Z)+0.3;
% vp = linterp(trymod.z,trymod.VP,Z);



%% ===================  LAYERISE PROFILE  ===================
[zlayt,zlayb,Vslay] = ...
    layerise(TRUEmodel.z,TRUEmodel.vs,par.forc.mindV/3,0); 
nlay = length(Vslay);

% S to P and rho structure
xs = 1:find(zlayb==TRUEmodel.zsed); if TRUEmodel.zsed ==0, xs = []; end
xc = find(zlayt==TRUEmodel.zsed):find(zlayb==TRUEmodel.zmoh);
xm = find(zlayt==TRUEmodel.zmoh):nlay;
Vplay = [sed_vs2vp(Vslay(xs));...
         TRUEmodel.vpvs*Vslay(xc);...
         mantle_vs2vp(Vslay(xm),mean([zlayt(xm),zlayb(xm)],2))];
rholay = [sed_vs2rho(Vslay([xs,xc]));...
          mantle_vs2rho(Vslay(xm),mean([zlayt(xm),zlayb(xm)],2))];
xilay = [ones(length(xs),1);...
         TRUEmodel.cxi*ones(length(xc),1);...
         TRUEmodel.mxi*ones(length(xm),1)]; % S radial anisotropy
philay = ones(nlay,1); % P radial anisotropy
% warning('Changing philay')
% philay2 = 1./xilay; 
% philay = (philay + philay2) ./2; 
% disp(philay);
% philay = philay2; 
warning('NO p wave anisotropy!!!')
etalay = ones(nlay,1); % eta anisotropy

TLM = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay,'xi',xilay,'phi',philay,'eta',etalay);
if any(isnan(TLM.rho))
    error('NaN densities')
end

%% ====== Calculate average VS in crust. For h-kappa stacks === brb2022.02.08

z2 = TRUEmodel.z; % Need to operate on z a little bit. 
indLay = z2(2:end) == z2(1:end-1); % Layers indecies. 
indLay = find(indLay)+1; % Bottom of layers. 
z2(indLay) = z2(indLay) + 0.01; % Need unique points. 
timeToDep = cumtrapz(z2, 1./TRUEmodel.vs); % Time it would take a vertically travelling ray to go from surface to some depth. 

TRUEmodel.vsAvCrust = h_crust / interp1(z2, timeToDep, h_crust); % "Average" travel time in crust. 
% brb2022.02.09: This assumes a station with zero elevation. It should be fine for most synthetic tests. TODO we do have an optional elevation value in the synthetic model. I should verify that z includes z<0 - if so, no code change is needed. 

%% PLOT FINAL MODEL
% warning('if plot set true'); ifplot = 1; 
if ifplot
figure(95); clf; set(gcf,'pos',[120 151 920 947])
subplot(131), hold on;
plot(TRUEmodel.vp,TRUEmodel.z,'-b','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);
subplot(132), hold on;
plot(TRUEmodel.vs,TRUEmodel.z,'-b','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);
subplot(133), hold on;
plot(TRUEmodel.rho,TRUEmodel.z,'-b','Linewidth',1.5);
set(gca,'ydir','reverse','fontsize',14);
end

model = TRUEmodel;
laymodel = TLM;


end

