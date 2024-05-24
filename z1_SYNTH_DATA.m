function [ data ] = z1_SYNTH_DATA(par,ifplot)

% synthetic data parameters
samprate = par.synth.samprate;

for id = 1:length(par.inv.datatypes)
    pdtyps(id,:) = parse_dtype(par.inv.datatypes{id}); 
end

% [predata.PsRF.Vp_surf] = deal(mean([predata.PsRF.Vp_surf]));
% [predata.PsRF.Vs_surf] = deal(mean([predata.PsRF.Vs_surf]));


global TRUEmodel TLM

if isempty(par.synth.surf_Vp_Vs)
    Vp_surf = TLM.Vp(1);
    Vs_surf = TLM.Vs(1);
else
    Vp_surf = par.synth.surf_Vp_Vs(1);
    Vs_surf = par.synth.surf_Vp_Vs(2);
end

%% ===================  CALCULATE RECEIVER FUNCTIONS  ===================
if any(string(pdtyps(:,1))=='BW') || any(string(pdtyps(:,1))=='RF') || any(string(pdtyps(:,1))=='HKstack'); % bb2021.12.07 test to stop calculating bw data when it's not wanted. 

gcarcs = par.synth.gcarcs;
np = 0;
ns = 0;
for ig = 1:length(gcarcs) %loop over synth gcarcs - the way to get multiple synth data
gcarc = gcarcs(ig);
tpP = tauptime('deg',gcarc,'phases','P'); raypp = tpP.rayparameter;
tpS = tauptime('deg',gcarc,'phases','S'); rayps = tpS.rayparameter;
P_inc = rayp2inc(raypp,TLM.Vp(end),6371-TLM.zlayb(end));
S_inc = rayp2inc(rayps,TLM.Vs(end),6371-TLM.zlayb(end));

[trudat_ps,tt_ps] = run_propmat_or_telewavesim(par.synth.propmat_or_telewavesim, TLM,'TegPs','Ps',samprate, P_inc, raypp, par.synth.synthperiod); % in R,T,Z
[trudat_sp,tt_sp] = run_propmat_or_telewavesim(par.synth.propmat_or_telewavesim, TLM,'TegSp','Sp',samprate, S_inc, rayps, par.synth.synthperiod); % in R,T,Z

% channel order
trudat_ps_ZRT = trudat_ps(:,[3,1,2]); % in Z,R,T
trudat_sp_ZRT = trudat_sp(:,[3,1,2]); % in Z,R,T

if strcmp(par.forc.PSVorZR,'PSV')
    % convert to P-SV
    [trudat_ps_PSV(:,1),trudat_ps_PSV(:,2)] = ...
        Rotate_XZ_to_PSV(trudat_ps_ZRT(:,2),-trudat_ps_ZRT(:,1),Vp_surf,Vs_surf,rayp_sdeg2skm(raypp,TLM.zlayb(end)));
    [trudat_sp_PSV(:,1),trudat_sp_PSV(:,2)] = ...
        Rotate_XZ_to_PSV(trudat_sp_ZRT(:,2),-trudat_sp_ZRT(:,1),Vp_surf,Vs_surf,rayp_sdeg2skm(rayps,TLM.zlayb(end)));
elseif strcmp(par.forc.PSVorZR,'ZR')
    % keep as ZR (but kill T; Z positive UP)
    trudat_ps_PSV = trudat_ps_ZRT(:,[1,2]);         
    trudat_sp_PSV = trudat_sp_ZRT(:,[1,2]);         
end
        
% norm
trudat_ps_PSV = trudat_ps_PSV./maxab(trudat_ps_PSV(:,1)); % normalise on parental max, make positive
trudat_sp_PSV = trudat_sp_PSV./maxab(trudat_sp_PSV(:,2)); % normalise on parental max, make positive

if any(isnan(trudat_sp_PSV)) % this will happen if S wave is inhomogeneous; should be killed
    trudat_sp_PSV = zeros(size(trudat_sp_PSV));
    fprintf('Inhomogeneous S wave at %.1f distance!!\n',gcarc); 
end 

tt_Par = mean(tt_ps(trudat_ps_PSV(:,1)==max(trudat_ps_PSV(:,1))));% estimate main P-wave arrival time from first big upswing
tt_Sar = mean(tt_sp(trudat_sp_PSV(:,2)==max(trudat_sp_PSV(:,2)))); % estimate main S-wave arrival time from first big upswing

tt_ps = tt_ps - tt_Par;
tt_sp = tt_sp - tt_Sar;
tt_ps = round_level(tt_ps,0.001);
tt_sp = round_level(tt_sp,0.001);

Ps_widewind = [-10 50];
Sp_widewind = [-50 10];


inwind_ps = (tt_ps >= Ps_widewind(1)) & (tt_ps < Ps_widewind(2)); 
inwind_sp = (tt_sp >= Sp_widewind(1)) & (tt_sp < Sp_widewind(2)); 

% crop
trudat_ps_PSV = trudat_ps_PSV(inwind_ps,:);
trudat_sp_PSV = trudat_sp_PSV(inwind_sp,:);
tt_ps = tt_ps(inwind_ps);
tt_sp = tt_sp(inwind_sp);

% taper
trudat_ps_PSV = flat_hanning_win(tt_ps,trudat_ps_PSV,Ps_widewind(1),Ps_widewind(2),3); % 3s taper
trudat_sp_PSV = flat_hanning_win(tt_sp,trudat_sp_PSV,Sp_widewind(1),Sp_widewind(2),3); % 3s taper

% normalise to small energy, flip to positive main pulse
% normf_ps = trudat_ps_PSV(:,1)'*trudat_ps_PSV(:,1) + trudat_ps_PSV(:,2)'*trudat_ps_PSV(:,2);
% trudat_ps_PSV = trudat_ps_PSV/sqrt(normf_ps/2)./sign(maxab(trudat_ps_PSV(:)));

% normf_sp = trudat_sp_PSV(:,1)'*trudat_sp_PSV(:,1) + trudat_sp_PSV(:,2)'*trudat_sp_PSV(:,2);
% trudat_sp_PSV = trudat_sp_PSV/sqrt(normf_sp/2)./sign(maxab(trudat_sp_PSV(:)));

% add noise
if strcmp(par.synth.noisetype,'gaussian') || strcmp(par.synth.noisetype,'gauss')
    trudat_ps_PSV = trudat_ps_PSV + random('norm',0,par.synth.noise_sigma_BW_Ps,size(trudat_ps_PSV));
    trudat_sp_PSV = trudat_sp_PSV + random('norm',0,par.synth.noise_sigma_BW_Sp,size(trudat_sp_PSV)); 
end

if ifplot
	figure(58),clf,set(gcf,'pos',[2 275 1047 830])
    %Ps
    subplot(411),plot(tt_ps,trudat_ps_PSV(:,1),'k','linewidth',2), xlim(par.datprocess.Ps.Twin.def), ylabel('Z tru','fontsize',18), title('Ps TRUE','fontsize',22)
    subplot(412),plot(tt_ps,trudat_ps_PSV(:,2),'r','linewidth',2), xlim(par.datprocess.Ps.Twin.def), ylabel('R tru','fontsize',18)
    %Sp
    subplot(413),plot(tt_sp,trudat_sp_PSV(:,1),'k','linewidth',2), xlim(par.datprocess.Sp.Twin.def), ylabel('Z tru','fontsize',18), title('Sp TRUE','fontsize',22)
    subplot(414),plot(tt_sp,trudat_sp_PSV(:,2),'r','linewidth',2), xlim(par.datprocess.Sp.Twin.def), ylabel('R tru','fontsize',18)
    
    rfSave = struct('tt_ps', tt_ps, 'trudat_ps_PSV', trudat_ps_PSV, ...
                    'ps_twin', par.datprocess.Ps.Twin, ...
                    'tt_sp', tt_sp, 'trudat_sp_PSV', trudat_sp_PSV, ...
                    'sp_twin', par.datprocess.Sp.Twin); % Just saving this for making some figures.... Not necessary. brb2022.03.21
    save([par.res.resdir '/tru_receiver_functions.mat'],'rfSave'); 
        
end

% save
np = np+1;
PsRF(np,1) = struct('PSV',trudat_ps_PSV,'tt',tt_ps,'gcarc',gcarc,'rayp',raypp,'inc',P_inc,'samprate',samprate,'nsamp',length(tt_ps),'sigma',par.mod.data.prior_sigma.BW.Ps.def,'Vp_surf',Vp_surf,'Vs_surf',Vs_surf);
if gcarc>65 % don't do if too close - inhomogeneous
ns = ns+1;
SpRF(ns,1) = struct('PSV',trudat_sp_PSV,'tt',tt_sp,'gcarc',gcarc,'rayp',rayps,'inc',S_inc,'samprate',samprate,'nsamp',length(tt_sp),'sigma',par.mod.data.prior_sigma.BW.Sp.def,'Vp_surf',Vp_surf,'Vs_surf',Vs_surf);
end
% clear trudat_ps_PSV trudat_ps_ZRT 
% clear trudat_sp_PSV trudat_ps_ZRT

end % loop over gcarcs

%% ===================  MAKE DATA STRUCTURE  ===================
% data = struct('BW_Ps',PsRF,'BW_Sp',SpRF);
data = struct(); 
if any(string(par.inv.datatypes)=='BW_Ps') || any(string(par.inv.datatypes)=='RF_Ps') || any(string(par.inv.datatypes)=='HKstack_P')
    data.BW_Ps = PsRF; 
end
if any(string(par.inv.datatypes)=='BW_Sp') || any(string(par.inv.datatypes).contains('RF_Sp'))
    data.BW_Sp = SpRF; 
end

% %% ==================== Get H-kappa stack using P wave receiver function ==
% if any(string(pdtyps(:,1))=='HKstack'); % Use the Ps receiver function to get h-kappa stack
%     phase_wts = [.7, .2, .1]; % Zhu and Kanamori 2000 weights. Ordered as w1, w2, w3, but I'm not sure I Zach used the same ordering in HKstack.m bb2021.12.08
% %     [HK_A, HK_H, HK_K] = HKstack(PsRF.PSV(:,2), PsRF.tt, PsRF.rayp/111, phase_wts, ...
% %         PsRF.Vs_surf, linspace(25, 70, 200)',linspace(1.6, 2.1, 200)); 
%     [HK_A, HK_H, HK_K] = HKstack(PsRF.PSV(:,2), PsRF.tt, PsRF.rayp/111, phase_wts, ...
%         PsRF.Vs_surf, linspace(10, 70, 200)',linspace(1.5, 2.1, 201)); 
% %     warning('10 to 60 km h kappa')
%     
%     figure(1); clf; hold on; 
%     subplot(1,1,1); 
%     thing = surf(HK_K, HK_H, HK_A, 'EdgeAlpha', 0)
%     xlabel('kappa'); ylabel('H'); title('Synthetic H-kappa stack'); 
%     set(gca, 'ydir', 'reverse'); 
%     colorbar(); 
%     
%     % TODO what is Vs_surf? Average vs? Why named surf? 
%     % TODO change H vector. Is 25 to 70 in IRIS ears It seems... Might just
%     % load a real HK stack and use those same values. 
%     % PsRF.inc is P_inc. 
% %     S_inc0 = rayp2inc(rayps,TLM.Vs(1),6371-TLM.zlayb(1))/ 
%     
% %     HKstack_P = struct('H', HK_H, 'K', HK_K', 'Esum', HK_A, ...
% %         'Nobs', length(gcarcs), 'dof', length(gcarcs) ); 
%     HKstack_P = struct('H', HK_H, 'K', HK_K', 'Esum', HK_A, ...
%         'Nobs', 7, ...
%         'dof' , 7 ); 
%     warning('bb2021.12.20 artifically setting 7 h-kappa observations')
%     % Not sure what dof is. Not sure how to establish Nobs, or if it even
%     % matters. 
%     % For IRIS ears data loaded from Zach's functions, K is first
%     % dimension, H is second
%     data.HKstack_P = HKstack_P; 
% end

end % doing body waves

% % % %% Get Sp CCP stacks
% % % dz = 0.5; 
% % % zz = [-15:dz:300]';
% % % gcarcs = par.synth.gcarcs(1); 
% % % if length(par.synth.gcarcs)>1; error('Only have Sp CCP ready for 1 gcarc. You need to build a loop over more, or just use one gcarc.'); end
% % % tpS = tauptime('deg',gcarc,'phases','S'); 
% % % rayps = tpS.rayparameter;
% % % dof_per_z=1/par.datprocess.CCP.parent_zw; 
% % % data.RF_Sp_ccp = struct('zz',zz,'rayp',rayps,'dz',dz,'nsamp',length(zz),'dof_per_z',dof_per_z); 
% % % warning('Fake degree of freedom'); 
% % % 
% % % % par.inv.datatypes{6} = 'RF_Sp_ccp'
% % % [data] = b3_FORWARD_MODEL_RF_ccp(TRUEmodel, TLM, par, data, 'temporary', 1)
% % % %!%! Probably add noise here. 
% % % %%


%% ===================  CALCULATE PHASE VELOCITIES  ===================
% % Use Menke phV solver method on layered model 
% phinmod = [TLM.zlayb-TLM.zlayt,TLM.Vp,TLM.Vs,TLM.rho];
% trudat_cph = Calc_Ray_dispersion(SWperiods,phinmod,1,2000,0);
% % Use MINEOS
if ifplot, figure(13); clf; hold on; end
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id};
    pdtyp = parse_dtype(dtype);
    if ~strcmp(pdtyp{1},'SW'), continue, end
    nstr = ['noise_sigma_',pdtyp{1},'_',pdtyp{2}];
    SWperiods = par.synth.([dtype,'_periods']);
    if strcmp(pdtyp{2},'HV')
        truSWdat.HVr = run_HVkernel(TRUEmodel,SWperiods,'initmod',1,0,par.inv.verbose);
    else
        par_mineos = struct('R_or_L',pdtyp{2},'phV_or_grV',pdtyp{3},'ID','synthmod');
        [truSWdat.phV,truSWdat.grV] = run_mineos(TRUEmodel,SWperiods,par_mineos,1,0,par.inv.verbose); %MINEOS_REPLACE
    end

    % add noise
    truSWdat.(pdtyp{3}) = truSWdat.(pdtyp{3}) + random('norm',0,par.synth.(nstr),size(truSWdat.(pdtyp{3})));

    if ifplot
        hold on; plot(SWperiods,truSWdat.phV,'o')
    end

    data.(dtype) = struct('periods',SWperiods,pdtyp{3},truSWdat.(pdtyp{3}),'sigma',par.mod.data.prior_sigma.(pdtyp{1}).(pdtyp{2}).(pdtyp{3}));

    data.(dtype).for_mod_info = struct('n_periods_calc', length(SWperiods),...
        'all_periods', SWperiods,...
        'min_period',min(SWperiods),'max_period',max(SWperiods),...
        'periods_calc',SWperiods); % Information on what forward modelling needs to be done. We don't want to do forward modelling once for each seperate author. Instead, forward model once, accounting for every author simultaneously. 

    
end




end

