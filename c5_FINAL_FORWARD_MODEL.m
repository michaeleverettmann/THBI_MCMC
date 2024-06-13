function [ final_predata ] = c5_FINAL_FORWARD_MODEL( final_model,par,data,posterior )
% Do forward model from final median model. This script should follow the
% same procedure as forward modelling during the inversion. 

ID = [par.data.stadeets.nwk '.' par.data.stadeets.sta '_finalmod'];

%% ===================  PUT INTO DATA STRUCTURE  ===================
final_predata = data;

for id = 1:length(par.inv.datatypes)
    pdtyps(id,:) = parse_dtype(par.inv.datatypes{id}); 
end

%% ===================  LAYERISE PROFILE  ===================

if isfield(par.datprocess.CCP, 'layerise_version') && strcmp(par.datprocess.CCP.layerise_version, 'no_sed'); % Remove sediment for propmat. 
    [zlayt,zlayb,Vslay,...
        Vplay,rholay,xilay,philay] = ...
        layerise_test_sediment_sp_ccp(final_model.Z,final_model.VSav,par.forc.mindV,0,...
        final_model.VPav, final_model.rhoav, ...
        final_model.Sanisav./100+1,final_model.Panisav./100+1); 
else
    [zlayt,zlayb,Vslay,...
        Vplay,rholay,xilay,philay] = ...
        layerise(final_model.Z,final_model.VSav,par.forc.mindV,0,...
        final_model.VPav, final_model.rhoav, ...
        final_model.Sanisav./100+1,final_model.Panisav./100+1); 
end
nlay = length(Vslay);
etalay = ones(nlay,1); % eta anisotropy TODO get eta from model. eta is not in model right now, so can't yet. 

laymodel = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay,'xi',xilay,'phi',philay,'eta',etalay);
if any(isnan(laymodel.rho))
    error('NaN densities')
end

%% ===================  PS RFs FROM PROPAGATOR MATRIX  ====================

if any(strcmp(pdtyps(:,2),'Ps'))
    Psdat = par.inv.datatypes{find(strcmp(pdtyps(:,2),'Ps'),1,'first')};
    [ unique_rayps_P,irayps_P ] = rayp_vals( [final_predata.(Psdat).rayp] );
    for ir = 1:length(unique_rayps_P)
        rayp = unique_rayps_P(ir);
        samprate = unique([final_predata.(Psdat)(irayps_P==ir).samprate]);
        P_inc = rayp2inc(rayp,laymodel.Vp(end),6371-laymodel.zlayb(end));
        [predat_ps,tt_ps] = run_propmat_or_telewavesim(par.synth.propmat_or_telewavesim, laymodel,ID,'Ps',samprate, P_inc, rayp, par.forc.synthperiod,par.forc.nsamps);
        % pad with zeros
        tt_ps = [tt_ps(1) + [-1000:-1]'./samprate; tt_ps ;tt_ps(end) + [1:1000]'./samprate];
        predat_ps = [zeros(1000,3);predat_ps;zeros(1000,3)];
        %correct corrdinate order
        predat_ps_ZRT = predat_ps(:,[3,1,2]); % in Z,R,T
        if strcmp(par.forc.PSVorZR,'PSV')
            clear predat_ps_PSV;
            % convert to P, SV
            [predat_ps_PSV(:,1),predat_ps_PSV(:,2)] = ...
                Rotate_XZ_to_PSV(predat_ps_ZRT(:,2),-predat_ps_ZRT(:,1),...
                mean([final_predata.(Psdat).Vp_surf]),mean([final_predata.(Psdat).Vs_surf]),...
                rayp_sdeg2skm(rayp,laymodel.zlayb(end)));
        elseif strcmp(par.forc.PSVorZR,'ZR')
            % keep as ZR (but kill T; Z positive UP)
            predat_ps_PSV = predat_ps_ZRT(:,[1,2]);         
        end
        predat_ps_PSV = predat_ps_PSV./maxab(predat_ps_PSV(:,1)); % normalise on parental max, make positive
        tt_ps_Par = mean(tt_ps(predat_ps_PSV(:,1)==max(predat_ps_PSV(:,1))));% estimate main P-wave arrival time from first big upswing
        tt_ps = tt_ps - tt_ps_Par;
        tt_ps = round_level(tt_ps,0.001);

        Ps_widewind = [-10 50];
        inwind = (tt_ps >= Ps_widewind(1)) & (tt_ps < Ps_widewind(2)); 
        % crop
        predat_ps_PSV = predat_ps_PSV(inwind,:);
        tt_ps = tt_ps(inwind);  
        % taper
        predat_ps_PSV = flat_hanning_win(tt_ps,predat_ps_PSV,Ps_widewind(1),Ps_widewind(2),3); % 3s taper
        % normalise to unit energy
        %normf_ps = predat_ps_PSV(:,1)'*predat_ps_PSV(:,1) + predat_ps_PSV(:,2)'*predat_ps_PSV(:,2);
        %predat_ps_PSV = predat_ps_PSV/sqrt(normf_ps);
        % -----------------  PUT INTO DATA STRUCTURE  -----------------
        inds = find(irayps_P==ir);
        for iir = 1:length(inds)
            final_predata.(Psdat)(inds(iir),1).PSV=predat_ps_PSV; 
            final_predata.(Psdat)(inds(iir),1).tt=tt_ps;
            final_predata.(Psdat)(inds(iir),1).nsamp = length(final_predata.(Psdat)(inds(iir)).PSV);
        end
    end
end

%% ===================  SP RFs FROM PROPAGATOR MATRIX  ====================
if any(strcmp(pdtyps(:,2),'Sp'))
    if any(~strcmp(pdtyps(strcmp(pdtyps(:,2),'Sp'),3),'ccp'))
        Spdat = par.inv.datatypes{find(strcmp(pdtyps(:,2),'Sp'),1,'first')};
        [ unique_rayps_S,irayps_S ] = rayp_vals( [final_predata.(Spdat).rayp] );
        for ir = 1:length(unique_rayps_S)
            rayp = unique_rayps_S(ir);
            samprate = unique([final_predata.(Spdat)(irayps_S==ir).samprate]);
            S_inc = rayp2inc(rayp,laymodel.Vs(end),6371-laymodel.zlayb(end));
            
            P_inc = rayp2inc(rayp,laymodel.Vp(end),6371-laymodel.zlayb(end));
            % find layers where S to P conversion will not go inhomogeneous
            Play_incs = asind(laymodel.Vp*sind(P_inc)./laymodel.Vp(end));
            if any(~isreal(Play_incs))
                nimagplay = [1:find(imag(Play_incs),1,'first')-1];
                fns = fieldnames(laymodel);
                laymodel_Suse = laymodel; 
                laymodel_Suse.nlay = length(nimagplay);
                for jj = 1:length(fns)
                    if length(laymodel.(fns{jj}))==1, continue; end
                    laymodel_Suse.(fns{jj}) = laymodel_Suse.(fns{jj})(nimagplay);
                end
                S_inc = rayp2inc(rayp,laymodel_Suse.Vs(end),6371-laymodel_Suse.zlayb(end));
            else
                laymodel_Suse = laymodel;
            end
            
            [predat_sp,tt_sp] = run_propmat_or_telewavesim(par.synth.propmat_or_telewavesim, laymodel_Suse,ID,'Sp',samprate, S_inc, rayp, par.forc.synthperiod,par.forc.nsamps);
            % pad with zeros
            tt_sp = [tt_sp(1) + [-1000:-1]'./samprate; tt_sp ;tt_sp(end) + [1:1000]'./samprate];
            predat_sp = [zeros(1000,3);predat_sp;zeros(1000,3)];
            %correct corrdinate order
            predat_sp_ZRT = predat_sp(:,[3,1,2]); % in Z,R,T
            if strcmp(par.forc.PSVorZR,'PSV')
                clear predat_sp_PSV;
                % convert to P, SV
                [predat_sp_PSV(:,1),predat_sp_PSV(:,2)] = ...
                    Rotate_XZ_to_PSV(predat_sp_ZRT(:,2),-predat_sp_ZRT(:,1),...
                    mean([final_predata.(Spdat).Vp_surf]),mean([final_predata.(Spdat).Vs_surf]),...
                    rayp_sdeg2skm(rayp,laymodel_Suse.zlayb(end)));
            elseif strcmp(par.forc.PSVorZR,'ZR')
                % keep as ZR (but kill T; Z positive UP)
                predat_sp_PSV = predat_sp_ZRT(:,[1,2]);         
            end
            
            predat_sp_PSV = predat_sp_PSV./maxab(predat_sp_PSV(:,2)); % normalise on parental max, make positive
            tt_sp_Sar = mean(tt_sp(predat_sp_PSV(:,2)==max(predat_sp_PSV(:,2)))); % estimate main S-wave arrival time from first big upswing
            tt_sp = tt_sp - tt_sp_Sar;
            tt_sp = round_level(tt_sp,0.001);
    
            Sp_widewind = [-50 10];
            inwind = (tt_sp >= Sp_widewind(1)) & (tt_sp < Sp_widewind(2)); 
            % crop
            predat_sp_PSV = predat_sp_PSV(inwind,:);
            tt_sp = tt_sp(inwind);
            % taper
            predat_sp_PSV = flat_hanning_win(tt_sp,predat_sp_PSV,Sp_widewind(1),Sp_widewind(2),3); % 3s taper
            % normalise to unit energy
            normf_sp = predat_sp_PSV(:,1)'*predat_sp_PSV(:,1) + predat_sp_PSV(:,2)'*predat_sp_PSV(:,2);
            predat_sp_PSV = predat_sp_PSV/sqrt(normf_sp);
    % -----------------  PUT INTO DATA STRUCTURE  -----------------
            inds = find(irayps_S==ir);
            for iir = 1:length(inds)
                final_predata.(Spdat)(inds(iir),1).PSV=predat_sp_PSV; 
                final_predata.(Spdat)(inds(iir),1).tt=tt_sp;
                final_predata.(Spdat)(inds(iir),1).nsamp = length(final_predata.(Spdat)(inds(iir)).PSV);
            end
        end
    end % Sp not ccp
end

%% distribute data for different processing (e.g. _lo, _cms)
% for idt = 1:length(par.inv.datatypes)
%     dtype = par.inv.datatypes{idt}; 
%     pdt = parse_dtype( dtype ); 
%     if strcmp(pdt{1},'BW') && (~strcmp(pdt{3},'def') || ~strcmp(pdt{4},'def'))
%         final_predata.(dtype) = final_predata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
%     end
% end

%% ===================  CALCULATE SURFACE WAVE VELOCITIES  ===================
if any(strcmp(pdtyps(:,1),'SW'))
    modminrun = struct('z',final_model.Z,...
                       'VS',final_model.VSav,...
                        'VP',final_model.VPav,...
                        'rho',final_model.rhoav,...
                        'Sanis',final_model.Sanisav,...
                        'Panis',final_model.Panisav);
                    
    if any(strcmp(pdtyps(:,2),'Ray')), itp = par.inv.datatypes(find(strcmp(pdtyps(:,2),'Ray'),1,'first'));
        par_mineos = struct('R_or_L','R','ID',ID);
        [SW.Ray.phV,SW.Ray.grV] = run_mineos(...
            modminrun,data.(itp{1}).for_mod_info.all_periods,par_mineos,1,0,par.inv.verbose); %MINEOS_REPLACE
    end
    if any(strcmp(pdtyps(:,2),'Lov')), itp = par.inv.datatypes(find(strcmp(pdtyps(:,2),'Lov'),1,'first'));
        par_mineos = struct('R_or_L','L','ID',ID);
        [SW.Lov.phV,SW.Lov.grV] = run_mineos(...
            modminrun,data.(itp{1}).for_mod_info.all_periods,par_mineos,1,0,par.inv.verbose); %MINEOS_REPLACE
    end
    if any(strcmp(pdtyps(:,2),'HV')), itp = par.inv.datatypes(find(strcmp(pdtyps(:,2),'HV'),1,'first'));
        SW.HV.HVr = run_HVkernel(...
            modminrun,data.(itp{1}).for_mod_info.all_periods,'final',1,0,par.inv.verbose);
    end
    for id = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{id}; pdtyp=parse_dtype(dtype); if ~strcmp(pdtyp{1},'SW'),continue; end
        final_predata.(dtype).(pdtyp{3}) = interp1(data.(dtype).for_mod_info.all_periods, ...
            SW.(pdtyp{2}).(pdtyp{3}), data.(dtype).periods);
    end
end

%% ===================  SP RFs FROM PROPAGATOR MATRIX  ====================
if any(strcmp(pdtyps(strcmp(pdtyps(:,2),'Sp'),3),'ccp'))
    rayp = final_predata.RF_Sp_ccp.rayp;
    % declare samprate of 10 - moot as migrated to depth and interpolated later
    samprate = 10;

    S_inc = rayp2inc(rayp,laymodel.Vs(end),6371-laymodel.zlayb(end));
    P_inc = rayp2inc(rayp,laymodel.Vp(end),6371-laymodel.zlayb(end));
    % check if S inhomogeneous
    if isreal(asind(laymodel.Vs*sind(S_inc)./laymodel.Vs(end))) % 

        % find layers where S to P conversion will not go inhomogeneous
        Play_incs = asind(laymodel.Vp*sind(P_inc)./laymodel.Vp(end));
        if any(~isreal(Play_incs))
            nimagplay = [1:find(imag(Play_incs),1,'first')-1];
            fns = fieldnames(laymodel);
            laymodel_Suse = laymodel; 
            laymodel_Suse.nlay = length(nimagplay);
            for jj = 1:length(fns)
                if length(laymodel.(fns{jj}))==1, continue; end
                laymodel_Suse.(fns{jj}) = laymodel_Suse.(fns{jj})(nimagplay);
            end
            % set appropriate S_inc for the actual base
            S_inc = rayp2inc(rayp,laymodel_Suse.Vs(end),6371-laymodel_Suse.zlayb(end));
        else
            laymodel_Suse = laymodel;
        end
        [predat_sp,tt_sp] = run_propmat_or_telewavesim(par.synth.propmat_or_telewavesim, laymodel_Suse,ID,'Sp',samprate, S_inc, rayp, par.forc.synthperiod,par.forc.nsamps);
        % pad with zeros
        tt_sp = [tt_sp(1) + [-1000:-1]'./samprate; tt_sp ;tt_sp(end) + [1:1000]'./samprate];
        predat_sp = [zeros(1000,3);predat_sp;zeros(1000,3)];
        %correct corrdinate order
        predat_sp_ZRT = predat_sp(:,[3,1,2]); % in Z,R,T

        % convert to P, SV - use LAYMODEL Vp,Vs to do this precisely
        [predat_sp_PSV(:,1),predat_sp_PSV(:,2)] = ...
            Rotate_XZ_to_PSV(predat_sp_ZRT(:,2),-predat_sp_ZRT(:,1),... % z positive down
            laymodel.Vp(1),laymodel.Vs(1),...
            rayp_sdeg2skm(rayp,laymodel_Suse.zlayb(end)));

        predat_sp_PSV = predat_sp_PSV./maxab(predat_sp_PSV(:,2)); % normalise on parental max, make positive
        tt_sp_Sar = mean(tt_sp(predat_sp_PSV(:,2)==max(predat_sp_PSV(:,2)))); % estimate main S-wave arrival time from first big upswing
        tt_sp = tt_sp - tt_sp_Sar;
        tt_sp = round_level(tt_sp,0.001);
        
        %%% Introduce a fake synthetic parent pulse. To eliminate sediment reverberations from the parent. 
        %%% The parent pulse used in propmat needs to be exactly the same as I'm defining here. 
        if ~isfield(par.datprocess.CCP, 'simple_parent_pulse'); 
            warning('simple_parent_pulse was not an option for these results. Using propmat SV as parent pulse.'); 
        elseif par.datprocess.CCP.simple_parent_pulse; 
            % Make waveform parent pulse that propmat makes in sourc1.gauss.f, subroutine incid1
            std_fact_in_propmat_sourc1_gauss = 1/6; % This factor of 6 is hardwired into propmat!!! If pulse widths start seeming wrong, look in sourc1.gauss (or maybe sourc1?) at the subroutine incid1. npstd=peri*npps/6, that's where 6 comes from. brb2022.08.3 
            sigma = par.forc.synthperiod * std_fact_in_propmat_sourc1_gauss; % Width of impulse.
            amp = max(predat_sp_PSV(:,2)); % Magnitude of direct S arrival. 
            p_simp = amp * exp((-1/(2*sigma^2))*tt_sp.^2); % parent pulse, simplified. 
% % %             figure(1); clf; hold on; 
% % %             plot(tt_sp, predat_sp_PSV./(max(max(predat_sp_PSV))) )
% % %             plot(tt_sp, p_simp)
            predat_sp_PSV(:,2) = p_simp; 
        end
        %%%

        Sp_widewind = [-50 10];
        inwind = (tt_sp >= Sp_widewind(1)) & (tt_sp < Sp_widewind(2)); 
        % crop
        predat_sp_PSV = predat_sp_PSV(inwind,:);
        tt_sp = tt_sp(inwind);

        if ~any(inwind)
            tt_sp = []; predat_sp_PSV = nan; % NO GOOD DATA
        end

        %% migrate from time to depth
        %     fprintf(' Migrating RF(z) to RF(t) using %s, rayp = %.3f s/deg\n',par.datprocess.CCP.migratemodel,par.datprocess.CCP.rayp_S)
        simp_final_mod = struct('z',final_model.Z,'VS',final_model.VSav,'VP',final_model.VPav);
        zz_mig = migrate_PorS_conv(tt_sp,simp_final_mod,rayp_sdeg2skm(rayp),0,'Sp');

        % fill in zz_mig "above" surface, corresponding to positive time...
        zz_mig(max(tt_sp)>=tt_sp & tt_sp>0) = -flipud(zz_mig(0>tt_sp & tt_sp>=-max(tt_sp)));
        % discard remaining nans
        kill = isnan(zz_mig);
        zz_mig(kill) = [];
        predat_sp_PSV(kill,:) = [];

        % convert RF_t to RF_z 
        RF_P  = interp1(zz_mig,predat_sp_PSV(:,1),final_predata.RF_Sp_ccp.zz,'linear',0);
        RF_SV  = interp1(zz_mig,predat_sp_PSV(:,2),final_predata.RF_Sp_ccp.zz,'linear',0);

        % taper off daughter = P component
        Zwin = par.datprocess.CCP.Zwin.def;
        taperz = par.datprocess.CCP.taperz;
        RF_Pw = flat_hanning_win(final_predata.RF_Sp_ccp.zz,RF_P,Zwin(1)-taperz/2,Zwin(2)+taperz/2,taperz);    


    else
        final_predata.RF_Sp_ccp.PSV = NaN;    
    end
    % -----------------  PUT INTO DATA STRUCTURE  -----------------
    final_predata.RF_Sp_ccp.PSV = [RF_Pw,RF_SV]; 
    
end

%% ===================  HK-stack: put in values  ====================
if any(strcmp(pdtyps(:,1),'HKstack'))
    
    fm = final_model; % Temporary final model. Put some posterior values in for calculating HK stack with hk_forward_model
    fm.vpvs = median(posterior.vpvs); 
    fm.rho = fm.rhoav; 
    fm.zmoh = fm.zmohav; % median(posterior.zmoh); 
    fm.VS = final_model.VSav; 
    fm.VP = final_model.VPav; 
    fm.z = fm.Z; % this is why we should always use cammel case or something else consistent! 
    fm.vpvs = fm.vpvsav;
    fm.Panis = fm.Panisav; 
    fm.Sanis = fm.Sanisav; 
    
    [final_predata, par] = hk_forward_model(...
        par, fm, final_predata, pdtyps, ...
        'insistRerun', true); 
    
    plot_HK_stack(final_predata.HKstack_P.Hgrid,...
                  final_predata.HKstack_P.Kgrid,...
                  final_predata.HKstack_P.Esum ,...
                  'model_vpvs',final_predata.HKstack_P.K,...
                  'model_zmoh',final_predata.HKstack_P.H,...
                  'title','HK stack based on inverted model',...
                  'saveStr',[par.res.resdir '/HK_from_inverted_model.pdf']); 
                          

end



end

