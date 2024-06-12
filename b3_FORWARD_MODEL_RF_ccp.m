function [ predata ] = b3_FORWARD_MODEL_RF_ccp( model,laymodel,par,predata,ID,ifplot )
%[ predata ] = b3_FORWARD_MODEL_RFCCP( laymodel,par,predata,ID,ifplot )
% 
%   Do forward model to calculate predicted data.
% 
% INPUTS
%   model   - model structure
%   laymodel - layered model
%   par     - parameters structure
%   predata - data structure with all datatypes 
%   ID      - unique ID for the propmat script to avoid overwriting files
%             if running in parallel.
%   ifplot  - flag with option to plot (1) or not (0)
% 
% OUTPUTS
%   predata - structure identical to input data structure, but with
%             predicted data, rather than observed data
%%
% An important component is the layerising of the model - conversion of
% continuous model into a bunch of layers, with coarseness partly
% determined by the minimum dVs in any layer (specified as an input). The
% layerised 1D model is also output from this function.

%% ===================  PREPARE DATA STRUCTURE  ===================
indtypes = fieldnames(predata); % look at all dtypes in predata
iccp = find(strcmp(indtypes,'RF_Sp_ccp'),1);
if isempty(iccp), return; end

%% ===================  LAYERISE PROFILE  ===================
% Should be already done - if not, do here
if isempty(laymodel)
    warning('brb2022.07.06 Why didnt layer model get passed to this function?'); 
    [zlayt,zlayb,Vslay,...% brb2022.07.06 Version that prevents us from needing to make assumptions here about xi, rho, eta, ...
    Vplay,rholay,xilay,philay] = ...
        layerise(model.z,model.VS,par.forc.mindV,0,...
        model.VP, model.rho, ...
        model.Sanis./100+1,model.Panis./100+1); 
    nlay = length(Vslay);
    etalay = ones(nlay,1); % eta anisotropy TODO get eta from model. eta is not in model right now, so can't yet. 
    
    laymodel = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay,'xi',xilay,'phi',philay,'eta',etalay);
    if any(isnan(laymodel.rho))
        error('NaN densities')
    end
end

% figure(1); clf; hold on; 
% % plot(laymodel.zlayt, laymodel.zlayb); 
% % laymodel.zlayt - laymodel.zlayb
% scatter(laymodel.zlayt,1+ones(size(laymodel.zlayt))); 
% scatter(laymodel.zlayb,ones(size(laymodel.zlayb))); 
% ylim([-30, 10]); 

%% ===================  SP RFs FROM PROPAGATOR MATRIX  ====================
rayp = predata.RF_Sp_ccp.rayp;
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

    [predat_sp,tt_sp] = run_propmat_or_telewavesim(par.synth.propmat_or_telewavesim, laymodel_Suse,ID,'Sp',samprate, S_inc, rayp,...
        par.forc.synthperiod,par.forc.nsamps);
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
    zz_mig = migrate_PorS_conv(tt_sp,model,rayp_sdeg2skm(rayp),0,'Sp');
    
    % fill in zz_mig "above" surface, corresponding to positive time...
    zz_mig(max(tt_sp)>=tt_sp & tt_sp>0) = -flipud(zz_mig(0>tt_sp & tt_sp>=-max(tt_sp)));
    % discard remaining nans
    kill = isnan(zz_mig);
    zz_mig(kill) = [];
    predat_sp_PSV(kill,:) = [];
    
    % convert RF_t to RF_z 
    try 
        RF_P  = interp1(zz_mig,predat_sp_PSV(:,1),predata.RF_Sp_ccp.zz,'linear',0);
        RF_SV  = interp1(zz_mig,predat_sp_PSV(:,2),predata.RF_Sp_ccp.zz,'linear',0);
    catch e 
        fprintf('\nProblem with receiver function calculation. psv waveforms disp below:\n')
        disp(predat_sp_PSV)
        error(['\nReceiver function error. Maybe nan receiver functions? ',...
            'If so, the model probably should count as junk. ',...
            'Error message was this: \n',...
            '%s\n'],getReport(e)); 
    end

    % taper off daughter = P component
	Zwin = par.datprocess.CCP.Zwin.def;
    taperz = par.datprocess.CCP.taperz;
	RF_Pw = flat_hanning_win(predata.RF_Sp_ccp.zz,RF_P,Zwin(1)-taperz/2,Zwin(2)+taperz/2,taperz);    
    
    
else
    predata.RF_Sp_ccp.PSV = NaN;    
end

% -----------------  PUT INTO DATA STRUCTURE  -----------------
predata.RF_Sp_ccp.PSV = [RF_Pw,RF_SV]; 

%% ifplot....
% ifplot = true; if ifplot; warning('Setting ifplot = true'); end; 
if ifplot
    zz  = predata.RF_Sp_ccp.zz;
    RFz = predata.RF_Sp_ccp.PSV;
   
    figure(2), clf, set(gcf,'position',[-1358 247 776 796], 'color', 'white');
    tiledlayout(1,2,'TileSpacing', 'compact'); 
        
	nexttile; cla; hold on; box on; set(gca, 'linewidth', 1.5); 
    plot(RFz(:,1),zz,'Linewidth',2)
    title('P comp (= CCP)','fontsize',19,'fontweight','bold')
    set(gca,'ydir','reverse','xlim',0.3*[-1 1],...
        'ylim',[-par.datprocess.CCP.parent_zw/2 par.datprocess.CCP.Zwin.def(2)])

	nexttile; cla; hold on; box on; set(gca, 'linewidth', 1.5); 
    plot(RFz(:,2),zz,'Linewidth',2)
    xlim(par.datprocess.Sp.Twin.def);
    set(gca,'ydir','reverse','xlim',[-1 1],...
        'ylim',[-par.datprocess.CCP.parent_zw/2 par.datprocess.CCP.Zwin.def(2)])
    title('SV comp (= parent)','fontsize',19,'fontweight','bold')
end % on ifplot

end

