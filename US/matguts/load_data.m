function [ trudata,zeroDstr,par ] = load_data( par )
%[ trudata,cheatstr ] = load_data( par )
paths = getPaths(); 

if nargin < 6
    baz = [];
end

P_win = [-5 25];
S_win = [-35 5] ;
tapertime = 2;
samprate = 10;
zeroDstr = []; % whether or not the daughter component is zeroed at the time of the parent

wd = pwd;

BWavardir = par.data.avardir;
sta = par.data.stadeets.sta;
nwk = par.data.stadeets.nwk;

seismoddir = [paths.models_seismic '/']; 

%% work out which data types to grab
allpdytp = parse_dtype_all(par);

%% ------------------------------------------------------------------
%% --------------------  BODY WAVE STACKED DATA  --------------------
if any(strcmp(allpdytp(:,1),'BW'))
    cd(BWavardir);
    %% Seismogram data
    datfiles = dir(sprintf('avar_dat_%s_%s_*%s.mat',sta,nwk,zeroDstr));

    if isempty(datfiles)
        trudata=[];
        return
    end

    clear('allavars','seazmean');
    for id = 1:length(datfiles)
        load(datfiles(id).name);
        allavars(id) = avar;
        seazmean(id) = mean(avar.seaz);
    end
    % find clusters of backazimuths
    clusts = eqcluster(seazmean,90*ones(size(seazmean)),0.5,30,10,0);

    np = 0;
    ns = 0;
    for ii = 1:length(allavars)
        avar = allavars(ii);
        clust = clusts(ii);
        Pind = strcmp(avar.phases,'P') | strcmp(avar.phases,'Ps') ;
        Sind = strcmp(avar.phases,'S') | strcmp(avar.phases,'Sp') ;

        %% BW_Ps ==> flip Z to 'up', taper, downsample, window

        if any(Pind)
            Pdat = avar.dataPSVSH(:,:,Pind); 
            Pdat_t  = flat_hanning_win(avar.tt(:,Pind),Pdat,P_win(1),P_win(2),tapertime);
            Pdat_td = downsamp(Pdat_t,unique(round(1./diff(avar.tt(:,Pind)))),samprate);
            tt_d = avar.tt(1,Pind) + [0:size(Pdat_td,1)-1]'/samprate;
            tt_w = P_win(1) + [0:(diff(P_win)*samprate-1)]'./samprate;
            Pdat_tdw = interp1(tt_d,Pdat_td,tt_w);
            np = np+1;
            BW_Ps(np,1) = struct('PSV',Pdat_tdw(:,1:2),'tt',tt_w,...
                                 'clust',clust,'rayp',avar.rayp(Pind),'seaz',avar.seaz(Pind),'gcarc',avar.gcarc(Pind),...
                                 'samprate',samprate,'nsamp',size(Pdat_tdw,1),...
                                 'Vp_surf',avar.Vp_Vs_surf(1),'Vs_surf',avar.Vp_Vs_surf(2));
        end

        %% BW_Sp ==> flip Z to 'up', taper, downsample, window
        if any(Sind)
            Sdat = avar.dataPSVSH(:,:,Sind); 
            Sdat_t  = flat_hanning_win(avar.tt(:,Sind),Sdat,S_win(1),S_win(2),tapertime);
            Sdat_td = downsamp(Sdat_t,unique(round(1./diff(avar.tt(:,Sind)))),samprate);
            tt_d = avar.tt(1,Sind) + [0:size(Sdat_td,1)-1]'/samprate;
            tt_w = S_win(1) + [0:(diff(S_win)*samprate-1)]'./samprate;
            Sdat_tdw = interp1(tt_d,Sdat_td,tt_w);
            ns = ns+1;
            BW_Sp(ns,1) = struct('PSV',Sdat_tdw(:,1:2),'tt',tt_w,...
                                 'clust',clust,'rayp',avar.rayp(Sind),'seaz',avar.seaz(Sind),'gcarc',avar.gcarc(Sind),...
                                 'samprate',samprate,'nsamp',size(Sdat_tdw,1),...
                                 'Vp_surf',avar.Vp_Vs_surf(1),'Vs_surf',avar.Vp_Vs_surf(2));
        end

    end % end loop on data files

    cd(wd);
end

%% ------------------------------------------------------------------
%% --------------------  RECEIVER FUNCTION DATA  --------------------
if any(strcmp(allpdytp(:,1),'RF'))
    for idt = find(strcmp(allpdytp(:,1),'RF'))
    % should we be pulling in CCP data?
        if strcmpi(allpdytp(idt,3),'ccp')
            fprintf('\n Extracting %s "RF" from CCP stack\n',allpdytp{idt,2})
            % addpath to CCP data
            addpath(seismoddir); % brb2022.06.21 Should not have addpath anywhere except in a0_STARTUP_BAYES.m
            slat = par.data.stadeets.Latitude;
            slon = par.data.stadeets.Longitude;

            % grab RF(z) from CCP
            % negative P polarity means positive velocity gradient (e.g.
            % moho), consistent with the forward models
            [RF_dz,zz_d] = load_CCP_RFz( [slat,slon],0,'version',par.datprocess.dat_version); 
            zz_d = double(zz_d);
            dz = unique(diff(zz_d));
            
            % window to desired depth range
            Zwin = par.datprocess.CCP.Zwin.def;
            taperz = par.datprocess.CCP.taperz;
            RF_dzw = flat_hanning_win(zz_d,RF_dz,Zwin(1)-taperz/2,Zwin(2)+taperz/2,taperz);
            
            % make synthetic parent
            RF_p = synthtrace(par.datprocess.CCP.parent_zw/2 + max(Zwin),...
                              par.datprocess.CCP.parent_zw,...
                              1,dz,'gauss',0);
            zz = [-par.datprocess.CCP.parent_zw/2:dz:max(Zwin)-dz]';
            RF_d = interp1(zz_d,RF_dzw,zz,'linear',0);
                                      % 

%             fprintf('NOTE WHAT TO DO ABOUT SURFV FOR CCP??\n');
            RF_struct = struct('PSV',[RF_d,RF_p],'zz',zz,...
                 'rayp',par.datprocess.CCP.rayp_S,...
                 'dz',dz,'nsamp',size(zz,1),...
                 'dof_per_z',1/par.datprocess.CCP.parent_zw); % say parent pulse width is 1 dof
            eval(sprintf('RF_%s_ccp = RF_struct;',allpdytp{idt,2}));
        else
            error('expecting CCP...')   ; 
        end
    end
end

%% ------------------------------------------------------------------
%% ----------------------  SURFACE WAVE DATA  ----------------------
SW_Ray_phV_structures = struct(''); % capture any rayleigh phv datasets. Might want to modify to not just be phv and Rayleigh. 
SW_Lov_phV_structures = struct(''); % capture any Love phv datasets. 
if any(strcmp(allpdytp(:,1),'SW'))

    %% Phase velocity data
%     if ~exist(seismoddir,'dir'); error('bb2021.09.27 where is seismodir?'); end % , seismoddir = regexprep(seismoddir,'~','/Volumes/zeilon'); end 
    addpath(seismoddir);

    slat = par.data.stadeets.Latitude;
    slon = par.data.stadeets.Longitude;

    %% -------- Rayleigh waves
    if any(strcmp(allpdytp(strcmp(allpdytp(:,1),'SW'),2),'Ray'))
        
        % Determine which Rayleigh wave datasets to use. 
        phv_modifiers = {}; 
        for idata = 1:size(allpdytp,1); 
            isRayPhv = strcmp(allpdytp{idata,1},'SW') && ...
                strcmp(allpdytp{idata,2},'Ray') && ...
                strcmp(allpdytp{idata,3},'phV'); % Indecies of data where we care about authors.  
            if isRayPhv; 
                phv_modifiers{end+1} = allpdytp{idata,4}; 
            end
        end 
        % phv_modifiers should either be empty of have cells with modifiers.
            
        
%         [Rperiods,RphV]  = Rph_dispcurve_latlon( slat,slon); % grab composite of AN and EQ
        % err = error_curve_EQ_latlon( 1./periods,avar.slat,avar.slon,phVerrordir);
        
        % Procedure for adding in another phase velocity (or similar)
        % dataset. 1. Put data in models_seismic folder. 2: make a script
        % that loads and plots all the data. 3. Make script that loads one
        % stations data (e.g. load_EQphV_data_lynhelm.m). 4. Modify
        % Rph_dispcurve_latlon_single_auth to accept the relevant "dataset"
        % variable and run your loading script. 5. Use Rph_dispcurve_latlon_single_auth on it below -- this is basically done for you. 
        all_periods = []'; % keep track of all periods that any data set used. 
        for idtype = [1:length(par.inv.datatypes)]; 
            thisdtype = par.inv.datatypes{idtype}; 
            if contains(thisdtype, 'SW_Ray'); 
                [Rperiods_i,RphV,SW_Ray_phV_struct]=...
                    Rph_dispcurve_latlon_single_auth(slat,slon,'dataset',thisdtype); 
                if isempty(SW_Ray_phV_struct); % Don't keep this dataset if there is no data for the station. 
                    warning(['Could not get data for %s. ',...
                        '\n Not using this datatype! ',...
                        'Removing from par.inv.datatypes.'],thisdtype); 
                    par.inv.datatypes{idtype} = []; % Can't keep datatype in par.inv.datatypes if we don't have that datatype. par.inv.datatypes stays same size throughout the loop. Later, we will clear the empty cells from par.inv.datatypes. 
                    continue; 
                end
                all_periods = [all_periods; Rperiods_i]; 
                SW_Ray_phV_structures(1).(thisdtype) = SW_Ray_phV_struct; 
            end
        end
        
        all_periods = unique(sort(all_periods)); 
        min_period = min([all_periods]); 
        max_period = max([all_periods]); 
        periods_calc = all_periods; % Just run calculations at the periods where authors made their measurements. Could optionally calculate at predefined periods and interpolate from those periods for individual datasets, but I don't thin kwe gain anything. . 
        n_periods_calc = length(periods_calc); 
        for_mod_info = struct('n_periods_calc', n_periods_calc,...
            'all_periods', all_periods,...
            'min_period',min_period,'max_period',max_period,...
            'periods_calc',periods_calc); % Information on what forward modelling needs to be done. We don't want to do forward modelling once for each seperate author. Instead, forward model once, accounting for every author simultaneously. 
        
        % Attach forward modelling structure to each Rayleigh wave structure. 
        fns = fieldnames(SW_Ray_phV_structures); 
        for ifn = [1:length(fns)]; 
            SW_Ray_phV_structures.(fns{ifn}).for_mod_info = for_mod_info; 
        end
        
    end

    %% -------- Love waves
    if any(strcmp(allpdytp(strcmp(allpdytp(:,1),'SW'),2),'Lov'))
        
        % Determine which Love wave datasets to use.
        phv_modifiers = {}; 
        for idata = 1:size(allpdytp,1); 
            isLovPhv = strcmp(allpdytp{idata,1},'SW') && ...
                strcmp(allpdytp{idata,2},'Lov') && ...
                strcmp(allpdytp{idata,3},'phV'); % Indices of data where we care about authors.  
            if isLovPhv; 
                phv_modifiers{end+1} = allpdytp{idata,4}; 
            end
        end 
        % phv_modifiers should either be empty or have cells with modifiers.

        all_periods = []'; % keep track of all periods that any dataset used. 
        for idtype = [1:length(par.inv.datatypes)]; 
            thisdtype = par.inv.datatypes{idtype}; 
            if contains(thisdtype, 'SW_Lov'); 
                [Lperiods_i,LphV,SW_Lov_phV_struct]=...
                    Lph_dispcurve_latlon_single_auth(slat,slon,'dataset',thisdtype); 
                if isempty(SW_Lov_phV_struct); % Don't keep this dataset if there is no data for the station. 
                    warning(['Could not get data for %s. ',...
                        '\n Not using this datatype! ',...
                        'Removing from par.inv.datatypes.'],thisdtype); 
                    par.inv.datatypes{idtype} = []; % Can't keep datatype in par.inv.datatypes if we don't have that datatype. par.inv.datatypes stays same size throughout the loop. Later, we will clear the empty cells from par.inv.datatypes. 
                    continue; 
                end
                all_periods = [all_periods; Lperiods_i]; 
                SW_Lov_phV_structures(1).(thisdtype) = SW_Lov_phV_struct; 
            end
        end
        
        all_periods = unique(sort(all_periods)); 
        min_period = min([all_periods]); 
        max_period = max([all_periods]); 
        periods_calc = all_periods; % Just run calculations at the periods where authors made their measurements. Could optionally calculate at predefined periods and interpolate from those periods for individual datasets, but I don't think we gain anything. 
        n_periods_calc = length(periods_calc); 
        for_mod_info = struct('n_periods_calc', n_periods_calc,...
            'all_periods', all_periods,...
            'min_period',min_period,'max_period',max_period,...
            'periods_calc',periods_calc); % Information on what forward modelling needs to be done. We don't want to do forward modelling once for each separate author. Instead, forward model once, accounting for every author simultaneously. 
        
        % Attach forward modelling structure to each Love wave structure. 
        fns = fieldnames(SW_Lov_phV_structures); 
        for ifn = [1:length(fns)]; 
            SW_Lov_phV_structures.(fns{ifn}).for_mod_info = for_mod_info; 
        end
        
    end

    %% -------- Rayleigh HV ratios
	if any(strcmp(allpdytp(strcmp(allpdytp(:,1),'SW'),2),'HV'))
        [ HVratios,HVstds,HVperiods] = disp_curve_HV( round_level([slat,slon],0.1),0,paths.models_seismic,10,1.5);

        [HVperiods,iT] = sort(HVperiods);
        HVratios = HVratios(iT);
        HVstds = HVstds(iT);

        gdT = ~isnan(HVratios);
        HVperiods = HVperiods(gdT);
        HVratios = HVratios(gdT);
        HVstds = HVstds(gdT);

        if ~isempty(HVratios)
            SW_HV = struct('periods',HVperiods,'HVr',HVratios,'sigma',HVstds);
        else
            SW_HV=[];
        end
        
        
        all_periods = unique(sort(HVperiods)); 
        min_period = min([all_periods]); 
        max_period = max([all_periods]); 
        periods_calc = all_periods; % No need for extra periods here unless we are bringing in more HV datasets. 
        n_periods_calc = length(all_periods); 
        for_mod_info = struct('n_periods_calc', n_periods_calc,...
            'all_periods', all_periods,...
            'min_period',min_period,'max_period',max_period,...
            'periods_calc',periods_calc); % Information on what forward modelling needs to be done. We don't want to do forward modelling once for each seperate author. Instead, forward model once, accounting for every author simultaneously. 
        SW_HV.for_mod_info = for_mod_info; 
    end

end

%% ------------------------------------------------------------------
%% --------------------  H-K stack  --------------------
if any(strcmp(allpdytp(:,1),'HKstack'))
	fprintf('\n pulling out EARS H-K stack\n')
    EARSdir = [seismoddir,'US_EARS'];
    try; 
        % Get hkstack initially used for EARS. 
        hkstack = load(sprintf('%s/EARS_HKStack_%s_%s.mat',EARSdir,nwk,sta));
        hkstack = hkstack.hkstack; 
        
        plot_HK_stack(hkstack.H, hkstack.K, hkstack.Esum,...
            'title', 'HK stack from EARS',...
            'saveString', sprintf('%s/HK_stack_EARS.pdf',par.res.resdir)); 
        
        % Get wavesforms used in EARS to make KK stack. 
        rfWaves = load(sprintf('%s/IRIS_EARS//Ears/gauss_2.5/%s.%s/rfArr.mat',EARSdir,nwk,sta)); 
        rfWaves.tt = rfWaves.tt'; 
        
        fprintf('Number of receiver functions loaded: %1.0f\n',size(rfWaves.rf,2)); 
        
        %%% Make my own HK stacks. Seems like phase weighting in IRIS EARS make inversion unstable... brb2022.04.04
        
        fprintf('\nbrb2022.02.04: TODO synthetic starting HK stack: making first hk stack using 0.5, 0.3, 0.2 weights. Also using 3.5 as average s crustal velocity, similar to Zhu and Kanamori 2000. These should be parameters. \n')
        new_h = [par.mod.crust.hmin-1:0.25:par.mod.crust.hmax+1]'; % No need to go outside prior bounds. Those models will be rejected anyway. However, go just slightly outside model bounds to increase odds of code stability. Go outside prior bounds by at least a few steps, or the prior maximum value might not ever be reached! brb2022.07.27.  
        new_k = [par.mod.crust.vpvsmin-0.02:0.005:par.mod.crust.vpvsmax+0.02]'; 
        Esum2 = zeros(length(new_k), length(new_h)); 
        Nobs = 0; 
        for irf = [1:size(rfWaves.rf,2)]; 
            [Esum2Irf, h2, k2] = HKstack(rfWaves.rf(:,irf), rfWaves.tt', ...
                rfWaves.rayParmSecDeg(irf)/111.1949, [0.5, 0.3, 0.2], ...
                3.5,new_h,new_k'); 
            Esum2 = Esum2 + Esum2Irf'; 
            Nobs = Nobs + 1; % Could just use length of rfWaves, but we might start excluding some waveforms at some point. 
        end
        
        plot_HK_stack(new_h, new_k, Esum2,'title',...
            'HK stack with vs=3.5, xi=0. Using to determine starting H and K'); 
        %%% End making my own HK stacks for starting model. 
        
    catch e 
        fprintf('\n%s\n',getReport(e)); 
        error('bb2021.10.01 Dont have h-k stack OR MAYBE WAVEFORMS for station %s.%s. Should make code to handle this situation', nwk, sta)
    end    
    if isempty(hkstack.Nobs)
        fprintf('\t Not sure how many EQ in EARS obs - assigning 100\n');
        hkstack.Nobs = 100;
    end
    HKstack_P           = struct(''); 
    HKstack_P(1).Esum      = Esum2          - mingrid(Esum2         );
    HKstack_P(1).Esum_orig = HKstack_P.Esum - mingrid(HKstack_P.Esum); % Duplicate HK stack. Might end up replacing other info with my own HK stack. 
    HKstack_P(1).H_orig    = new_h; % HKstack_P.H; 
    HKstack_P(1).K_orig    = new_k; % HKstack_P.K; 
    HKstack_P(1).H         = new_h; % HKstack_P.H; 
    HKstack_P(1).K         = new_k; % HKstack_P.K; 
    HKstack_P(1).waves     = rfWaves; 
    HKstack_P(1).Nobs      = Nobs; 
    
    HKstack_P(1).E_by_Esuper_max = ...
        hk_maximum_possible_value(rfWaves.rf, rfWaves.tt); % Not finished. Try to estimate the highest possible value that could ever be achieved in HK stack. 
    
end



%% collate
if ~exist('BW_Ps','var'),       BW_Ps = struct([]); end
if ~exist('BW_Sp','var'),       BW_Sp = struct([]); end
if ~exist('SW_Ray_phV','var'),  SW_Ray_phV = []; end
if ~exist('SW_Lov_phV','var'),  SW_Lov_phV = []; end
if ~exist('SW_HV','var'),       SW_HV = []; end
if ~exist('RF_Ps','var'),       RF_Ps = struct([]); end
if ~exist('RF_Sp','var'),       RF_Sp = struct([]); end
if ~exist('RF_Ps_ccp','var'),   RF_Ps_ccp = struct([]); end
if ~exist('RF_Sp_ccp','var'),   RF_Sp_ccp = struct([]); end
if ~exist('HKstack_P','var'),   HKstack_P = struct([]); end


trudata = struct('BW_Ps',BW_Ps,'BW_Sp',BW_Sp,...
                 'RF_Ps',RF_Ps,'RF_Sp',RF_Sp,...
                 'RF_Ps_ccp',RF_Ps_ccp,'RF_Sp_ccp',RF_Sp_ccp,...
                 'HKstack_P',HKstack_P,...
                 'SW_Ray_phV',SW_Ray_phV,'SW_Lov_phV',SW_Lov_phV,'SW_HV',SW_HV);

% Loop over SW Rayleigh phV datasets and add them to trudata if they exist.
fns = fieldnames(SW_Ray_phV_structures); 
for ifn = [1:length(fns)]; 
    thisfn = fns{ifn}; 
    trudata.(thisfn) = SW_Ray_phV_structures.(thisfn); 
end

fns = fieldnames(SW_Lov_phV_structures); 
for ifn = [1:length(fns)]; 
    thisfn = fns{ifn}; 
    trudata.(thisfn) = SW_Lov_phV_structures.(thisfn); 
end

% If we could not use some datatype, clear that from par.inv.datatypes. 
par.inv.datatypes = par.inv.datatypes(~cellfun(@isempty, par.inv.datatypes)); % Remove datatypes where we didn't have data (those should be empty now for rayleigh waves). 

cd(wd);

