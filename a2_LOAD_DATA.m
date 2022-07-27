function [trudata,par] = a2_LOAD_DATA(par,options)
    arguments
        par
        options.nwk = [] % Must provide if loading real data
        options.sta = [] % Must provide if loading real data
    end
% function to load data (differs depending on which synth or real)

fprintf('LOADING data\n')
global TRUEmodel %#ok<NUSED>
paths = getPaths(); 

%% -----------------------------------------------------------------
%% EXAMPLE
if any(strcmp(par.proj.name,{'EXAMPLE'}))
    
    fprintf(' > Loading example data and target model\n')
    
    global TRUEmodel TLM trudata
    
    load('target_model.mat')
    load('trudata.mat')
    
    %% -----------------------------------------------------------------
%% SYNTHETICS
elseif any(strcmp(par.proj.name,{'SYNTHETICS','test_RF_vs_BW'}))
    
    fprintf(' > Creating custom model and synthetic data\n')
    
    % make synth model
    [model,laymodel,par] = z0_SYNTH_MODEL_splinemoduse(par,0);  %close(95)
    save([par.res.resdir,'/trumodel'],'TRUEmodel');

    % make synth data
%     [ junkdata ] = z1_SYNTH_DATA(par,0);
%     [predata,laymodel1] = b3__INIT_PREDATA(TRUEmodel,par,junkdata,0 );
%     predata = b3_FORWARD_MODEL_BW(TRUEmodel,laymodel1,par,predata,'test_bb2021_12_07',0 ); 
%     predata = b3_FORWARD_MODEL_RF_ccp(   TRUEmodel,laymodel1,par,predata,'test_bb2021_12_07',0 );
    [ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format
    
    if strcmp(par.synth.noisetype,'real')
        [ trudata,par ] = z2_NOISIFY_SYNTH( trudata, par,par.synth.noise_sta_deets );
    end

    % distribute data for different processing (e.g. _lo, _cms)
    trudata = duplicate_data_distribute(trudata,par);
   
    % conv to receiver functions
    if any(strcmp(par.proj.name,{'SYNTHETICS','test_RF_vs_BW'}))
        trudata = BW_2_RF(trudata,par);
    end
    
    % Use receiver functions to make h-kappa stacks. 
    [ trudata ] = z1_SYNTH_DATA_h_kappa(par,trudata,model,0);
    
%% -----------------------------------------------------------------
%% LAB TESTING
elseif strcmp(par.proj.name,'LAB_tests')
	z0_SYNTH_MODEL_LAB_TEST(par,par.synth.model.zsed,par.synth.model.zmoh,par.synth.model.zlab,par.synth.model.wlab,par.synth.model.flab,1) ;
	save([par.res.resdir,'/trumodel'],'TRUEmodel');

	% make synth data
	[ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format
	trudata_noiseless = trudata;

	trudata = trudata_noiseless;
	if strcmp(par.synth.noisetype,'real')
	%     [ trudata,par ] = z2_NOISIFY_SYNTH_makestack( trudata, par,noise_sta_deets );
		[ trudata,par ] = z2_NOISIFY_SYNTH( trudata, par,par.synth.noise_sta_deets );
    end
    
    trudata = duplicate_data_distribute(trudata,par);
    par.data.stadeets.stadeets = struct('Latitude',[],'Longitude',[]);

%% -----------------------------------------------------------------
%% REAL DATA
else
% % % 	try 
% % %         try % Janky try catch in try catch... I cant use internet with slurm on CNSI computers? so have to have a simple solution for loading saved .mat file instead of using irisetch. 
% % %             stadeets = irisFetch.Stations('station',par.data.stadeets.nwk,par.data.stadeets.sta,'*','*'); 
% % %             save([par.proj.rawdatadir,'/stadeets.mat'], 'stadeets'); 
% % %         catch e
% % %             warning(['brb2022.04.11 Could not irisFetch stadeets. Loading saved file instead. ',...
% % %               'For this to work, you have to send the stadeets.mat file ',...
% % %               'to this cluster!!! I did this MANUALLY before.'])
% % %             stadeets = load([par.proj.rawdatadir,'/stadeets.mat']); 
% % %             stadeets = stadeets.stadeets; 
% % %         end    

        if isempty(options.nwk) || isempty(options.sta); 
            error('Provide nwk and sta to a2_LOAD_DATA.m so we know what to load if you are using real data'); 
        end
        
        fname = [par.proj.rawdatadir,'stadeets_' options.nwk '_' options.sta '.mat']; 
        
        if ~ java.io.File(fname).exists(); 
            error('stadeets file for %s.%s does not exist. Probably you need to run RUN_prep_data on this machine.'); 
        end
        
        stadeets = load(fname); 
        stadeets = stadeets.stadeets_netsta; 
            
        fns = fieldnames(stadeets);
        for ii = 1:length(fns)
            par.data.stadeets.(fns{ii}) = stadeets.(fns{ii});
        end
        
% % % 	catch e 
% % %         warning('Could not get stadeets.mat. Might need to send that directly to the cluster computer. Looking for stainfo_master.mat instead - might not work...')
% % %         fprintf('\nError report was: %s\n',getReport(e)); 
% % %         load([par.proj.rawdatadir,'/stainfo_master.mat']); 
% % %         par.data.stadeets.Latitude = stainfo(strcmp({stainfo.StationCode},par.data.stadeets.sta)).Latitude;
% % %         par.data.stadeets.Longitude = stainfo(strcmp({stainfo.StationCode},par.data.stadeets.sta)).Longitude;
% % % 	end
    
    %% LOAD THE DATA!!
    [trudata,zeroDstr,par] = load_data(par);
    par.data.stadeets.sta = [par.data.stadeets.sta,zeroDstr];
    
    if isempty(trudata), return; end
    
    % distribute data for different processing (e.g. _lo, _cms)
    trudata = duplicate_data_distribute(trudata,par);
    % set prior sigma as the mean of the data uncertainty if available
    par = prior_sigma_distribute(par,trudata);
    
    % get rid of data that wont be used in inversion - NB NEED EXACT DATA MATCH
    trudtypes = fieldnames(trudata);
    for idt = 1:length(trudtypes)
        if all(~strcmp(trudtypes{idt},par.inv.datatypes)) % no match
            fprintf('WARNING - removing %s data from trudata\n',trudtypes{idt})
            trudata = rmfield(trudata,trudtypes{idt});
            continue
        end
        if isempty(trudata.(trudtypes{idt}))
            fprintf('WARNING - No %s data in trudata\n',trudtypes{idt})
            trudata = rmfield(trudata,trudtypes{idt});
        end
    end
    
    % get rid of remaining data that is not in the right cluster
    trudtypes = fieldnames(trudata);
    for idt = 1:length(trudtypes)    
        if par.inv.BWclust~=0 && any(regexp(trudtypes{idt},'BW'))
            if isempty(trudata.(trudtypes{idt}))
                trudata = rmfield(trudata,trudtypes{idt});
                continue
            end
            fprintf('WARNING - removing %s data not in cluster %.0f\n',trudtypes{idt},par.inv.BWclust)
            trudata.(trudtypes{idt}) = trudata.(trudtypes{idt})([trudata.(trudtypes{idt}).clust]==par.inv.BWclust);
            % if getting rid of that cluster killed the data type, then
            % delete the data fieldname
            if isempty(trudata.(trudtypes{idt}))
                fprintf('THAT WAS ALL THE %s data\n',trudtypes{idt})
                trudata = rmfield(trudata,trudtypes{idt});
            end
        end
    end
    
    % get rid of datatypes from par if not (any longer) in real data
    kill = false(length(par.inv.datatypes),1);
    for idt = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{idt}; pdt = parse_dtype( dtype ); 
        % since we added in all data types above, if any not left now, it
        % is because they don't exist in the real data - should not use
        if ~isfield(trudata,par.inv.datatypes{idt}) %&& strcmp(pdt{1},'BW') 
            kill(idt) = true;
        end
    end
    par.inv.datatypes(kill) = [];

end

% save data
% save([par.res.resdir,'/trudata_ORIG'],'trudata');

%% PROCESS - window, filter data 
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    [ trudata ] = predat_process( trudata,dtype,par);
end

%% DEGREES OF FREEDOM
fprintf('Computing degrees of freedom (N) for each data type\n');
trudata = dofTHBI(trudata, par);

% plot_TRU_WAVEFORMS(trudata);
% plot_TRUvsPRE(trudata,trudata);
save([par.res.resdir,'/trudata_USE'],'trudata');

if ~exist('sta','var') %TODOEXIST bb2021.11.22 exist is SUPER slow - It's ok since you don't load data often. 
        sta = 'SYN';
end 

end

