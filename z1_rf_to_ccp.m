function [ data ] = z1_rf_to_ccp(par,data,model,ifplot); 

tt_sp = data.RF_Sp.tt; 
rayp = data.RF_Sp.rayp; 
predat_sp_PSV = data.RF_Sp.PSV; 

% brb2022.08.17 The hard coded variables come from an example from trudata variable that came out of an inversion on real data. 
dz = 0.5; % 
zz = [-15:dz:300]';
data.RF_Sp_ccp.zz = zz; 
dof_per_z=1/par.datprocess.CCP.parent_zw; 
data.RF_Sp_ccp = struct('zz',zz,'rayp',data.RF_Sp.rayp,...
    'dz',dz,'nsamp',length(zz),'dof_per_z',dof_per_z); 

%% migrate from time to depth
% brb copy pasted from b3_FORWARD_MODEL_RF_ccp. Should re-write into function. 
zz_mig = migrate_PorS_conv(tt_sp,model,rayp_sdeg2skm(rayp),0,'Sp');

% fill in zz_mig "above" surface, corresponding to positive time...
zz_mig(max(tt_sp)>=tt_sp & tt_sp>0) = -flipud(zz_mig(0>tt_sp & tt_sp>=-max(tt_sp)));
% discard remaining nans
kill = isnan(zz_mig);
zz_mig(kill) = [];
predat_sp_PSV(kill,:) = [];

% convert RF_t to RF_z 
try 
    RF_P   = interp1(zz_mig,predat_sp_PSV(:,1),data.RF_Sp_ccp.zz,'linear',0);
    RF_SV  = interp1(zz_mig,predat_sp_PSV(:,2),data.RF_Sp_ccp.zz,'linear',0);
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
RF_Pw = flat_hanning_win(data.RF_Sp_ccp.zz,RF_P,Zwin(1)-taperz/2,Zwin(2)+taperz/2,taperz);  

data.RF_Sp_ccp.PSV = [RF_Pw,RF_SV]; 

if ifplot; 
    
    figure(1); clf; hold on; set(gcf, 'color', 'white'); 
    tiledlayout(1,3,'TileSpacing', 'compact'); 
    
    nexttile(); hold on; grid on; box on; 
    set(gca,'ydir', 'reverse');
%     ylim(zaxlim); 
    plot(model.VS, model.Z); 
    xlabel('VS'); 
    ylabel('Depth (km)'); 
    
    nexttile(); hold on; grid on; box on; 
    set(gca,'ydir', 'reverse');
%     ylim(zaxlim); 
    plot(data.RF_Sp_ccp.PSV,data.RF_Sp_ccp.zz); ; 
    xlabel('Sp ccp'); 
    ylabel('Depth (km)'); 
    
    nexttile(); hold on; grid on; box on; 
    plot(data.RF_Sp.PSV, data.RF_Sp.tt); 
    xlabel('Sp'); 
    ylabel('Time (s)'); 
end

dtps = string(par.inv.datatypes); 
fns = string(fieldnames(data)); 

if ~any(dtps=='RF_Sp') && any(fns=='RF_Sp'); % If we had 'RF_Sp' in our data but we only wanted it for making ccp stacks, remove this from the data. 
    data = rmfield(data, 'RF_Sp'); 
end

end
