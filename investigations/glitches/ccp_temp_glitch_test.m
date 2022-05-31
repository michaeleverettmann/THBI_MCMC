% clc; clear; 
% run('../../a0_STARTUP_BAYES.m'); 
% load('/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_eri/S57A_TA_dat1/all_002/chain_with_broken_something_TA.S57A_A.mat')
% load('/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_eri/S57A_TA_dat1/all_002/chain_with_broken_something_TA.S57A_L.mat')
% fig_path = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/figures/testing/forward_model_breaking/'; 
% model = model1; 

%% Use this code to save only what's necessary for testing out the broken model. Saves file space. 
% saveToPath = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/workspaces/broken_chains_for_zach/'; 
% save([saveToPath 'breaksPropmat_1.mat'], 'model1', 'par', 'trudata', 'laymodel1', 'predata', 'ID', 'predataPrev', 'Kbase', 'model'); 
%% 

nFails = 0; 
for imod = [1:10]; 
    
    if imod ~= 1; 
        [model1,ptbnorm,ifpass,p_bd,Pm_prior1,...
            ptb,modptb,nchain,breakTrue,non_acceptk]...
            = delay_reject(...
                model, Pm_prior, ptb, ii, par, temp, Kbase,nchain,...
                model1,ptbnorm,p_bd,Pm_prior1k,non_acceptk); 
    end
    
    try
%         par.forc.synthperiod = 2; 
%         par.forc.nsamps = 2048 * 1
        [predata,laymodel1] = b3__INIT_PREDATA(model1,par,trudata,1 );
    %     [predata,par] = b3_FORWARD_MODEL_BW(model1,laymodel1,par,predata,ID,0,predataPrev);
        predata = b3_FORWARD_MODEL_RF_ccp(   model1,laymodel1,par,predata,ID,1 );
    %     predata = b3_FORWARD_MODEL_SW_kernel(model1,Kbase,par,predata );
    
% Increase velocity of deep stuff. 
%         model1.VS = model1.VS + linspace(0,1,length(model1.VS))'; 
%         model1.VP = model1.VP + linspace(0,1,length(model1.VS))'; 
%         model1.rho= model1.rho+ linspace(0,1,length(model1.VS))'; 

% Try reducing velocity of shallowest sediment. 
%         nc = 10; % Number of changed indicies. 
%         model1.VS (1:nc) = model1.VS (1:nc) - linspace(.5,0,nc)'; 
%         model1.VP (1:nc) = model1.VP (1:nc) - linspace(.5,0,nc)';
%         model1.rho(1:nc) = model1.rho(1:nc) - linspace(.5,0,nc)';

    
        [ predata,SW_precise ] = b3_FORWARD_MODEL_SW_precise( model1,par,predata,ID )
        [Kbase,predata] = b7_KERNEL_RESET(model,Kbase,predata,ID,ii,par,0,SW_precise)
        success_model = true; 
    catch 
        success_model = false; 
        nFails = nFails + 1; 
        fprintf('\nFail!\n'); 
    end
        
    figure(1); clf; hold on; set(gcf, 'color', 'white'); 
    nRow = 1; nCol = 3; 

    subplot(nRow, nCol, 1); cla; hold on; set(gca,'Ydir','reverse'); box on; grid on; 
    xlabel('Vs,p'); 
    plot(laymodel1.Vs, laymodel1.zlayt); 
    plot(laymodel1.Vp, laymodel1.zlayt); 

    subplot(nRow, nCol, 2); cla; hold on; set(gca,'Ydir','reverse'); box on; grid on; 
    xlabel('Rho'); 
    plot(laymodel1.rho, laymodel1.zlayt); 

    subplot(nRow, nCol, 3); cla; hold on; set(gca,'Ydir','reverse'); box on; grid on; 
    xlabel('Xi'); 
    plot(laymodel1.xi, laymodel1.zlayt); 

%     subplot(nRow, nCol, 4); cla; hold on; set(gca,'Ydir','reverse'); box on; grid on; 
%     xlabel('phi'); 
%     plot(laymodel1.phi, laymodel1.zlayt); 
% 
%     subplot(nRow, nCol, 5); cla; hold on; set(gca,'Ydir','reverse'); box on; grid on; 
%     xlabel('eta'); 
%     plot(laymodel1.eta, laymodel1.zlayt); 
    
    fold_break = {'yes_break', 'no_break'}; 
    sub_fold = fold_break{int16(success_model)+1}; 

    sgtitle(sprintf('elastic (z): success = %1.0f',int16(success_model))); 
    
%     exportgraphics(gcf, ...
%         sprintf('%s%s/model_broke_forward_code_idmod%1.0f.pdf',...
%         fig_path,sub_fold,imod)); 

end
fprintf('\nNumber of fails: %1.0f\n',nFails)