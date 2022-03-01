function [predata,predata0,laymodel1,fail_chain,breakTrue,newK,SW_precise] = forward_model(...
    ID,ptb,ii,predataOrig,model1,par,TD,laymodel1,fail_chainOrig,ptbnorm,newKOrig,SW_preciseOrig)

% Mostly these are variables that might be redefined in this function. 
predata=predataOrig;
predata0=predata; % save orig.
newK=newKOrig; 
SW_precise=SW_preciseOrig; 
breakTrue=false; 
fail_chain=fail_chainOrig; 

%% ===========================  FORWARD MODEL  ===========================
	% don't re-calc if the only thing perturbed is the error, or if there
	% is zero probability of acceptance!
    if ~strcmp('sig',ptb{ii}(1:3)) || isempty(predata)
        % make random run ID (to avoid overwrites in parfor)
% 		ID = [chainstr,num2str(ii,'%05.f'),num2str(randi(99),'_%02.f')];

        try
            [predata,laymodel1] = b3__INIT_PREDATA(model1,par,TD.Value,0 );
            predata = b3_FORWARD_MODEL_BW(       model1,laymodel1,par,predata,ID,0 );
            predata = b3_FORWARD_MODEL_RF_ccp(   model1,laymodel1,par,predata,ID,0 );
            predata = b3_FORWARD_MODEL_SW_kernel(model1,Kbase,par,predata );
        catch
            fail_chain=fail_chain+1;
            fprintf('Forward model error, failchain %.0f\n',fail_chain);  breakTrue=true; return
        end

        % continue if any Sp or PS inhomogeneous or nan or weird output
        if ifforwardfail(predata,par)
            fail_chain=fail_chain+1;
            fprintf('Forward model error, failchain %.0f\n',fail_chain);  breakTrue=true; return
        end

        % process predata - filter/taper etc.
        for idt = 1:length(par.inv.datatypes)
            predata = predat_process( predata,par.inv.datatypes{idt},par);
        end

		% Explicitly use mineos + Tanimoto scripts if ptb is too large
        if ptbnorm/par.inv.kerneltolmax > random('unif',par.inv.kerneltolmed/par.inv.kerneltolmax,1,1) % control chance of going to MINEOS
            newK = true;
            [ predata,SW_precise ] = b3_FORWARD_MODEL_SW_precise( model1,par,predata,ID );
        end

    end % only redo data if model has changed

%      plot_TRUvsPRE(TD.Value,predata);

    % continue if any Sp or PS inhomogeneous or nan or weird output
    if ifforwardfail(predata,par)
        fail_chain=fail_chain+1; ifpass=0;
        fprintf('Forward model error, failchain %.0f\n',fail_chain);  breakTrue=true; return
    else
        fail_chain = 0;
    end
    
    
end