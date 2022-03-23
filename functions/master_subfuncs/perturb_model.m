function [model1,ptbnorm,ifpass,p_bd,Pm_prior1,ptb,modptb,nchain,breakTrue]...
    = perturb_model(model, Pm_prior, ptbOrig, ii, par, temp, Kbase,nchainOrig)
%% ===========================  PERTURB  ===========================

ptb = ptbOrig; 
nchain = nchainOrig; 

%% ===========================  PERTURB  ===========================
    if ii==1 % no perturb on first run
        model1 = model; ptbnorm = 0; ifpass = 1; p_bd = 1; Pm_prior1 = Pm_prior; 
        modptb = nan; breakTrue = false; % In case not assigned in else. 
        ptb{ii,1} = 'start';
    else
        [model1, ptb{ii,1}, p_bd        ] = b2_PERTURB_MODEL(model,par,temp);

        % brb2022.03.14 Temporary code
%         if strcmp(ptb{ii,1},'Moho_h'); 
%             fprintf('\nPerturbed moho = %2.1f km\n',model1.zmoh)
%         end
        
        ifpass = a1_TEST_CONDITIONS( model1, par, par.inv.verbose  );
        breakTrue=false; 
		if p_bd==0, if par.inv.verbose, fprintf('  nope\n'); end; breakTrue=true; end % No need to return if breakTrue=true. We still want other variables to be defined. 
		if ~ifpass, if par.inv.verbose, fprintf('  nope\n'); end; breakTrue=true; end

        Pm_prior1 = calc_Pm_prior(model1,par);  % Calculate prior probability of new model
		[ modptb ] = calc_Vperturbation( Kbase.modelk,model1);

        ptbnorm = 0.5*(norm(modptb.dvsv) + norm(modptb.dvsh)) + norm(modptb.dvpv);
		if par.inv.verbose, fprintf('    Perturbation %.2f\n',ptbnorm); end
    end

    % quickly plot model (if not in parallel and if verbose)
    plot_quickmodel(par,model,model1)

    nchain  = kchain_addcount( nchain,ptbnorm,par );
    
end