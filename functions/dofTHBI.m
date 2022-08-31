function trudata = dofTHBI(trudata,par)
% trudata = dofTHBI(trudata)
% 
% Function to compute "N", or degrees of freedom, for each datatype:
%  # of periods if SW, 
%  # of degrees of freedom if BW or RF time series
%  some function of depth range if RF ccp
% 
dtypes = fieldnames(trudata);
for id = 1:length(dtypes)
    dtype = dtypes{id};
    pdt = parse_dtype(dtype);
    
    for itr = 1:length(trudata.(dtype))
        if strcmp(pdt{1},'SW') 
            trudata.(dtype)(itr).dof = length(trudata.(dtype)(itr).periods);
% % %             if max(trudata.(dtype)(itr).periods) > 100; 
% % %                 trudata.(dtype)(itr).dof = 6; 
% % %             else; 
% % %                 trudata.(dtype)(itr).dof = 5; 
% % %             end
% % %             warning('Testing low DOF surface waves!'); 
        elseif strcmp(pdt{1},'HKstack') 
%             trudata.(dtype)(itr).dof = trudata.(dtype)(itr).Nobs;
            trudata.HKstack_P.Nobs = par.mod.data.deg_of_freedom.h_kappa; 
            trudata.HKstack_P.dof  = par.mod.data.deg_of_freedom.h_kappa; 
            warning('Manually assigning HK stack degree of freedom'); 
        elseif any(strcmp({'BW','RF'},pdt{1}))           
            % Remove starting and ending parts of time series where they == 0. Those are artificial. 
            
            if strcmp(pdt{2}, 'Sp'); 
                daughter = trudata.(dtype)(itr).PSV(:,1); 
                parent   = trudata.(dtype)(itr).PSV(:,2); 
            elseif strcmp(pdt{2}, 'Ps'); 
                daughter = trudata.(dtype)(itr).PSV(:,2); 
                parent   = trudata.(dtype)(itr).PSV(:,1); 
            end

            if strcmp({'RF'},pdt{1}); 
                neq0d = find(daughter ~= 0); % Daughter is not zero
                start_ind = neq0d(1); 
                end_ind   = neq0d(end); 
                daughter  = daughter(start_ind:end_ind,:); 
%                 multip = ones(size(daughter)); 
%                 multip(1:80) = linspace(3,1,80); 
%                 daughter = daughter .* multip; 
                dofd = scdofcalc(daughter); 
                trudata.(dtype)(itr).dof = dofd;
            end

            if strcmp({'BW'},pdt{1}); 
                error('brb2022.08.23 Need to decide how to handle degrees of freedom for body waves if they arent receiver functions'); 
                trudata.(dtype)(itr).dof = mean([scdofcalc(daughter),...
                               scdofcalc(parent)]); % brb2022.08.23 .This should help get started. Its what Zach originally used. The problem with this is that the parent pulse isnt so important, yet contributes so many degrees of freedom!
% % %                 neq0p = find(parent   ~= 0); % Daughter is not zero
% % %                 start_ind = neq0p(1); 
% % %                 end_ind   = neq0p(end); 
% % %                 parent    = parent  (start_ind:end_ind,:); 
% % %                 dofp = scdofcalc(parent); 
            end
            
            
            
            set_dof_manually = strcmp(pdt{2}, 'Sp') && ...
                isfield(par.mod.data.deg_of_freedom, 'Sp') && ...
                ~isnan(par.mod.data.deg_of_freedom.Sp); % Sort of temporary, for using Hopper ccp dataset. brb2022.08.30. 
            
            if set_dof_manually
                dof_old = trudata.(dtype)(itr).dof; 
                dofd = par.mod.data.deg_of_freedom.Sp; 
                trudata.(dtype)(itr).dof = dofd; 
                fprintf(['\nSetting Sp degree of freedom manually to %1.2f from bayes_inv_parms.m.\n',...
                    'If you arent using the Hoper dataset, reconsider this number. \n',...
                    'For this stations specific receiver function, the Silver and Chan method\n',...
                    'Would have given you degree of freedom = %1.3f.\n',...
                    'You can manually pick a number that seems representative of all your stations.\n'], dofd, dof_old); 
            end % 
                

        end
    end
    fprintf('Degrees of freedom = %7.3f for datatype %s\n',...
        trudata.(dtype)(itr).dof, dtype)
end

end
