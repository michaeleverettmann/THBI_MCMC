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
        elseif strcmp(pdt{1},'HKstack') 
            trudata.(dtype)(itr).dof = trudata.(dtype)(itr).Nobs;
        elseif any(strcmp({'BW','RF'},pdt{1}))
            if strcmp(pdt{3},'ccp')
                trudata.(dtype)(itr).dof = length(trudata.(dtype)(itr).zz)*trudata.(dtype)(itr).dof_per_z;
            else
                trudata.(dtype)(itr).dof = mean([scdofcalc(trudata.(dtype)(itr).PSV(:,1)),...
                               scdofcalc(trudata.(dtype)(itr).PSV(:,2))]);
            end
        end
    end
    
    
end

if isfield(trudata, 'HKstack_P')
    trudata.HKstack_P.Nobs = par.mod.data.deg_of_freedom.h_kappa; 
    trudata.HKstack_P.dof  = par.mod.data.deg_of_freedom.h_kappa; 
    warning('bb2021.12.22 setting hk structure: Nobs, dof to a constant from par file.')
end

end

