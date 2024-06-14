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
%             trudata.(dtype)(itr).dof = length(trudata.(dtype)(itr).periods);
            if max(trudata.(dtype)(itr).periods) > 100; 
                trudata.(dtype)(itr).dof = 6; 
            else; 
                trudata.(dtype)(itr).dof = 5; 
            end
            warning('Testing low DOF surface waves!'); 
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
        end
    end
    fprintf('Degrees of freedom = %7.3f for datatype %s\n',...
        trudata.(dtype)(itr).dof, dtype)
end

end

%%% Scratch code for some tests, if needed. 
% pnew = linspace(6, 167, 10)'; 
% vnew = interp1(trudata.('SW_Ray_phV')(1).periods, trudata.('SW_Ray_phV')(1).phV, pnew);      
% scdofcalc(vnew)
% pnew, vnew
% p = [    6.0000
%     7.0298
%     8.2362
%     9.6498
%    11.3059
%    13.2463
%    15.5197
%    18.1833
%    21.3040
%    24.9603
%    29.2442
%    34.2632
%    40.1437
%    47.0333
%    55.1055
%    64.5630
%    75.6436
%    88.6260
%   103.8365
%   121.6575
%   142.5370
%   167.0000]; 
% v = [2.6956
% 2.7483
% 2.8123
% 2.8953
% 2.9229
% 2.9920
% 3.0446
% 3.1129
% 3.2402
% 3.3833
% 3.5604
% 3.7055
% 3.8511
% 3.8948
% 3.9973
% 4.0343
% 4.0993
% 4.1259
% 4.1924
% 4.2146
% 4.2869
% 4.3986];
