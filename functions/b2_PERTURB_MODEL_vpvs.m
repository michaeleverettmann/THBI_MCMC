% Comes from b2_PERTURB_MODEL.m. 
% Because this code gets used twice. 
% But the code is different enough each time that it's not easy to use it
% as a funciton. And the return and continue doesn't work inside the
% function. I haven't deleted this just because it could be useful later. 

% % % % function [std, ptb, V0, V1, model, vma, vmi, should_continue, should_return...
% % % %     ] = b2_PERTURB_MODEL_vpvs(temp, par, model)
% % % % % Factored out of b2_PERTURB_MODEL.m because this comes up multiple times
% % % % % in that script. brb2022.08.01
% % % % should_continue = false; 
% % % % should_return   = false; 
% % % % 
% % % % std = temp.*par.mod.crust.vpvsstd; % get std of perturbation
% % % % % if std==0, should_continue = true; ; end % don't perturb if no perturbation
% % % % ptb = 'crust_vpvs';
% % % % 
% % % % V0 = model.crustmparm.vpvs;
% % % % V1 = V0 + random('norm',0,std,1); % calc. random perturbation
% % % % 
% % % % model.crustmparm.vpvs = V1; % insert perturbed val
% % % % if par.inv.verbose, fprintf('    Changed crustal vpvs from %.2f to %.2f\n',V0,V1); end
% % % % 
% % % % % within bounds?
% % % % vma = par.mod.crust.vpvsmax;
% % % % vmi = par.mod.crust.vpvsmin;
% % % % if V1>vma || V1<vmi, P_bd = 0; should_return = true, end
% % % % 
% % % % end
% % % % 
% % % % 
% % % % % % %                     [~, ptb, ~, ~, model, ~, vmi, ~, should_return...
% % % % % % %                         ] = b2_PERTURB_MODEL_vpvs(temp, par, model)
% % % % % % %                     if should_continue; continue; end
% % % % % % %                     if should_return  ; return  ; end