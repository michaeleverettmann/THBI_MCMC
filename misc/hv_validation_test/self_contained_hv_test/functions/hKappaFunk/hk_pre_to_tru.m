function [trudata] = hk_pre_to_tru(predata,trudata)
% Does nothing on brb2022.03.01
% Trivial function to take estimated hk stack and put that in the trudata. 
% Because the best HK stack is probably the one that we have inverted...
%% Extract some info
     
% if isfield(predata, 'HKstack_P'); 
%     trudata.HKstack_P.HRecalc    = predata.HKstack_P.H    ; 
%     trudata.HKstack_P.KRecalc    = predata.HKstack_P.K    ; 
%     trudata.HKstack_P.EsumRecalc = predata.HKstack_P.Esum ; 
% end
% warning('Figure it out.')

end