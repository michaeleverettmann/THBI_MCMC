function c0_SAVE_OUTPUT(resdir,misfits_perchain,allmodels_perchain,par) %#ok<INUSD>
%  c0_SAVE_OUTPUT(resdir,misfits_perchain,allmodels_perchain)
% 
%  Function to save the important outputs from the inversion into the
%  resuts directory

save([resdir,'/misfits_perchain_orig'],'misfits_perchain');
save([resdir,'/allmodels_perchain_orig'],'allmodels_perchain');
save([resdir,'/par'],'par'); % Par can get modified. 


end
