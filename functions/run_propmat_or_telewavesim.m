function [traces,tt,status,cmdout] = run_propmat_or_telewavesim(p_or_t,LAYmodel,ID,ph,samprate,inc,ray_parm,synthperiod,nsamps,cutf,sourc)
arguments
    p_or_t {mustBeMember(p_or_t, {'propmat', 'telewavesim'})}
    LAYmodel
    ID string = "example"
    ph string = "Ps"
    samprate = [] 
    inc = [] 
    ray_parm = [] 
    synthperiod = 1
    nsamps = [] 
    cutf = []
    sourc string = "gauss"
end

if strcmp(p_or_t, 'propmat'); 
    [traces,tt,status,cmdout] = run_propmat(LAYmodel,ID,ph,samprate,inc,synthperiod,nsamps,cutf,sourc); 
elseif strcmp(p_or_t, 'telewavesim'); 
    [traces,tt,status,cmdout] = run_telewavesim(LAYmodel,ID,ph,samprate,ray_parm,synthperiod,nsamps,cutf,sourc); 
end 

end