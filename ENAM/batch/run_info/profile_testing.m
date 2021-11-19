mpiprofile on; 
% profile ON
parfor i=[1,2]; 
    junk()
end
mpiprofile off; 
% profile OFF
stats = mpiprofile('info')
% p = profile('info'); 
save myprofiledata stats
% mpiprofile viewer
mpiprofile('viewer',stats) 

function junk()

for j=1:10000000; 
    if mod(j,2); 
        jk = j+1; 
    else
        jk = j + 1 ;
    end
end
    
end