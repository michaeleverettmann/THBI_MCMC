nWorkers = 10;
% % % 
% % % delete(gcp('nocreate'))
% % % myCluster = parcluster('local');
% % % myCluster.NumWorkers = nWorkers; 
% % % % maxWorkers = myCluster.NumWorkers; 
% % % parpool(maxWorkers); 
% % % 


allWorkers = [1:nWorkers]; 
mpiprofile on; 
parfor iwrkr = allWorkers; 
    
    for irep = [1:10]; 
%     disp(iwrkr)
%     [a,b]=system(sprintf('%s 0.3 echo Hi',timeoutPath));
    timeoutPath = '/usr/local/bin/timeout'; 
    [a,b]=system(sprintf('%s 2 sleep .2; echo Done sleeping', timeoutPath));
    fprintf('\n%s\n',b)
%     fprintf('\nFinished for worker %1.0f\n',iwrkr)
    end
end
mpiprofile off; 
% mpiStats = mpiprofile('info');
% mpiprofile('viewer', mpiStats);
% What happens if we time just the system command? 
% fhand = @()system('echo ""')
% timeit(fhand)
% I got about 0.007 seconds. Should be fast enough for what we are doing. 