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
% % % 
% % % pb = java.lang.ProcessBuilder(["echo Hello world!"]);
% % % proc = pb.start;  % now proc is a java.lang.Process object
% % % stdin = proc.getOutputStream;  % Here's how you send commands to the process
% % % stdout = proc.getInputStream;  % And here's how you get its output
% % % stderr = proc.getErrorStream;
% % % 
% % % 
% % % 
% % % % % runtime=java.lang.Runtime.getRuntime();
% % % % % process=runtime.exec('echo hello');
% % % % % status=process.waitFor();
% % % 
p = java.lang.ProcessBuilder({'/bin/bash', '-c', [' ./' execfile] }).start();
reader = java.io.BufferedReader(java.io.InputStreamReader(p.getInputStream()));
str = char(java.util.Scanner(reader).useDelimiter('\A').next())

% % % addpath('/Users/brennanbrunsvik/MATLAB/shell/jsystem'); 
% % % 
% % % 
% % % % global jsystem_path 
% % % % jsystem_path = '/Applications/MATLAB_R2021a.app/bin:/Users/brennanbrunsvik/anaconda3/bin:/Users/brennanbrunsvik/anaconda3/bin:/Users/brennanbrunsvik/anaconda3/condabin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/opt/X11/bin:/Applications/GMT-6.3.0.app/Contents/Resources/bin:/usr/local/sac/bin';
% % % 
execfile = 'TegPs.cmd'; 
[res ,cmdout ,err ] = jsystem([...
    'cd /Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/SYNTHETICS ; '...
    './' execfile], '/bin/bash -c') 

