function [results]=ram_to_HD(paths, mainDir, nwk, sta);
% Copy relevant inversion results from RAM to main hard disk. 
% Paths? is structure containing relevant paths. 
% chainstr: is string used for indicating one chains ID in file names. 
% mainDir: is where you want to move things to. Should be folder where you
% started the inversion. bb2021.12.03

% Could use rsync, but probably slow since we anticipate that the files
% have changed. Comparing files might be slow. 
% TODO change to copyfile if the system call is slow. Shouldn't matter
% since this is ran only every few hundred iterations. 
% rsync command use -W to avoid doing diff on files. 
[status, resultVel] = system(  ['rsync ' paths.ramDrive '/' nwk '_' sta '/' '*vel_profile* ' mainDir]); % Includes finalmodRvel_profile.mat, the 
[status, resultDiary] = system(['rsync ' paths.ramDrive '/' nwk '_' sta '/diary_* '          mainDir]); 
results = struct('resultvel', resultVel, 'resultDiary', resultDiary); 

% How full is the ram drive compared to how big we made it? 
% If this is slow, comment it out. We can check manually. 
[~, res_du] = system(['du -h -d 0']); 
sprintf('Ram currently in use in %s: %s', paths.ramDrive, res_du)

% nwk = 'US'; sta = 'CEH'

end