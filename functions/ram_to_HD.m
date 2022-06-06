function [results]=ram_to_HD(paths, mainDir, nwk, sta);
% Copy relevant inversion results from RAM to main hard disk. 
% Paths is structure containing relevant paths. 
% chainstr: is string used for indicating one chain's ID in file names. 
% mainDir: is where you want to move things to on hard disk. For example, the folder where you started the inversion. 

% bb2021.12.07 I tried scp, rsync, some other options. Doing rsync once at
% the end of the inversion gets the job done. Copying several times
% throughout the inversion gets slow, so probably avoid that. 
[status, resultVel] = system(  ['rsync ' paths.ramDrive '/' nwk '_' sta '/' '*vel_profile* ' mainDir]); % Includes finalmodRvel_profile.mat
[status, resultDiary] = system(['rsync ' paths.ramDrive '/' nwk '_' sta '/diary_* '          mainDir]); 
results = struct('resultvel', resultVel, 'resultDiary', resultDiary); 

% How full is the ram drive compared to how big we made it? 
[~, res_du] = system(['du -h -d 0']); 
sprintf('Ram currently in use in %s: %s', paths.ramDrive, res_du)

% brb2022.06.03 Don't remove filres from ram. diary and inv files need to
% be copied to resdir. Then the files are later deletd off ram. 
% % % [status, resultVel] = system(  ['rm ' paths.ramDrive '/' nwk '_' sta '/' '*vel_profile* ' ]); % Includes finalmodRvel_profile.mat
% % % [status, resultDiary] = system(['rm ' paths.ramDrive '/' nwk '_' sta '/diary_* '          ]); 


end