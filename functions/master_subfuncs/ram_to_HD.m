function [results]=ram_to_HD(paths, mainDir, nwk, sta);
% Copy relevant inversion results from RAM to main hard disk. 
% Then remove this stations ram folder and cd into this stations results
% folder. 
% Paths is structure containing relevant paths. 
% chainstr: is string used for indicating one chain's ID in file names. 
% mainDir: is where you want to move things to on hard disk. For example, the folder where you started the inversion. 

% How full is the ram drive compared to how big we made it? 
[~, res_du] = system(['du -h -d 0']); 
fprintf('\nRam in use before clearing ram in %s: %s', paths.ramDrive, res_du)

% Move important things that we want to keep. 
[status, resultVel] = system(  [...
    'mv ' paths.ramDrive '/' nwk '_' sta '/' '*vel_profile* ' mainDir]); % Includes finalmodRvel_profile.mat
[status, resultDiary] = system([...
    'mv ' paths.ramDrive '/' nwk '_' sta '/diary_* '          mainDir]); 

% The rest of things on ram should be junk. But keep them in a folder on hard drive just in case. 
[~, res_du] = system(['du -h -d 0']); 
sprintf('Will transfer this much remaining junk from ram %s: %s', pwd, res_du); 

system(['mkdir -p ' mainDir '/excess_junk_from_ram']); 
system([...
    'mv ' paths.ramDrive '/' nwk '_' sta '/* ' mainDir ...
    '/excess_junk_from_ram']); 

% Now there shouldn't be anything left on this stations ram folder. 
% Verify that. 
[~,duRam] = system(sprintf('du -h %s',paths.ramDrive)); 
fprintf('\n---------\nDisk usage of ram drive after moving files off of it: \n%s\n---------\n', duRam )

fprintf('\n\nCHANGING DIRECTORY FROM\n    %s \n    TO\n    %s\n\n',pwd, mainDir);
cd(mainDir); 

results = struct('resultvel', resultVel, 'resultDiary', resultDiary); 

end