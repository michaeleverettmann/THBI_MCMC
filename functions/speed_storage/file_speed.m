function [a] = file_speed(fnum)
% if isempty(fnum); 
%     fnum = ''; 
% end
if nargin<1; 
    fnum = ''; 
end

[fid] = fopen(['/Volumes/RAMDisk/Junk/junk_' fnum '.txt'], 'w'); 
fprintf(fid, 'Junk') ; 
a = fclose(fid); 

% /Volumes/RAMDisk/Junk

% Simple command to make a 1 mb ram disk
% diskutil erasevolume HFS+ "RAMDisk" `hdiutil attach -nomount ram://2048`
ramMb = 200; % integer number of Mb to dedicate to ram volume. 
mountPoint = ['RAMDisk'] % Folder where ram volume will be placed. 
sysStr = ['diskutil erasevolume HFS+ "' mountPoint '" `hdiutil attach -nomount ram://' num2str(round(2048*ramMb)) '`'] % command to make ram disk

end