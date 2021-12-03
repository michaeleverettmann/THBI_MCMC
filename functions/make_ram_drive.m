function [ramDir]=make_ram_drive(options); 
    arguments
        options.linuxRam = '/dev/shm';  
        options.macRam = '/Volumes/ramDriveTHBI';
        options.ramMb = 1000; % How many mb to use for one stations ram. Should only apply on Mac. bb2021.12.03. 
    end
%% Function to check whether a ram drive is present, make one if it is not, and return that ram drives position. 
% Needed for using temp files really fast throughout THBI. 
% linuxRam: If you are on a linux computer, this is where we should look for the ram drive.
% macRam: If you are on a Mac, this is where we will look for OR MAKE a ram drive. 
% options.ramMb: on a Mac only (bb2021.12.03), your ram drive will be this large. 
% Partially tested as of 2021.12.03

if exist(options.linuxRam); 
    sprintf('Linux ram drive %s found. Using this for storing temp files', options.linuxRam)
    ramDir = options.linuxRam; 
elseif ~exist(options.linuxRam); 
    sprintf('Linux ram drive %s not present. Assuming we are using mac.', options.linuxRam) 
    
    % Do things to make mac ram drive. 
    if exist(options.macRam); 
        sprintf('Mac ram drive %s found. Reusing it (might want to clear it also...)', options.macRam)
        ramDir = options.macRam; 
    else 
        sprintf('Mac ram drive not %s found. Mounting ram drive. ', options.macRam)
        mountPoint = erase(options.macRam, '/Volumes/'); % Get rid of the volumes part. It seems to come up automatically. 
        sysStr = ['diskutil erasevolume HFS+ "' mountPoint ...
            '" `hdiutil attach -nomount ram://' ...
            num2str(round(2048*options.ramMb)) '`'];
        [status,cmdout] = system(sysStr); 
        ramDir = options.macRam; 

        if a==0; 
            sprintf('Mac ram drive mounted (I think)'); 
        else; 
            error('bb2021.12.03: Mac ram drive did not mount. system status %s', num2str(status))
        end
    end
end
   
end