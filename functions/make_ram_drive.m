function [ramDir]=make_ram_drive(options); 
    arguments
        options.linuxRamBase = '/dev/shm'; 
        options.linuxRam = '/dev/shm/ramDriveTHBI';  
        options.macRam = '/Volumes/ramDriveTHBI';
        options.ramMb = 1000; % How many mb to use for one stations ram. Should only apply on Mac. bb2021.12.03. 
    end
%% Function to check whether a ram drive is present, make one if it is not, and return that ram drives position. 
% Needed for using temp files really fast throughout THBI. 
% linuxRamBase: Linux computers might have /dev/shm which is always a ram
%   drive. If it's there, we can make a folder in it for THBI. 
% linuxRam: If you are on a linux computer, this is where we should look for the ram drive.
% macRam: If you are on a Mac, this is where we will look for OR MAKE a ram drive. 
% options.ramMb: on a Mac only (bb2021.12.03), your ram drive will be this large. 
%   Partially tested as of 2021.12.03

if exist(options.linuxRamBase); 
    
    sprintf('Linux ram drive %s found. Making or using folder there for THBI.', options.linuxRamBase)
    if exist(options.linuxRam); 
        sprintf('Linux ram THBI drive %s found. Using it for THBI. ', options.linuxRam)
        ramDir = options.linuxRam; 
    else 
        sprintf('Making THBI ram drive %s.', options.linuxRam); 
        mkdir(options.linuxRam); 
        ramDir = options.linuxRam;
    end
        
elseif ~exist(options.linuxRamBase); 
    
    sprintf('Linux ram drive %s not present. Assuming we are using mac.', options.linuxRamBase) 
    if exist(options.macRam); 
        sprintf('Mac ram drive %s found. Reusing it (might want to clear it also...)', options.macRam)
        ramDir = options.macRam; 
    else 
        sprintf('Mac ram drive not %s found. Mounting ram drive. ', options.macRam)
        mountPoint = erase(options.macRam, '/Volumes/'); % Get rid of the volumes part. It seems to come up automatically. 
        sysStr = ['diskutil erasevolume HFS+ "' mountPoint ... % String to execute with system which mounts a mac drive. 2048 has to be multiplied to mb to get ram blocks (according to... the internet. bb2021.12.03). 
            '" `hdiutil attach -nomount ram://' ...
            num2str(round(2048*options.ramMb)) '`'];
        [status,cmdout] = system(sysStr); 
        ramDir = options.macRam; 
        if status==0; 
            sprintf('Mac ram drive %s mounted (I think)', options.macRam) 
        else; 
            error('bb2021.12.03: Mac ram drive did not mount. system status %s', num2str(status))
        end
    end
    
end
   
end