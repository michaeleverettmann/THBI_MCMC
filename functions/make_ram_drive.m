function [ramDir]=make_ram_drive(options); 
    arguments
        options.linuxRamBase = '/dev/shm'; 
        options.linuxRam = '/dev/shm/ramDriveTHBI';  
        options.macRam = '/Volumes/ramDriveTHBI';
        options.ramMb = 1000; 
    end
%% Function to check whether a ram drive is present, make one if it is not, and return that ram drives position. 
% Needed for using temp files really fast throughout THBI. 
% linuxRamBase: Linux computers might have /dev/shm which is always a ram
%   drive. If it's there, we can make a folder in it for THBI. 
% linuxRam: For a specific station, use this folder. 
% macRam: If you are on a Mac, this is where we will look for OR MAKE a ram drive. 
% options.ramMb: on a Mac only (bb2021.12.03), your ram drive will be this large. 
%   Partially tested as of 2021.12.03

if exist(options.linuxRamBase); % bb2021.12.07 Fine to use slow exist() since this is called very infrequently. 
    
    fprintf('\nLinux ram drive %s found. Making or using folder there for THBI.\n', options.linuxRamBase)
    if exist(options.linuxRam); 
        fprintf('\nLinux ram THBI drive %s found. Using it for THBI. \n', options.linuxRam)
        ramDir = options.linuxRam; 
    else 
        fprintf('\nMaking THBI ram drive %s.\n', options.linuxRam); 
        mkdir(options.linuxRam); 
        ramDir = options.linuxRam;
    end
        
elseif ~exist(options.linuxRamBase); 
    
    fprintf('\nLinux ram drive %s not present. Assuming we are using mac. \n', options.linuxRamBase) 
    if exist(options.macRam); 
        fprintf('\nMac ram drive %s found. Reusing it (might want to clear it also...)\n', options.macRam)
        ramDir = options.macRam; 
    else 
        fprintf('\nMac ram drive not %s found. Mounting ram drive. \n', options.macRam)
        mountPoint = erase(options.macRam, '/Volumes/'); % Get rid of the volumes part. It seems to generate automatically. 
        sysStr = ['diskutil erasevolume HFS+ "' mountPoint ... % String to execute with system which mounts a mac drive. 2048 has to be multiplied to mb to get ram blocks (according to... the internet. bb2021.12.03). 
            '" `hdiutil attach -nomount ram://' ...
            num2str(round(2048*options.ramMb)) '`'];
        [status,cmdout] = system(sysStr); 
        ramDir = options.macRam; 
        if status==0; 
            fprintf('\nMac ram drive %s mounted (I think)\n', options.macRam) 
        else; 
            error('bb2021.12.03: Mac ram drive did not mount. system status %s', num2str(status))
        end
    end
    
end
   
end