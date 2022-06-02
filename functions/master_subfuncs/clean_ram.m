function [clearStat,clearMessage]=clean_ram(ii,par,options)
    arguments
        ii
        par
        options.clearRamInterval=10
        options.verbose=false
    end

if mod(ii,options.clearRamInterval)==0 &&... % If on iteration where we want to clear the ram. 
        contains(par.res.chainExecFold, 'ramDrive'); % If we know we are in ramDrive
    
    if options.clearRamInterval > 70; 
        [~,diskUsageMes]=system('du -h -d 1'); 
        fprintf(['\nOnly clearning ram every %1.0f Iterations.\n',...
            '  Clearing now. In %s. \n',...
            '  If folder has many files, system() gets very slow. \n',...
            '  Disk usage was: %s\n'], options.clearRamInterval, pwd, diskUsageMes)
    end
    
    % Decide what to delete. brb2022.03.30: Everything not matching *vel_profile or *diary*
    % Plug in par.res.chainExecFold to be sure we don't clear the wrong folder!
    clearRamCmd=sprintf(...
        "shopt -s extglob; find %s -type f -not \\( -name '*vel_profile' -or -name '*diary*' \\) -delete",...
        par.res.chainExecFold); 
    
    [clearStat,clearMessage]=system(clearRamCmd); % Actually clear ram. Scary, in case we are not in the right folder. 
    
    if options.verbose
        fprintf('\nClearing ram. \n%s\n',clearMessage)
    end
    
    if ~clearStat==0; 
        warning('Problem clearing ram. Might just already be empty. Message: \n%s',clearMessage); 
    end
    
elseif ~ contains(par.res.chainExecFold, 'ramDrive'); % You might be trying to clear a NON ram folder. 
    error('Not in executing folder! might clear everything in the NOT ram folder? \n  It does not have "ramDrive" in folder name. \n  Folder was:\n  %s',par.res.chainExecFold); 
else % Not trying to clear
    clearStat = nan; 
    clearMessage = 'Did not clear.'; 
end



end