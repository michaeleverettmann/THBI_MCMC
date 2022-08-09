function [errorInfo]=assess_timeout_results(status, cmdout); 
% Function to see results of [status, cmdout]=system('timeout tlim do_stuff.exec')
% Status of 0 means no known errors (probably depends on do_stuff.exec). 
% Status of 124 means timeout limit reached. Very important to know if this occurs!
% Other status values are... something else. Maybe a do_stuff.exec problem. 
% brb2022.04.12

if (status == 124) % 124 means timeout limit reached. 
    warning(sprintf(['\nError with a system(timeout) call. Timeout limit reached.',...
             '  This is often a cause of a computer running highly parallel, or moving to a new computer.',...
             '  Try increasing timeout duration just before this line. \ncmdout was: \n%s \nbrb2022.04.12'],cmdout)); 
elseif (status ~= 124) && (status ~= 0)% Don't know what this means. 
    warning('Error with a system(timeout) call. Status=%1.0f. Error occured BEFORE timeout limit. Message shown below. brb2022.04.12',status); 
    fprintf('\n%s\n',cmdout)
end

errorInfo = struct(); % For future, could put info in here if you get something besides 124 or 0. 

end