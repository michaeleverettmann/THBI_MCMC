function [mismatch] = hKappaError(scale_error, E_by_Emax, min_error, options); 
    arguments
        scale_error
        E_by_Emax
        min_error
        options.HK_Esum = 1; 
    end
%bb2022.01.27
% Function to go from h-kappa energy stack (E_by_Emax) to mismatch.

% E_by_Emax = E_by_Emax ./ (max(max(E_by_Emax))); % Double check that it is normalized. 
% Nevermind, we would need the E ARRAY to do that. 
E_by_Emax = E_by_Emax ./ maxgrid(options.HK_Esum); 

mismatch = scale_error .* (1-E_by_Emax + min_error); 

end