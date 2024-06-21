function [SW_V_kernels] = f_mineos_kernels(parm, swperiods)

% parameter_FRECHET_save(parm);  
parameter_FRECHET
% a1_run_mineos_check
a2_mk_kernels

% Get correct frechet variable. 
if ( TYPE == 'S') 
    frech = FRECH_S; 
elseif ( TYPE == 'T')
    frech = FRECH_T; 
end 

% Parameter names 
parmso = ["Z"  , "Vsv", "Vsh", "Vpv", "Vph", "eta", "rho"]; % Parameter names Zach's code
parmsn = ["rad", "vsv", "vsh", "vpv", "vph", "eta", "rho"]; % Parameter names Josh's code
n_per = length(frech); 

% Initialize kernel cell array of structures. 
SW_V_kernels = {}; 
for ip = 1:n_per; 
    SW_V_kernels{ip} = struct();
end

% Convert from Josh's kernel format to Zach's. 
for iv = 1:length(parmso); % For each velocity-like parameter 
    for ip = 1:n_per; % For each period

        % Get input and output parameter names. 
        outp = parmso(iv); % Output parameter, Zach's format. 
        inp  = parmsn(iv); % Input parameter, Josh's format. 

        % Extract input kernel value. If it is not there, it is zeros. 
        if (TYPE=='T') && any(inp == ["vpv", "vph", "eta"]); % No kernel for Love vph, vph
            ker = zeros( size(frech(ip).rad) ); 
        else 
            ker = frech(ip).(inp); % The kernel component we want 
        end

        % Convert between radius and depth.      Some notes: Zach's Z and Josh's rad have the same order (0 to 6371000). But the other kernel values are backwards (radius versus depth). Flip them. Zach's Z: 6370999 is last value. Last Vsv values are 0, top are not. Josh's rad: 6371000 is last value. Last Vsv values are non zero, first are 0. 
        if ~ strcmp(outp, "Z"); 
            ker = flip(ker); 
        end 

        % Assign input kernel to output
        SW_V_kernels{ip}.(outp) = ker; 

        % Insert period. TODO could do only if ip==n_per
        per_out = [frech.per]'; 
        SW_V_kernels{ip}.period = per_out(ip); 
    end
end
SW_V_kernels = transpose(SW_V_kernels);


end