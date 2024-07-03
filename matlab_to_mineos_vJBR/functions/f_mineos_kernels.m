function [SW_V_kernels] = f_mineos_kernels(phV, grV, parm, swperiods)

%% Wrap modified version of Josh's Matlab wrapper to Mineos. 
parameter_FRECHET_save(parm, swperiods); % Always re-save parameters. If we had just ran Rayleigh, and are now going to Love, then we would be loading the wrong parameters. 
parameter_FRECHET
a2_mk_kernels 

%% Get kernels into MCMC format

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

        % Direction. convert between radius and depth. 
        if strcmp(outp, "Z"); 
            ker = 6371000-ker; 
        end
        ker = flip(ker); 

        % Units. Josh's kernels were 1000 times smaller than Zach's. 
        if ~ strcmp(outp, "Z"); 
            ker = ker * 1000; 
        end 

        % Josh's version of Mineos has different units. To undo that, we
        % need to multiply kernels by model_value(z)/phase_velocity(frequency). brb20240702
        if ~strcmp(outp, "Z"); % This block handles units. 

            % Model values. 
            rad_mod = mod_in.('rad'); 
            mod_mult = mod_in.(inp); 
            z_mod = flip(6371000-rad_mod)./1000; 
            mod_mult = flip(mod_mult); 

            % Sensitivities. ker was made above. 
            z_frech = flip(6371000-frech(ip).('rad'))./1000; % Get z. Process consistent with how ker was made. 
            
            % Fix units. eta was unitless, so leave it alone.  
            if ~ strcmp(outp, "eta"); 
                mod_mult = mod_mult./1000; 
            end

            % Phase velocity at this period. Will normalize sensitivity to this. 
            phv_grv_mult = phV(ip); 

            % Remove duplicate points so we can interpolate to the same z basis. Love waves are sometimes cut off at the outer core, so have less z indices. Do a moving average. This returns kernels that are virtually identifical to if there is no interpolateion. 
            terps = [0.999998, 1-0.999998]; % Moving average weights, for interpolating. 
            z_mod(2:end-1) = terps(1)*z_mod(2:end-1) + terps(2)*z_mod(3:end) + terps(2)*z_mod(1:end-2); 
            z_frech(2:end-1) = terps(1)*z_frech(2:end-1) + terps(2)*z_frech(3:end) + terps(2)*z_frech(1:end-2); 

            % Interpolate, to handle cases where number of Zs is different
            mod_mult = interp1(z_mod, mod_mult, z_frech, 'linear'); 
            
            % % Plot kernel and model parameter to make sure their z's are consistent. 
            % figure(1); clf; hold on; 
            % plot(ker*max(mod_mult)/max(ker), z_frech, 'DisplayName', 'Frechet'); 
            % plot(mod_mult, z_mod, 'DisplayName', 'Model'); 
            % set(gca, "YDir", 'reverse'); 
            % % ylim([0, 300]); 
            % legend(); 

            % Make sure the units are consistent between phase velocity and our model parameter. 
            unit_check = max(mod_mult)/max(phv_grv_mult); 
            if (unit_check < 0.1) || (unit_check > 10) % Check that mod_mult and phv_grv_mult have consistent units. Doesn't matter if we use km, m, whatever, because the units should cancel out. Just have to be consistent. brb20240702
                sprintf('phv or grv is factor of %f different from model parameter. They should have consistent units, but seem not to. ',unit_check); % Quick unit check
                disp(outp)
                disp(max(mod_mult))
            end

            % Multiply kernel by model value / phase velocity for units. 
            ker = ker .* mod_mult / phv_grv_mult; % TODO handle grV case. % Here use mod_in
        end 

        % Assign input kernel to output
        SW_V_kernels{ip}.(outp) = ker; 

        % Insert period. TODO could do only if ip==n_per
        per_out = [frech.per]'; 
        SW_V_kernels{ip}.period = per_out(ip); 
    end
end

% Loaded
SW_V_kernels = transpose(SW_V_kernels);

% figure(1); clf; hold on; 
% kt = targ{1}; 
% plot(kt.Vsv, kt.Z, 'DisplayName', 'Example');
% kn = SW_V_kernels{1}; 
% plot(kn.Vsv, kn.Z, 'DisplayName', 'New'); 
% legend(); 

% % Can compare to an example
% SW_V_kernels = targ; 

ifplot = false; 
if ifplot; 
    figure(1008); clf; hold on; 
    for ip = 1:n_per; 
        k = SW_V_kernels{ip}; 
        plot(k.Vsh, k.Z*1e-6); 
    end
end


end