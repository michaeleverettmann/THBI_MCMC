% Goals: (1) Check if we have results for each desired station. (2) Copy
% results from different computers to a common area (if they haven't been
% coppied). 
%
% Assumes the desired stamp is "standard" (for now). 

clc; clear;
fprintf('Running examine_results from %s\n', pwd); 

f_sta_list_all = '/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/ENAM/batch/staList_all.txt'; 
comps = {'STASinv_pod', 'STASinv_eri', 'STASinv_cnsi'}; % All computers that might have good results. 
comp_dest = 'STASinv_collate'; % Where we want to copy results to, from individual computers. 

% String of folder that stores results. Use with sprintf. Give it computer folder, station, the network.  
sta_dir = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/%s/%s_%s_dat1/standard/'; 
res_str = [sta_dir 'allmodels_perchain_orig.mat']; 
copy_res_str = '/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/%s/%s_%s_dat1/standard'; 

desired_chains = 12; 
desired_iter   = 16000; 

%% Read all stations we want to have results for. 
sta_list_all = string(table2cell(readtable(...
    f_sta_list_all, 'ReadVariableNames', 0))); 
sta_list_all = sta_list_all(:,1) + ' ' + sta_list_all(:,2); % Combine net and sta
sta_list_all = sta_list_all(1:end-1,:); % Remove non station final line

have_res = logical(zeros(length(sta_list_all), length(comps))); % Build up array saying where we have results. 
res_qual = logical(zeros(length(sta_list_all), length(comps))); % How much we might trust results from that computer for that station. 

%% Figure out which stations we have results for, and on which computer. 
% FYI Due to glitches with computers, I had to switch from one computer to another after part of the inversion was done. 
for istn = 1:length(sta_list_all); 
for icomp = 1:length(comps); 
    comp = comps{icomp}; 
    sta_str = sta_list_all(istn); 
    sta_str_splt = split(sta_str, ' '); 
    net = sta_str_splt(1); 
    sta = sta_str_splt(2); 


    file_res = sprintf(res_str, comp, sta, net); 
    fpar = sprintf([sta_dir '%s'], comp, sta, net, 'par.mat'); 
    fmod = sprintf([sta_dir '%s'], comp, sta, net, 'final_model.mat'); 
    fmis = sprintf([sta_dir '%s'], comp, sta, net, 'final_misfit.mat'); 

    have_files = exist(file_res, 'file') * ...
        exist(fpar, 'file') * exist(fmod, 'file') * exist(fmis, 'file'); 

    if ~ have_files; continue; end

    % Make sure these are the results we wanted. 
    par=load(fpar); par=par.par; 
    have_files = have_files *...
        ((par.inv.niter==desired_iter) && (par.inv.nchains==desired_chains)); % Don't plot results if we didn't run the inversion with enough chains or iterations. Might have been test runs..
    
    if have_files; 
%         fprintf('%s %s exists in %s.\n', net, sta, comp); 
        have_res(istn, icomp) = 1; 
        
        %%% Now copy results to generic folder. 
        res_to   = sprintf(copy_res_str, comp_dest, sta, net); 
        res_from = sprintf(copy_res_str, comp     , sta, net); 
        if ~ exist([res_to '/allmodels_perchain_orig.mat']); % If we haven't already copied results. 
            fprintf('Copying from \n   %s\n   to\n   %s\n',res_from, res_to); 
            cp_cmd = sprintf('mkdir -p %s; cp -a %s/* %s', res_to, res_from, res_to); 
            [cp_status, cp_message] = system(cp_cmd); 
            system(sprintf(...
                'echo "From %s. \n\nThese results copied from %s" !> %s/original_path.txt',...
                comp, res_from, res_to)); % Save a file saying where we got these results. 
        end
        %%% END copy results
    end

end
end

have_a_result = (sum(have_res,2) > 0); 

dup_stas = sta_list_all(sum(have_res,2)>1 ); % Warn if have useable results at multiple stations. 
if length(dup_stas) > 0; 
    fprintf('Have at multiple stations: %s\n',  dup_stas); 
end

%% Save list of stations we still need 
have_res_stn = any(have_res,2); % We have results for this station already, from some computer. 
no_res_stn = ~ have_res_stn; 
stn_lst_need = sta_list_all(no_res_stn); % List of all stations where we don't have results. Or, at least where we don't have the allmodels_perchain_orig.mat file. 
stn_lst_need = unique(stn_lst_need); % Remove any duplicates. 

system('rm station_list_still_need.txt'); 
for iline = 1:length(stn_lst_need); 
    system(sprintf('echo "%s" >> station_list_still_need.txt', stn_lst_need(iline))); 
end
system('echo "# Adding strange necessary end line" >> station_list_still_need.txt'); 

