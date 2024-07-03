%%% Merge net sta combinations from several station.txt files. 
%%% Removes duplicates. 
%%% brb2022.09.26

sta_lists = {'staList_all.txt' ...
    'stations_ne_cirlce2.txt' 'stations_se_circl2.txt'}; % Manually choose files to include. 

sta_list_new = ""; 
for ilst = 1:length(sta_lists); % Loop over each file and combine their contents. 

    fsta_lst_i = sta_lists{ilst} ; 

    sta_lst_i = string(table2cell(readtable(...
        fsta_lst_i, 'ReadVariableNames', 0))); 
    sta_lst_i = sta_lst_i(:,1) + ' ' + sta_lst_i(:,2); % Combine net and sta
    sta_lst_i = sta_lst_i(1:end-1,:); % Remove non station final line

    sta_list_new = [sta_list_new; sta_lst_i]; 
end; 

% Minor formatting fixes. 
sta_list_new = unique(sta_list_new); 
sta_list_new = sta_list_new(sta_list_new.strlength>0); 

% Write a version for using with IRIS ears retrieval code. 
sta_list_ears = sta_list_new; 
sta_list_ears = sta_list_ears.replace(' ', ','); 
writecell(sta_list_ears.cellstr(), 'stationList.csv',...
    "QuoteStrings", 'none')

% Write a version to be compatible with the MCMC slurm submission
sta_list_new = [sta_list_new; ...
    "# Extra line to make sure bash script reads properly."]; 
writecell(sta_list_new.cellstr(), 'sta_list_merged.txt')
