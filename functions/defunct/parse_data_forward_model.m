function [ datap ] = parse_data_forward_model( data )

% % Loop over fields in data. Make new data structures that instruct how we
% % should do forward modelling. 
% datap = struct('BW_Ps',[],'BW_Sp',[],...
%                  'RF_Ps',[],'RF_Sp',[],...
%                  'RF_Ps_ccp',[],'RF_Sp_ccp',[],...
%                  'HKstack_P',[],...
%                  'SW_Ray_phV',[],'SW_Lov_phV',[],'SW_HV',[]); % List of anything that might require forward modelling. brb2022.06.24 I'm only thinking through surface waves right now...
% datap_names = string(fieldnames(datap)); 
%              
% 
% 
% 
% fn = fieldnames(data);
% for k=1:numel(fn)
%     field_name = fn{k}; 
%     dat_inf = split(field_name, '_'); 
%     
%     base_dat = [dat_inf{1} '_' dat_inf{2} '_' dat_inf{3}]; 
%     if any(base_dat == datap_names); 
%         datap.(base_dat) = 
%     
%     data.(fn{k})
% end

warning('brb2022.06.27 I dont think this script should be working. It was a test. If you find it in use, see whats happening. Otherwise, delete this script.')

datap = struct(); 

fn = fieldnames(data);
for k=1:numel(fn)
    field_name = fn{k}; 
    dat_inf = split(field_name, '_'); 
    
    if ~isfield(datap, dat_inf{1});
        datap.(dat_inf{1}) = struct(); 
    end
    
%     base_dat = [dat_inf{1} '_' dat_inf{2} '_' dat_inf{3}]; 
%     if any(base_dat == datap_names); 
%         datap.(base_dat) = 
    
%     data.(fn{k})
end

end