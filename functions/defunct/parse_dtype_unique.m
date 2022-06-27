function [ parsed_dtype ] = parse_dtype_unique( datatypes )

dtype_str = string(datatypes)'; 
% for idtype = 1:length(datatypes); 
%     this_dtype = datatypes{idtype}
%     contains(this_dtype, 'SW_HV'); 
% end
any(contains(dtype_str, 'SW_HV')); 

% 
% dlm = [0,strfind(dtype,'_'),length(dtype)+1]';
% pdt = cell(1,length(dlm)-1);
% for iparse = 1:length(pdt)
%     pdt{iparse} = dtype(dlm(iparse)+1:dlm(iparse+1)-1);
% end


if strcmp(pdt{1},'BW')
    % defaults
    parsed_dtype = cell(1,4); 
    parsed_dtype{1} = 'BW';
    parsed_dtype{2} = '';
    parsed_dtype{3} = 'def';
    parsed_dtype{4} = 'def';
    % PHASE
    if any(strcmp(pdt,'Ps'))
        parsed_dtype{2} = 'Ps';
    elseif any(strcmp(pdt,'Sp'))
        parsed_dtype{2} = 'Sp'; 
    end
    % TWIN
    if any(strcmp(pdt,'cms'))
        parsed_dtype{3} = 'cms';
    end        
    % FILTF
    if any(strcmp(pdt,'lo'))
        parsed_dtype{4} = 'lo';
    end
    
elseif strcmp(pdt{1},'RF')
    % defaults
    parsed_dtype = cell(1,4); 
    parsed_dtype{1} = 'RF';
    parsed_dtype{2} = '';
    parsed_dtype{3} = 'def';
    parsed_dtype{4} = 'def';
    % PHASE
    if any(strcmp(pdt,'Ps'))
        parsed_dtype{2} = 'Ps';
    elseif any(strcmp(pdt,'Sp'))
        parsed_dtype{2} = 'Sp'; 
    end
    % TWIN
    if any(strcmpi(pdt,'ccp'))
        parsed_dtype{3} = 'ccp';
    end        
    % FILTF
    if any(strcmpi(pdt,'lo'))
        parsed_dtype{4} = 'lo';
    end
    
elseif strcmp(pdt{1},'HKstack')
    % defaults
    parsed_dtype = cell(1,4); 
    parsed_dtype{1} = 'HKstack';
    parsed_dtype{2} = 'P';
    % PHASE
    if any(strcmp(pdt,'_P'))
        parsed_dtype{2} = 'P';
    elseif any(strcmp(pdt,'_S'))
        parsed_dtype{2} = 'S'; 
    end

elseif strcmp(pdt{1},'SW')
    % defaults
    parsed_dtype = cell(1,4); 
    parsed_dtype{1} = 'SW';
    parsed_dtype{2} = 'Ray';
    parsed_dtype{3} = 'phV';
    % WAVETYPE
    if any(strcmp(pdt,'Ray'))
        parsed_dtype{2} = 'Ray';
    elseif any(strcmp(pdt,'Lov'))
        parsed_dtype{2} = 'Lov'; 
    elseif any(strcmp(pdt,'HV'))
        parsed_dtype{2} = 'HV'; 
        parsed_dtype{3} = 'HVr';
    end
    % VELTYPE
    if any(strcmp(pdt,'ph'))
        parsed_dtype{3} = 'phV';
    elseif any(strcmp(pdt,'gr'))
        parsed_dtype{3} = 'grV'; 
    end
    % Modifications
    if any(strcmpi(pdt,'eks')); % brb2022.06.21. Add in author or other modifications... not a great system, since strings in different positions could match strings in this position. 
        parsed_dtype{4} = 'eks'; 
    elseif any(strcmpi(pdt,'dal')); 
        parsed_dtype{4} = 'dal'; 
    end
end






