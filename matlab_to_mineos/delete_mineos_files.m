function delete_mineos_files( ID,R_or_L )
%delete_mineos_files( ID )
%   Function to delete all files associted with mineos running

%% filenames
if ~ischar(ID), ID = num2str(ID);end
cardfile = [ID,'.model'];       %if exist(cardfile,'file')~=2, cardfile = ''; end
ID = [ID,R_or_L];
execfile = [ID,'.run_mineos'];  %if exist(execfile,'file')~=2, execfile = ''; end
eigfile = [ID,'*.eig'];          %if exist(eigfile,'file')~=2, eigfile = ''; end
eigfile_fix = [ID,'*.eig_fix'];          %if exist(eigfile,'file')~=2, eigfile = ''; end
qfile = [ID,'.q'];
logfile = [ID,'.log'];          %if exist(logfile,'file')~=2, logfile = ''; end
execfile_k = [ID,'.run_kernels'];%if exist(execfile_k,'file')~=2, execfile_k = ''; end
stripfile = [ID,'.strip'];      %if exist(stripfile,'file')~=2, stripfile = ''; end
tabfile = [ID,'.table'];        %if exist(tabfile,'file')~=2, tabfile = ''; end
tabfile_hdr = [tabfile,'_hdr'];        
branchfile = [tabfile,'_hdr.branch'];   
kernelfile = [ID,'.cvfrechet']; %if exist(kernelfile,'file')~=2, kernelfile = ''; end

% preamble
wd = pwd;

delete([ID,'_*.asc'])
delete([ID,'_*.eig'])
delete([ID,'_*.eig_fix'])

if java.io.File([pwd '/' execfile],'file').exists, delete(execfile); end
if java.io.File([pwd '/' cardfile],'file').exists, delete(cardfile); end
if java.io.File([pwd '/' eigfile],'file').exists, delete(eigfile); end
if java.io.File([pwd '/' eigfile_fix],'file').exists, delete(eigfile); end
if java.io.File([pwd '/' qfile],'file').exists, delete(qfile); end
if java.io.File([pwd '/' stripfile],'file').exists, delete(stripfile); end
if java.io.File([pwd '/' tabfile],'file').exists, delete(tabfile); end
if java.io.File([pwd '/' kernelfile],'file').exists, delete(kernelfile); end
if java.io.File([pwd '/' tabfile],'file').exists, delete(tabfile); end
if java.io.File([pwd '/' tabfile_hdr],'file').exists, delete(tabfile_hdr); end
if java.io.File([pwd '/' branchfile],'file').exists, delete(branchfile); end
% % % Alternate form of above :: I think these files are generally not present. The java file exists is so fast, and system call is fairly slow. I think it will save time to do the java test followed by delete... I left an alternative code below anyway. bb2021.11.24
% % system(['rm -f ' pwd '/' execfile      ]); % Using -f means will delete file if it exists and not return error if it doesn't exist. Easier to read, maybe faster. 
% % system(['rm -f ' pwd '/' cardfile      ]);
% % system(['rm -f ' pwd '/' eigfile       ]);
% % system(['rm -f ' pwd '/' eigfile_fix   ]);
% % system(['rm -f ' pwd '/' qfile         ]);
% % system(['rm -f ' pwd '/' stripfile     ]);
% % system(['rm -f ' pwd '/' tabfile       ]);
% % system(['rm -f ' pwd '/' kernelfile    ]);
% % % system(['rm -f ' pwd '/' tabfile       ]); % This was present twice...
% % system(['rm -f ' pwd '/' tabfile_hdr   ]);
% % system(['rm -f ' pwd '/' branchfile    ]);

if java.io.File([pwd '/' execfile_k]).exists %TODOEXIST bb2021.11.22 exist is SUPER slow
%     try
% %     delete(execfile_k,stripfile,tabfile,kernelfile,[tabfile,'_hdr'],[tabfile,'_hdr.branch']);
%     catch
%         fprintf('Tried to delete MINEOS files but some error');
%     end
    for ip = 1:length(swperiods)
        delete(ikernelfiles{ip}); 
    end
end
% % system(['rm -f ' pwd '/' execfile_k    ]);
% % system(['rm -f ' pwd '/' logfile       ]);
if java.io.File([pwd '/' execfile_k]).exists, delete(execfile_k); end 
if java.io.File([pwd '/' logfile]).exists, delete(logfile); end 

% postamble
cd(wd);

end

