function writeHVkernel_execfile( execfile,modelfile,ofile,logfile)
% writeHVkernel_execfile( execfile,modelfile,ofile,logfile)
%   
% Function to write execution file for HVkernel code
% 
% INPUTS:
%  execfile  - name of execution file to write
%  modelfile - name of file with model description
%  ofile     - name of output file with all phV, HVratio, kernel
%  logfile   - name of file to print screen output to



if exist(execfile,'file')==2 %TODOEXIST bb2021.11.22 exist is SUPER slow
    delete(execfile); % kill if it is there 
end

paths = getPaths(); 

%% write synth.in parameter file
fid = fopen(execfile,'w');
fprintf(fid,'#!/bin/csh\n');
%
fprintf(fid,'#\n');
%
fprintf(fid,'echo "Calculating HV ratios and kernels"\n');
%
fprintf(fid,'#\n');
%

% fprintf(fid,'set xdir=/Users/zeilon/Work/codes/HV_Tanimoto/bin\n'); % TODOpath
fprintf(fid,['set xdir=' paths.HV_ellipticity '/bin\n']); % TODOpath
fprintf(fid,'$xdir/HVkernel << ! > %s\n',logfile);
% disp('Temporarily commented writeHVkernel_execfile.m'); 
% disp('Temporarily commented writeHVkernel_execfile.m');
fprintf(fid,'%s\n',modelfile);
fprintf(fid,'%s\n',ofile);
fprintf(fid,'!\n');
% fprintf(fid, 'cp %s %s', logfile, ofile); disp('Adding nonsense to writeHVkernel_execfile.m'); 

fclose(fid);

end




