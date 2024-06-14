function writeHVkernel_modfile( ofile,Z,vp,vs,rho,Qp,Qmu,freqs ) %#ok<INUSL>
% writeHVkernel_modfile( ofile,Z,vp,vs,rho,Qp,Qmu,freqs ) %#ok<INUSL>
%   
% This function writes a model file for Toshiro Tanimoto's KV kernel
% script. The function expects that you provide rho,vp,vs, qp, qs at a set
% of crust/upper mantle depths.
% 
% Z. Eilon  07/2018


if isempty(ofile)
    ifwrite=false;
else
    ifwrite=true;
end
if nargin < 8 || isempty(freqs) 
    freqs = 1./[100 50 20 10];
end

  
Re = 6371;
dz = 1; % vertical resolution for output. brb2022.07.15 Doesn't seem to work if less than 1. But that can probably be fixed.  
nfreq = length(freqs);

%% Q values - Grab PREM
global prem_isotropic
if nargin<7 || isempty(Qmu) || nargin<6 || isempty(Qp)
    if isempty(prem_isotropic)
        prem_isotropic = prem;
    end

    prem_mod = prem_isotropic;

    if isempty(Qmu)
        Qmu = linterp(prem_mod.depth,prem_mod.qu,Z); %#ok<NASGU>
    end
    if isempty(Qp)
        Qk = linterp(prem_mod.depth,prem_mod.qk,Z); %#ok<NASGU>
        Qp =qkqu2qp(Qk,Qmu,vp,vs);
    end

end

%% m

%% Write model file
if ifwrite
N = length(Z);
% edit modelfile name so no periods
ofile_print = ofile;
ofile_print(regexp(ofile_print,'\.')) = '_';

fid = fopen(ofile,'w+');
fprintf(fid,'%u\n',nfreq); % frequencies
for ii = 1:nfreq
    fprintf(fid,'%.6f\n',freqs(ii)); % frequencies
end
fprintf(fid,'%.2f %u\n',dz,ceil(max(Z)./dz)); % depths for output (dz, Nz)
fprintf(fid,'%u\n',N);
for kk = 1:N % expect max radius first
    % output:
%      Radius (km), rho (g/cc), Vp (km/s), Vs (km/s), 1/Qp, 1/Qs
    fprintf(fid,'%9.7e  %9.7e  %9.7e  %9.7e  %9.7e  %9.7e\n',...
        Re - Z(kk),rho(kk),vp(kk),vs(kk),1./Qp(kk),1./Qmu(kk));
end
fclose(fid);
end

end

