function [HVr,HVK,phV,grV] = readHVkernel_ofile(ofile,swperiods,ifplot)
% [HVr,HVK,phV,grV] = readHVkernel_ofile(ofile,swperiods,[ifplot=false])
%  
%  Function to read KVkernel code output file

narginchk(2,3); 
if nargin < 3 || isempty(ifplot)
    ifplot = false;
end
Np = length(swperiods);

%% open file
fid = fopen(ofile,'r');

%% skip model recapitulation
line=fgetl(fid); % read line with # of lines of model
nskip=str2double(line); % will skip this many lines
for ii = 1:nskip
   fgets(fid);
end

%% how many freqs
a = textscan(fid,'%n frequencies');
Nf = a{1};
doperiods = nan(Nf,1);

%% Set up results
HVr_fs = zeros(Nf,1);
HVK_fs = struct('Z1',[],'Z2',[],...
              'Kph_Vs',[],'Kph_Vp',[],'Kph_rho',[],...
              'Kzh_Vs',[],'Kzh_Vp',[],'Kzh_rho',[], 'Kzh_d',[]);
phV_fs = zeros(Nf,1);
grV_fs = zeros(Nf,1);
fE_fs = zeros(Nf,1);


%% read for each freq
ik = 0;
while ~feof(fid)
    ik = ik+1;
    % grab header for this period
    b = textscan(fid,'frq= %n %n %n %n %n'); % f,error,grV,phV,ZH
    doperiods(ik) = 1./b{1};
    fE_fs(ik) = b{2};
    grV_fs(ik) = b{3};
    phV_fs(ik) = b{4};
    HVr_fs(ik) = b{5};
    % how many lines coming?
    line=fgetl(fid); % read line with # of lines of model
    nmod=str2double(line);
    c = textscan(fid,'%n %n %n %n %n %n %n %n %n',nmod); % ztop,zbot,Kph_rho,Kph_vp,Kph_vs,Kzh_rho,Kzh_vp,Kzh,vp
    HVK_fs(ik,1).Z1 = c{1};
    HVK_fs(ik,1).Z2 = c{2};
    HVK_fs(ik,1).Kph_rho = c{3};
    HVK_fs(ik,1).Kph_Vp = c{4};
    HVK_fs(ik,1).Kph_Vs = c{5};
    HVK_fs(ik,1).Kzh_rho = c{6};
    HVK_fs(ik,1).Kzh_Vp = c{7};
    HVK_fs(ik,1).Kzh_Vs = c{8};
    HVK_fs(ik,1).Kzh_d  = c{9}; 
    
    Kzh_d = HVK_fs(ik,1).Kzh_d; 
    Kzh_d(2:end-1) = Kzh_d(3:end)-Kzh_d(1:end-2); 
    Kzh_d(1) = Kzh_d(2); % Not sure what to do about top and bottom of model. 
    Kzh_d(end) = Kzh_d(end-1); 
    HVK_fs(ik,1).Kzh_d = Kzh_d; 
    
%     Kzh_d_diff = diff(Kzh_d); 
%     Kzh_d = max([[Kzh_d_diff;0], [0;Kzh_d_diff]], [], 2); % Don't know whether to associate with each point the diff to the point above or below. So just give the larger of the two values. 
%     HVK_fs(ik,1).Kzh_d = Kzh_d; 
    
    % account for possibilty of error
    fpos = ftell(fid);
    fgets(fid);
    if feof(fid), break; end
    line1=fgetl(fid); %
    line1parse = textscan(line1,'%s');
    ifer = strcmp(line1parse{1}{1}(1),'?');
    % rewind one if no error line
    if ~ifer, fseek(fid,fpos,-1); fgets(fid); end
    
    field_names = fieldnames(HVK_fs); 
    for ifn = 1:length(field_names); 
        if any(isnan([HVK_fs.(field_names{ifn})])); 
            error('Nan values in HV kernel! Not usuable. Is this fixable? brb2022.07.01'); % brb2022.07.01 I'm finding nan values at some indicies. This is sneaky and can break the rest of the code: An integral over kenel with any nan will be nan, so this kernel will be completely uselss as is. 
        end
    end
    
end
fclose(fid);

doperiods = round_level(doperiods,0.05);
rswperiods = round_level(swperiods,0.05);
if isequal(doperiods,rswperiods)
    HVr = HVr_fs;
    phV = phV_fs;
    grV = grV_fs;
    HVK = cell(Np,1); 
    fns = fieldnames(HVK_fs);
    for ip = 1:Np
        for ifn = 1:length(fns)
        HVK{ip}.(fns{ifn}) = HVK_fs(ip).(fns{ifn});
        end
    HVK{ip}.period = swperiods(ip);
    end

else
    
    %% interpolate to get kernels and velocities
    imeth = 'spline'; 
    HVr = interp1(doperiods,HVr_fs,swperiods,imeth);
    phV = interp1(doperiods,phV_fs,swperiods,imeth);
    grV = interp1(doperiods,grV_fs,swperiods,imeth);
    % interp kernels
    ind = interp1(doperiods,1:Nf,swperiods,imeth);
    fns = fieldnames(HVK_fs);
    HVK = cell(Np,1);
    for ip = 1:Np
        i1 = floor(ind(ip));
        i2 = ceil(ind(ip));
        f2 = rem(ind(ip),1);
        f1 = 1-f2;

        for ifn = 1:length(fns)
            HVK{ip}.(fns{ifn}) = f1*HVK_fs(i1).(fns{ifn}) + f2*HVK_fs(i2).(fns{ifn});
        end

        HVK{ip}.period = swperiods(ip);
    end
end

ifplot = 1; if(ifplot); warning('Forcing ifplot = 1'); end
if ifplot
    figure(3001);clf, 
    fns = {'Kzh_Vp','Kzh_Vs','Kzh_rho','Kzh_d'};
    for ifn = 1:length(fns)
        subplot(1,length(fns),ifn),hold on; box on; 
        for ip = 1:ik
            plot(HVK_fs(ip).(fns{ifn}),HVK_fs(ip).Z1,':');
        end
        for ip = 1:length(swperiods)
            plot(HVK{ip}.(fns{ifn}),HVK{ip}.Z1,'linewidth',2);
        end
        set(gca,'ydir','reverse','fontsize',15)
        title(regexprep(fns{ifn},'_','-'))

% This commented code changes limits to compare with Tanimoto 2009 GJI
%         ylim([0,20]); warning('Fake y lim in readHVkernel_ofile.m plotting'); 
%         if strcmp(fns{ifn},'Kzh_d'); 
%             xlim([-3e-3, 3e-3]); 
%         else
%             xlim([-.5, .5]); 
%         end
    end
end
