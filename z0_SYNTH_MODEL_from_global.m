% function [model,laymodel,par] = z0_SYNTH_MODEL_from_global(par,ifplot)

% if nargin < 2 || isempty(ifplot)
%     ifplot=false;
% end

global TRUEmodel TLM

ver = [par.stadeets.sta]'; % Treat the station as determining what model to load. 
% ver = 'cont_Arch'; % Example

% Make velocity model. 
% First make crust. 
% % %% Can run this commented code, after making TRUEmodel with splines, to get crustal values. 
% % % is_crust = find(TRUEmodel.z==45); 
% % % is_crust = is_crust(1); 
% % % vs_c = TRUEmodel.VS(1:is_crust); 
% % % z_c = TRUEmodel.z(1:is_crust); 
% % % vs_c = [3.3000    3.3067    3.3162    3.3284    3.3435    3.3612    3.3818    3.4051    3.4312 3.4600    3.4916    3.5260    3.5617    3.5955    3.6272    3.6567    3.6840    3.7093 3.7324    3.7534    3.7722    3.7889    3.8035    3.8100]'; 
% % % z_c = [0     2     4     6     8    10    12    14    16    18    20    22    24    26    28    30    32    34 36    38    40    42    44    45]'; 
% % % 
% % % 
% % % % Load mantle model. 
% % % SEMum2_avg = SEMum2_avgprofiles(); 
% % % z0  = SEMum2_avg.Z; 
% % % vs0 = SEMum2_avg.(ver);
% % % 
% % % % "continuous" model. 
% % % z = [0:1:300]'; 
% % % z = unique(sort([z; z_c(end)])); % Make sure Moho depth is in mantle model. 
% % % z = z(z>=z_c(end)); % Just keep mantle depths and vals from SEM model. 
% % % vs = interp1(z0, vs0, z, 'spline', 'extrap');
% % % 
% % % % Combine crust and mantle. 
% % % z = [z_c; z]; 
% % % vs = [vs_c; vs]; 
% % %     
% % % % % Quick plot to make sure interpolation and extrapolation aren't messing anything up too much. 
% % % figure(1); clf; hold on; set(gca, 'ydir', 'reverse'); 
% % % plot(vs, z); 
% % % plot(vs0, z0); 

vs_c = [3.3000    3.3067    3.3162    3.3284    3.3435    3.3612    3.3818    3.4051    3.4312 3.4600    3.4916    3.5260    3.5617    3.5955    3.6272    3.6567    3.6840    3.7093 3.7324    3.7534    3.7722    3.7889    3.8035    3.8100]'; 
z_c = [0     2     4     6     8    10    12    14    16    18    20    22    24    26    28    30    32    34 36    38    40    42    44    45]'; 


% Load mantle model. 
SEMum2_avg = SEMum2_avgprofiles(); 
z0  = SEMum2_avg.Z; 
vs0 = SEMum2_avg.(ver);

% "continuous" model. 
z_m = [0:par.mod.dz:300]'; 
z_m = unique(sort([z; z_c(end)])); % Make sure Moho depth is in mantle model. 
z_m = z(z>=z_c(end)); % Just keep mantle depths and vals from SEM model. 
vs_m = interp1(z0, vs0, z_m, 'spline', 'extrap');

% % Combine crust and mantle. 
% z = [z_c; z]; 
% vs = [vs_c; vs]; 
    
% % % Quick plot to make sure interpolation and extrapolation aren't messing anything up too much. 
% figure(1); clf; hold on; set(gca, 'ydir', 'reverse'); 
% plot(vs, z); 
% plot(vs0, z0); 

%% MAKE ALL PARAMETER STRUCTURES
sed = struct('h',0,'VS',[3.3 3.3]);
crust = struct('h',max(vs_c),'vpvs',1.75,'xi',1.05);
mantle = struct('xi',1);
model = struct('sedmparm',sed,'crustmparm',crust,'mantmparm',mantle,...
               'M', nan, 'datahparm', nan, 'selev',0);
           
%% TURN PARMS INTO REAL TARGET MODEL
TRUEmodel = make_mod_from_parms(model,par,'use_splines',false,...
    'vs', struct('crust', vs_c, 'mantle', vs_m),...
    'z' , struct('crust', z_c, 'mantle', z_m));