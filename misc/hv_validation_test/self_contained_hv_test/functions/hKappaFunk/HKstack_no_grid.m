function [Amap] = HKstack(...
    RF,tt,rayp,phase_wts,Vs_av,H_vec,K_vec,options)
% Amap = HKstack(RF,tt,Vs_av,H_grid,K_grid)
    arguments
        RF
        tt
        rayp
        phase_wts
        Vs_av
        H_vec = []
        K_vec = []
        options.kChoice = nan; 
        options.hChoice = nan; 
        options.plotPicks = false; 
        options.vpvsANIS_vpvs = []; 
        options.Vsv_av = []; 
    end
    
    

if ~isempty(options.vpvsANIS_vpvs); %!%!
    % Scale K and v to be anisotropic. 
    % Scale them back to isotropic values at the end. 
    K_vec = K_vec .* options.vpvsANIS_vpvs; % Effective K
    Vs_av_iso = Vs_av; 
    Vs_av = options.Vsv_av; % EFfective v
end


% compute predicted arr
Vp_av   = K_vec.*Vs_av;
kks     = sqrt(Vs_av.^-2 - rayp.^2);
kkp     = sqrt(Vp_av.^-2 - rayp.^2);
% using equations 2-4 from Zhu and Kanamori 2000
t_ps    = H_vec*(kks - kkp); % using outer product
t_ppps  = H_vec*(kks + kkp); % using outer product
t_ppss  = H_vec*(2*kks); % using outer product


% First attempt. Using loop because for each different time t index, we are
% sampling only of of the same index RF. Can it be vectorized? 
phase_wts = phase_wts/sum(phase_wts);
Amap = zeros(1,length(t_ps)); 
method = 'linear'; 
for iTrace = [1:length(t_ps)]; 
Amap(iTrace) =...
    phase_wts(1).*interp1(tt,RF(:,iTrace),t_ps  (iTrace), method) ...
  + phase_wts(2).*interp1(tt,RF(:,iTrace),t_ppps(iTrace), method) ...
  - phase_wts(3).*interp1(tt,RF(:,iTrace),t_ppss(iTrace), method,0); % negative phase!
end

Amap = sum(Amap); 


% % % sum weighted contributions from each phase type
% % % TODO likely need loop
% % phase_wts = phase_wts/sum(phase_wts);
% % Amap =  phase_wts(1).*interp1(tt,RF,t_ps) ...
% %       + phase_wts(2).*interp1(tt,RF,t_ppps) ...
% %       - phase_wts(3).*interp1(tt,RF,t_ppss,[],0); % negative phase!

%   
% if ~isempty(options.vpvsANIS_vpvs); %!%!
%     % Scale K and v back to isotropic values
%     K_vec = K_vec ./ options.vpvsANIS_vpvs;
%     Vs_av = Vs_av_iso; 
% end

% 
% if options.plotPicks; 
%     
%     % temp
%     % Dimensions are H size by K size. 
% %     kChoice = 1.8; % 
% %     hChoice = 45; % 
%     kChoice = options.kChoice; 
%     hChoice = options.hChoice; 
%         
%     t_ps_best   = interpn(H_vec, K_vec, t_ps  , hChoice, kChoice); 
%     t_ppps_best = interpn(H_vec, K_vec, t_ppps, hChoice, kChoice); 
%     t_ppss_best = interpn(H_vec, K_vec, t_ppss, hChoice, kChoice); 
% 
%     
%     figure(132); clf; hold on; 
%     plot(tt, RF); 
%     xlabel('t (s)'); 
%     ylabel('Amp'); 
%     grid on; 
%     box on; 
% 
%     plt_ps   = plot(ones(2,1) .* t_ps_best  , [0, max(RF)*1.2]); 
%     plt_ppps = plot(ones(2,1) .* t_ppps_best, [0, max(RF)*1.2]); 
%     plt_ppss = plot(ones(2,1) .* t_ppss_best, [0,-max(RF)*1.2]); 
%     
%     legend([plt_ps, plt_ppps, plt_ppss], ...
%         'ps', 'ppps', 'ppss');  
% 
%     title(sprintf('H = %2.1f, Vp/Vs = %1.2f', ...
%         hChoice, kChoice)); 
% end

end