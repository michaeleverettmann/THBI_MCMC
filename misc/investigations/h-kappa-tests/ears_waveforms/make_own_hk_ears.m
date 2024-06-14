addpath('/Users/brennanbrunsvik/MATLAB/seis_tools-master/'); % For HK stack
weights = [0.33, 0.33, 0.33]; 
% weights = [0.7, 0.2, 0.1]; 


rfArrPath = 'Ears/gauss_2.5/US.CEH/rfArr.mat';
rfArr = load(rfArrPath);

rf = rfArr.rf;
tt = rfArr.tt'; 
rayP = rfArr.rayParmSecDeg; 

% Load EARS HK stacks from  an example station
rfOrigPath = '/Users/brennanbrunsvik/Documents/repositories/data/models_seismic/US_EARS/EARS_HKStack_US_CEH.mat'; 
rfOrig = load(rfOrigPath); 
hkOrig = rfOrig.hkstack; 
HK_Kor = hkOrig.K; 
HK_Hor = hkOrig.H; 
HK_Aor = hkOrig.Esum; 


%% Look at the waveforms. 
figure(1); clf; hold on; set(gcf, 'color', 'white'); box on; 
title('RF waveforms'); 
xlabel('Time (tt)'); 
ylabel('rf'); 
ylim([-1, 1])
plot(tt, rf)

rfStack = sum(rf, 2); 
plot(tt, rfStack, 'k'); 

%% Make an HK stack using only that waveform stack, where they all have different ray parameters. 
[HK_A,HK_H,HK_K] = HKstack(rfStack, tt, 0.0553, weights, 3.762,...
    linspace(25, 70, 200)',linspace(1.6, 2.1, 201));

% The rest of the script expects transposed hk stack. For IRIS ears data loaded from Zach's functions, K is first dimension, H is second bb2021.12.31
HK_A = HK_A'; 
HK_H = HK_H'; 
HK_K = HK_K'; 

figure(199); clf; hold on; set(gcf,'color','white');
subplot(1,1,1); hold on; 
xlabel('kappa'); ylabel('H'); title('Hk stack from stacked waveforms with different ray parameters'); 
set(gca, 'ydir', 'reverse');         
sf = pcolor(HK_K, HK_H, HK_A'); %Can't use surf. It's 3d. It always covers other plotting objects. 
sf.EdgeAlpha = 0; 
colorbar(); 
xlim([min(HK_K), max(HK_K)]); 
ylim([min(HK_H), max(HK_H)]); 

% Simple wasy to estimate maximum h-k energy point. 
Emax = max(max(HK_A));
[xmax,ymax] = find(Emax==HK_A); 
kBest = HK_K(xmax); 
hBest = HK_H(ymax); 

% Plot position of max energy. 
scatter(gca, kBest, hBest, 50, 'red')
text(kBest, hBest, sprintf('H = %2.2f, k = %1.3f',...
    HK_H(ymax),HK_K(xmax)), 'verticalalignment', 'bottom' )





% % % % figure(198); clf; hold on; set(gcf,'color','white');
% % % subplot(1,2,2); hold on; 
% % % xlabel('kappa'); ylabel('H'); title('Original HK stack'); 
% % % set(gca, 'ydir', 'reverse');         
% % % sf = pcolor(HK_Kor, HK_Hor, HK_Aor'); %Can't use surf. It's 3d. It always covers other plotting objects. 
% % % sf.EdgeAlpha = 0; 
% % % colorbar(); 
% % % xlim([min(HK_K), max(HK_K)]); 
% % % ylim([min(HK_H), max(HK_H)]); 
% % % 
% % % % Simple wasy to estimate maximum h-k energy point. 
% % % Emax = max(max(HK_Aor));
% % % [xmax,ymax] = find(Emax==HK_Aor); 
% % % kBest = HK_Kor(xmax); 
% % % hBest = HK_Hor(ymax); 
% % % 
% % % % Plot position of max energy. 
% % % scatter(gca, kBest, hBest, 50, 'red')
% % % text(kBest, hBest, sprintf('H = %2.2f, k = %1.3f',...
% % %     HK_H(ymax),HK_K(xmax)), 'verticalalignment', 'bottom' )



%% RF stack for each individual trace. Can supply true ray parameters

HK_A_stack = zeros(size(HK_A)); 

figure(2021); clf; hold on; set(gcf, 'color', 'white'); 
title('HK stack on individual RFs'); 
for itrace=[1:size(rf,2)]; 
    subplot(3,3,itrace); hold on; 
    
    % For hk stack: convert from sec per degree to sec per kilometer
    % 0.0553 :: rayP(itrace)/111.1949
    [HK_A,HK_H,HK_K] = HKstack(rf(:,itrace), tt, rayP(itrace)/111.1949,...
        weights, 3.762,...
        linspace(25, 70, 200)',linspace(1.6, 2.1, 201));
    
%     HK_A = HK_A./max(max(HK_A)); 



    % Will have to migrate each to the same ray parameter
    % For now just using same from my synthetic model with 70 degree distances.
    % 

    % The rest of the script expects transposed hk stack. For IRIS ears data loaded from Zach's functions, K is first dimension, H is second bb2021.12.31
    HK_A = HK_A'; 
    HK_H = HK_H'; 
    HK_K = HK_K'; 

    xlabel('kappa'); ylabel('H');% title('Hk stack from (unmigrated) raw receiver functions'); 
    set(gca, 'ydir', 'reverse');         
    sf = pcolor(HK_K, HK_H, HK_A'); %Can't use surf. It's 3d. It always covers other plotting objects. 
    sf.EdgeAlpha = 0; 
    colorbar(); 
    xlim([min(HK_K), max(HK_K)]); 
    ylim([min(HK_H), max(HK_H)]); 

    % Simple wasy to estimate maximum h-k energy point. 
    Emax = max(max(HK_A));
    [xmax,ymax] = find(Emax==HK_A); 
    kBest = HK_K(xmax); 
    hBest = HK_H(ymax); 

    % Plot position of max energy. 
    scatter(gca, kBest, hBest, 50, 'red')
    text(kBest, hBest, sprintf('H = %2.2f, k = %1.3f',...
        HK_H(ymax),HK_K(xmax)), 'verticalalignment', 'bottom' )
    
    HK_A_stack = HK_A_stack + HK_A; 
end
HK_A_stack = HK_A_stack ./ size(rf,2); 

%% Now see what our HK stack looks like if stacking individual HK stacks. 
% HK_A_stack

sgtitle('Original versus redone receiver function'); 

figure(199); clf; hold on; set(gcf,'color','white');
subplot(1,2,1); hold on; 
xlabel('kappa'); ylabel('H'); title('Sum individual HK stacks'); 
set(gca, 'ydir', 'reverse');         
sf = pcolor(HK_K, HK_H, HK_A_stack'); %Can't use surf. It's 3d. It always covers other plotting objects. 
sf.EdgeAlpha = 0; 
colorbar(); 
xlim([min(HK_K), max(HK_K)]); 
ylim([min(HK_H), max(HK_H)]); 

% Simple wasy to estimate maximum h-k energy point. 
Emax = max(max(HK_A_stack));
[xmax,ymax] = find(Emax==HK_A_stack); 
kBest = HK_K(xmax); 
hBest = HK_H(ymax); 

% Plot position of max energy. 
scatter(gca, kBest, hBest, 50, 'red')
text(kBest, hBest, sprintf('H = %2.2f, k = %1.3f',...
    HK_H(ymax),HK_K(xmax)), 'verticalalignment', 'bottom' )

%% Now add original EARS HK stack for comparison
% figure(198); clf; hold on; set(gcf,'color','white');
subplot(1,2,2); hold on; 
xlabel('kappa'); ylabel('H'); title('Original HK stack'); 
set(gca, 'ydir', 'reverse');         
sf = pcolor(HK_Kor, HK_Hor, HK_Aor'); %Can't use surf. It's 3d. It always covers other plotting objects. 
sf.EdgeAlpha = 0; 
colorbar(); 
xlim([min(HK_K), max(HK_K)]); 
ylim([min(HK_H), max(HK_H)]); 

% Simple wasy to estimate maximum h-k energy point. 
Emax = max(max(HK_Aor));
[xmax,ymax] = find(Emax==HK_Aor); 
kBest = HK_Kor(xmax); 
hBest = HK_Hor(ymax); 

% Plot position of max energy. 
scatter(gca, kBest, hBest, 50, 'red')
text(kBest, hBest, sprintf('H = %2.2f, k = %1.3f',...
    HK_H(ymax),HK_K(xmax)), 'verticalalignment', 'bottom' )