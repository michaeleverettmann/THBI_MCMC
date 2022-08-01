%%% In development. Adjust as needed. 
% function plot_receiver_function_sp(); 

tlim = [-60, 60]; 

%%%%
figure(1); clf; hold on; 
set(gcf, 'pos', [-1070 494 651 358], 'color', 'white'); 
h=tiledlayout(3,1,'TileSpacing','compact'); 

nexttile(1); cla; hold on; box on; grid on; xlim(tlim); 
ylabel('Z'); 
plot(tt_sp, predat_sp(:,3))

nexttile(2); cla; hold on; box on; grid on; xlim(tlim); 
ylabel('R'); 
plot(tt_sp, predat_sp(:,1))

nexttile(3); cla; hold on; box on; grid on; xlim(tlim); 
ylabel('T'); 
plot(tt_sp, predat_sp(:,2))

xlabel('Time? (s)'); 


%%%%
figure(2); clf; hold on; 
set(gcf, 'pos', [-1070 494 651 358], 'color', 'white'); 
h=tiledlayout(3,1,'TileSpacing','compact'); 

nexttile(1); cla; hold on; box on; grid on; xlim(tlim); 
ylabel('P'); 
plot(tt_sp, predat_sp_PSV(:,1))

nexttile(2); cla; hold on; box on; grid on; xlim(tlim); 
ylabel('SV'); 
plot(tt_sp, predat_sp_PSV(:,2))

% nexttile(3); cla; hold on; box on; grid on; 
% ylabel('T'); 
% plot(tt_sp, predat_sp(:,2))

xlabel('Time? (s)'); 

% end