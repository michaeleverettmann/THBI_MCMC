%%% In development. Adjust as needed. 

figure(1); clf; hold on; 
h=tiledlayout(2,2,'TileSpacing','compact'); 

nexttile(1); cla; hold on; box on; grid on; 
title('ZRT', 'fontweight', 'normal'); 
plot(tt_sp, predat_sp_ZRT); 

nexttile(2); cla; hold on; box on; grid on; 
title('PSV', 'fontweight', 'normal'); 
plot(predat_sp_PSV); 