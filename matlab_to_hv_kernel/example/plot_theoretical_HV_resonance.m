b = linspace(0.5,5,100); 
z = linspace(0.1,50,201); 

[z,b] = ndgrid(z,b); 

f = b./(4.*z); 
t = 1./f; 

figure(1); clf; hold on; box on; 
title('Resonant period = 4 * Z / V_s', 'fontweight', 'normal');
xlabel('V_s'); 
ylabel('Z_{layer}'); 
contourf(b, z, t, 40, 'LineColor', 'none'); 
cbar = colorbar(); 
caxis([0.001, 200]); 