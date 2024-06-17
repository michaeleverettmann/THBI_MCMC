function axs = plot_TRUvsPRE_WAVEFORMS( trudata,predata,posterior,par,ifsave,ofile,ifnorm)
%plot_TRUvsPRE_WAVEFORMS( trudata,predata,ifsave,ofile )
%   
% function to plot predicted and true seismograms (Vertical and Radial)
% Assumes date in 3-column ZRT matrices with equal sample rate 

if nargin < 5 || isempty(ifsave)
    ifsave=false;
end

if nargin < 6 || isempty(ofile)
    ofile='true_vs_predicted_data';
end

if nargin < 7 || isempty(ifnorm)
    ifnorm=true;
end

boxLineWidth = 1.75; 
titleSize = 20; 
xlims = [-3 31;-31 3]; %[P;S]

dtypes = fieldnames(predata);
dtypes(strcmp(dtypes,'SW')) = [];


figure(58),clf,set(gcf,'pos',[44 150 1500 1100])
ax1  = axes('position',[0.05 0.69 0.20 0.26]); hold on % Ps? 
ax2  = axes('position',[0.05 0.39 0.20 0.26]); hold on % Sp? 
ax3  = axes('position',[0.28 0.69 0.20 0.15 ]); hold on % Sp, HK present I guess. 
ax4  = axes('position',[0.28 0.39 0.20 0.15 ]); hold on % Sp, parent pulse. 
ax5  = axes('position',[0.51 0.69 0.20 0.26]); hold on
ax6  = axes('position',[0.51 0.39 0.20 0.26]); hold on
ax7  = axes('position',[0.72 0.69 0.20 0.26]); hold on
ax8  = axes('position',[0.72 0.39 0.20 0.26]); hold on

ax9  = axes('position',[0.05 0.05 0.24 0.21]); hold on % Rayleigh 
ax10 = axes('position',[0.36 0.05 0.24 0.21]); hold on % Love
ax11 = axes('position',[0.67 0.05 0.24 0.21]); hold on % HV

axs=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11];
for iaxs = [1:length(axs)]; 
    set(axs(iaxs), 'Box', 'on'); 
    set(axs(iaxs), 'linewidth', boxLineWidth); 
end


if any(strcmp(fieldnames(trudata),'HKstack_P'))
    ax12  = axes('position',[0.05 0.39 0.20 0.3]); hold on % Ps? 
    delete([ax1,ax2]);
    set(ax12, 'Box', 'on'); 
    set(ax12, 'linewidth', boxLineWidth); 
    axs=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11];
end

for id = 1:length(dtypes)
dtype = dtypes{id};
pdtyp = parse_dtype( dtype );
switch pdtyp{1}

%% SW
    case 'SW'
        switch pdtyp{2}
            case 'Ray', 
                ax = ax9; 
                ylabstr = 'Phase Velocity (km/s)';
                set(ax, 'xscale', 'log'); 
                xticks(ax, [5, 10, 20, 40, 60, 100, 180]); 
            case 'Lov', 
                ax = ax9; 
                ylabstr = 'Phase Velocity (km/s)';
            case 'HV',  
                ax = ax11;
                ylabstr = 'H/V ratio'; 
                if ~ par.inv.synthTest; 
                    errorbar(ax,trudata.(dtype).periods,trudata.(dtype).HVr,2*trudata.(dtype).sigma.*ones(size(trudata.(dtype).periods)),'k')
                end
                set(ax, 'xscale', 'log'); 
                xticks(ax, [16, 20, 30, 40, 50, 70, 90]); 
        end
    hp(1) = plot(ax,trudata.(dtype).periods,trudata.(dtype).(pdtyp{3}),'ko-','linewidth',2,'markersize',10);
    hp(2) = plot(ax,predata.(dtype).periods,predata.(dtype).(pdtyp{3}),'r-x','linewidth',2,'markersize',10);
    hl = legend(ax,hp,'True','Pred','Location','SouthEast'); set(hl,'fontsize',15);
    set(ax,'fontsize',15)
    xlabel(ax,'Period (s)','fontsize',18)
    ylabel(ax,ylabstr,'fontsize',18)
    title(ax, 'Surface waves','fontsize',titleSize, 'fontweight', 'normal') ; 

%% BWs
    case {'BW','RF'}
        axn = 0;
        switch pdtyp{2}
            case 'Ps', axn = axn+1; 
            case 'Sp', axn = axn+2;
        end
        if strcmp(pdtyp{4},'lo'),axn = axn+2; end



    xa1 = axs(2*axn-1); % order [5,7,1,3]
    xa2 = axs(2*axn); % order [6,8,2,4]
    if strcmp(pdtyp{2}(1),'P'), ps=1;xp=1;xsv=0.1;elseif strcmp(pdtyp{2}(1),'S'), ps=2;xp=0.1;xsv=1; end
    
    if ~isempty(predata.(dtype)) && ~isempty(predata.(dtype)(1).PSV)
        trunrm = zeros(length(trudata.(dtype)),1);
        prenrm = zeros(length(predata.(dtype)),1);
        for itr = 1:length(trudata.(dtype))
            if strcmp(pdtyp{3},'ccp')
                trudata.(dtype)(itr).tt = trudata.(dtype)(itr).zz;
                predata.(dtype)(itr).tt = predata.(dtype)(itr).zz;
            end
          
            
            if ifnorm
                trunrm(itr) = norm(trudata.(dtype)(itr).PSV); % get the norm of the traces, to normalise power
                prenrm(itr) = norm(predata.(dtype)(itr).PSV); % get the norm of the traces, to normalise power
            else
                trunrm(itr) = 1;
                prenrm(itr) = 1;
            end
            plot(xa1,trudata.(dtype)(itr).tt,trudata.(dtype)(itr).PSV(:,1)./trunrm(itr),'k','linewidth',2.5)
            plot(xa1,predata.(dtype)(itr).tt,predata.(dtype)(itr).PSV(:,1)./prenrm(itr),'r','linewidth',1.5)
            plot(xa2,trudata.(dtype)(itr).tt,trudata.(dtype)(itr).PSV(:,2)./trunrm(itr),'k','linewidth',2.5)
            plot(xa2,predata.(dtype)(itr).tt,predata.(dtype)(itr).PSV(:,2)./prenrm(itr),'r','linewidth',1.5)
        end
        set(xa1,'xlim',xlims(ps,:),...
                'ylim',0.3*xp*[-1 1],...
                'fontsize',13)
        set(xa2,'xlim',xlims(ps,:),...
                'ylim',0.5*xsv*[-1 1],...
                'fontsize',13)
        title(xa1, regexprep(dtype,'_','-'),'fontsize',titleSize, 'fontweight', 'normal')
        xlabel(xa2, sprintf('Time from %s arrival',pdtyp{2}(1)),'fontsize',18)
        
        if strcmp(pdtyp{3},'ccp')
            set([xa1,xa2],'xlim',minmax(trudata.(dtype)(itr).zz')); 
            xlabel([xa1,xa2],'Depth of conversion (km)','fontsize',18)
        end
    
    end
    
%% HKstack        
    case {'HKstack'}         
        [~,contH] = contourf(ax12,predata.(dtype).Kgrid,predata.(dtype).Hgrid,...
            predata.(dtype).Esum',30,'linestyle','none'); % Use final_predata to get the HK stack estimated with our velocity model. 

        cbar = colorbar(ax12,'eastoutside'); 
        set(cbar, 'fontsize', 12); 
        try 
            colormap(viridis); 
        catch 
            warning('Missing colormap viridis. Should be in repositories somewhere. '); 
        end
        
        plot(predata.HKstack_P.K,predata.HKstack_P.H,'ok','linewidth',2,...
            'markerfacecolor','r','markersize',7)
        
        title(ax12, regexprep(dtype,'_','-'),'fontsize',titleSize, 'fontweight', 'normal')
        xlabel(ax12, 'Vp/Vs ratio','fontsize',16)
        ylabel(ax12, 'Moho depth','fontsize',16)
        set(ax12,'ydir','reverse', 'linewidth', 4)
  
end
end

pause(0.001)

if ifsave
    exportgraphics(gcf, ofile, 'resolution', 300); 
end

if any(string(par.inv.datatypes) == 'HKstack_P'); 
    dtype = 'HKstack_P'; 
    xlim_arr = [1.6, 2.1];
    ylim_arr = [15 , 70 ];
    fsetlims = @()set(gca(),'xlim', xlim_arr, 'ylim', ylim_arr, ...
        'ydir', 'reverse'); 


    figure(1058); clf; hold on; 
    set(gcf, 'pos', [405 698 320 383]); 
    axPdf = gca(); cla; box on; 
    fsetlims(); 
    ax12 = copyobj(axPdf, gcf); 
    set(ax12, 'Color', 'none'); 
    fsetlims(); 

    %%% Top 
    axes(ax12); 
    Esum = predata.(dtype).Esum'; 
    contour(ax12,predata.(dtype).Kgrid,predata.(dtype).Hgrid,...
        predata.(dtype).Esum', linspace(min(min(Esum)), max(max(Esum)),5), 'k'); % Use final_predata to get the HK stack estimated with our velocity model. 

    plot(predata.HKstack_P.K,predata.HKstack_P.H,'ok','linewidth',2,...
        'markerfacecolor','r','markersize',7)
    %%%

    ax12.Visible = 'on'; 
    axPdf.Visible = 'off'; 

    %%% Bottom
    axes(axPdf); 
    kgrd = predata.HKstack_P.Kgrid; 
    hgrd = predata.HKstack_P.Hgrid; 
    kgrd = linspace(min(kgrd), max(kgrd), 101); 
    hgrd = linspace(min(hgrd), max(hgrd), 100); 
    [hist_counts] = hist3(...
        [posterior.vpvs, posterior.zmoh],'ctrs',{kgrd,hgrd}); 
    hist_counts = imgaussfilt(hist_counts, 'FilterSize', 3); 
    hist_counts = log(hist_counts); 
    [pdf_hand] = contourf(axPdf, kgrd, hgrd', hist_counts', 15,...
        'linestyle', 'none'); 
    cbar = colorbar(); 
    cbar.Label.String = 'ln(times sampled)'; cbar.Label.FontSize = 12; 
    %%%

    linkaxes([axPdf,ax12 ]); 
    set(ax12, 'Position', axPdf.Position); 
    axes(axPdf); axes(ax12); 
    title(ax12, 'zmoh/vpvs sampling density','fontsize',titleSize, ...
        'fontweight', 'normal'); 
    xlabel(ax12, 'Vp/Vs ratio','fontsize',14, 'fontweight', 'normal'); 
    ylabel(ax12, 'Moho depth','fontsize',14, 'fontweight', 'normal'); 
    set(ax12,'linewidth', 2); 
    
    if ifsave
        exportgraphics(gcf, sprintf(...
            '%s/final_true_vs_pred_data_wavs_HK.png',par.res.resdir), ...
            'resolution', 300); 
    end
    
end
