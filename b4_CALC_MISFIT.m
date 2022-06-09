function [ misfit ] = b4_CALC_MISFIT( trudata,predata,par,ifplot,SWwt,options )
    arguments
        trudata
        predata
        par
        ifplot = []
        SWwt = []
        options.plotRFError = false; 
    end
  
% Calculate cross-convolution misfit between observed data and
% (unit-normalised) predicted data

if nargin < 4 || isempty(ifplot)
    ifplot = 0;
end
if nargin < 5 || isempty(SWwt)
    SWwt = struct([]);
    for id = 1:length(par.inv.datatypes)
        if regexp(par.inv.datatypes{id},'SW')
            SWwt(1).(par.inv.datatypes{id}) = 1;
        end
    end
end

%% Setup misfit 
% this is the SUM OF SQUARES MISFIT VECTOR
misfit = struct('E2',cell2struct(cell(size(par.inv.datatypes)),par.inv.datatypes,2));



%% BW misfit
for id = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{id};
    pdt = parse_dtype(dtype);
    
    % BODY WAVE 
    % Calculate cross-convolution misfit between observed data and
    % (unit-normalised) predicted data
    if any(strcmp({'BW','RF'},pdt{1}))
        misfit2    = zeros(length(trudata.(dtype)),1);
        stfpow_tru = zeros(length(trudata.(dtype)),1);
        stfpow_pre = zeros(length(trudata.(dtype)),1);
                
        for itr = 1:length(trudata.(dtype))
            
            if strcmp(dtype, 'RF_Sp_ccp') % Upweight the deeper parts of ccp relative to shallow
                rf_tru = trudata.(dtype)(itr).PSV(:,1); 
                rf_pre = predata.(dtype)(itr).PSV(:,1);  
                zz     = trudata.(dtype).zz; 

                weight_depth_val = par.datprocess.CCP.weight_depth_val; 
%                 weight_depth_val = [-9999999,0.3 ; 50,0.3 ; 100,1 ; 9999999,1]; % First collumn: Specify depths of interest. Second collumn: Ideal weight at those depths. 
                weight = interp1(weight_depth_val(:,1), weight_depth_val(:,2), zz); % Linear interpolation between depth-weight pairs. 
                depth_smooth = 80; % Smooth linearly interpolated weight vector with a Gaussian window over this duration
                [~,ind_depth_smooth] = min( abs(zz-zz(1) - depth_smooth) ); % Find how many indicies correspond to traversing the specified smoothing depth
                weight = smoothdata(weight, 'gaussian', ind_depth_smooth); % Make smoothed weight vector. 
                weight = weight / rms(weight); % Normalize, to avoid messing with sigma of the data type. 
                rf_tru_weight = rf_tru .* weight; 
                rf_pre_weight = rf_pre .* weight; 
 
                trudata.(dtype)(itr).PSV(:,1) = rf_tru_weight; 
                predata.(dtype)(itr).PSV(:,1) = rf_pre_weight; 

                try %^ Plot -- put in try loop. No need to stop inversion if this plot doesn't work. 
                    if par.inv.verbose  || (par.ii == par.inv.niter); 
                        figure(1); clf; hold on; 
                        box on; 
                        set(gcf, 'color', 'white', 'pos', [-742 527 530 387]);
                        set(gca, 'linewidth', 1); 
                        xlabel('Depth (km)'); 
                        title('Receiver function weighting', 'fontweight', 'normal', ...
                            'fontsize', 16); 
                        lineWidths = 2; 


                        yyaxis left; 
                        yticks([]); 
                        shiftVal = 2 * rms(rf_tru); 
                        rf_tru_hand = plot(zz, rf_tru+shiftVal       , ...
                            '-k','linewidth', lineWidths); 
                        rf_pre_hand = plot(zz, rf_pre-shiftVal       , ...
                            '-k','linewidth', lineWidths);  
                        rf_tru_weight_hand = plot(zz, rf_tru_weight+shiftVal, ...
                            '-b','linewidth', lineWidths);  
                        rf_pre_weight_hand = plot(zz, rf_pre_weight-shiftVal, ...
                            '-b','linewidth', lineWidths);  

                        yyaxis right; 
                        ylabel('Weight'); 
                        weight_hand = plot(zz, weight, ...
                            ':k','linewidth', lineWidths); 

                        legend([rf_tru_hand, rf_tru_weight_hand, weight_hand], ...
                            'Original', 'Weighted','Weight'); 
                        yyaxis left; 
                        text(gca, 0,  shiftVal*1.3, 'Observed' , 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); 
                        text(gca, 0, -shiftVal*1.3, 'Predicted', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); 

                        ax = gca;ax.YAxis(1).Color = 'k';ax.YAxis(2).Color = 'k';

                        if (par.ii == par.inv.niter) 
                            exportgraphics(gcf, sprintf('%s/rf_weighting_%s.pdf',par.res.resdir,par.res.chainstr));
                        end
                    end
                catch e
                    warning('Could not plot rf weighting scheme. Delete this try catch when you know this plotting works.\n%s',getReport(e))
                end %^ end plot
            end          
            
            misfit2(itr)    =  xconv_misfit(trudata.(dtype)(itr).PSV(:,1), ...
                                            trudata.(dtype)(itr).PSV(:,2), ...
                                            predata.(dtype)(itr).PSV(:,1), ...
                                            predata.(dtype)(itr).PSV(:,2));
                                           
            stfpow_tru(itr) = norm(trudata.(dtype)(itr).PSV(:,1)) +...
                                 norm(trudata.(dtype)(itr).PSV(:,2));
            stfpow_pre(itr) = norm(predata.(dtype)(itr).PSV(:,1)) +...
                                 norm(predata.(dtype)(itr).PSV(:,2));
                             
                            
        end
              
        misfit.E2.(dtype) = misfit2./stfpow_tru./stfpow_pre;      % SUM OF SQUARED MISFITS, NORMALISED
        
        %%% brb2022.02.18 Simple rms estimate for comparison. As of now,
        %%% not in use for inverison. 
        if strcmp(dtype, 'RF_Sp') || strcmp(dtype, 'RF_Sp_ccp'); 
            truSV = trudata.(dtype)(itr).PSV(:,1); 
            preSV = predata.(dtype)(itr).PSV(:,1); 
            
%             truSV = truSV ./ norm(truSV); 
%             preSV = preSV ./ norm(preSV); 
            difSV = truSV - preSV; 
                        
            misfit.RMS_simple.(dtype) = sqrt(norm(truSV - preSV) ./ norm(truSV)); 
            
            if par.inv.verbose || options.plotRFError; % brb2022.02.18 Optional plotting to better understand RF error calculation. 
                figure(1002); clf; hold on; set(gcf, 'color', 'white'); 
                title(sprintf('RMS(diff)/norm(true) = %0.4f, X conv err = %0.4f', ...
                        misfit.RMS_simple.(dtype), misfit.E2.(dtype) )); 
                xlabel('Index'); 
                ylabel('Sp receiver functions'); 
                box on; 
                LW = 3; 
                shiftVal = 2*rms(truSV);
                plot(difSV         , 'r', 'DisplayName', 'Difference', 'LineWidth', LW/2); 
                plot(truSV-shiftVal, 'k', 'DisplayName', 'truSV',      'LineWidth', LW); 
                plot(preSV+shiftVal, 'b', 'DisplayName', 'preSV',      'LineWidth', LW); 
                legend(); 
                exportgraphics(gcf, [par.res.resdir '/misfit_sp_rf.pdf']); 

                figure(1001); clf; hold on; set(gcf, 'color', 'white'); 
                sgtitle(['Error (unsquared) = ' num2str(sqrt(misfit.E2.(dtype)))] ); 
                subplot(4,1,1); hold on; box on; 
                title(['True SV, norm = ' ...
                    num2str(norm(trudata.(dtype)(itr).PSV(:,1)))]); 
                plot(trudata.(dtype)(itr).PSV(:,1))
                subplot(4,1,2); hold on; box on; 
                title(['True P, norm = ' ...
                    num2str(norm(trudata.(dtype)(itr).PSV(:,2)))]); 
                plot(trudata.(dtype)(itr).PSV(:,2))
                subplot(4,1,3); hold on; box on; 
                title(['Pre SV, norm = ' ...
                    num2str(norm(predata.(dtype)(itr).PSV(:,1)))]); 
                plot(predata.(dtype)(itr).PSV(:,1))
                subplot(4,1,4); hold on; box on; 
                title(['Pre P, norm = ' ...
                    num2str(norm(predata.(dtype)(itr).PSV(:,2)))]);  
                plot(predata.(dtype)(itr).PSV(:,2))
            end
            
        end    
        

        
    % SURFACE WAVE
    elseif strcmp(pdt{1},'SW')
        e = (trudata.(dtype).(pdt{3}) - predata.(dtype).(pdt{3}));
        misfit.E2.(dtype) = e'*(SWwt.(dtype).*e); % SUM OF SQUARED MISFITS, NORMALISED

    % HKstack
    elseif strcmp(pdt{1},'HKstack')
        
        % Use external function to go from E to error. Will use at other places throughout the inversion. 
        misfit.E2.(dtype) = hKappaError(par.datprocess.HKappa.scale_error,...
                                        predata.(dtype).E_by_Emax, ...
                                        par.datprocess.HKappa.min_error); 
                                        
    end % end on data type
end

if ifplot
    plot_TRUvsPRE( trudata,predata)
end

end
