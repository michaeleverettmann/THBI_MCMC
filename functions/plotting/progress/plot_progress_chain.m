function plot_progress_chain(absTimeIter, par, accept_info, ptb); 

if any(par.ii==[10 100 500 1000]) || (mod(par.ii, 2000) == 0) % Every thousand iterations or a few select iterations before that
    ifaccept  = [accept_info(1:par.ii).ifaccept       ]'; 
    loglik    = [accept_info(1:par.ii).log_likelihood ]'; 
    ptbnorm   = [accept_info(1:par.ii).ptbnorm        ]'; 
    iter      = [accept_info(1:par.ii).iter           ]'; 
    fail_chain= [accept_info(1:par.ii).fail_chain     ]'; 
    fail_total= [accept_info(1:par.ii).fail_total     ]'; 
%     ifpass    = [accept_info(1:par.ii).ifpass         ]'; 
 
    nRow = 5; nCol = 1; 
    figure(3001); clf; hold on; set(gcf, 'color', 'white', 'pos', [1000 673 800 664]);
    
    try 
        subplot(nRow, nCol, 1); hold on; box on; 
        scatter([absTimeIter.data], iter); 
        xlabel('Time'); 
        title('Iteration'); 
    catch
        subplot(nRow, nCol, 1); hold on; box on; 
        title('Problem making this iteration versus time plot...')
        warning('Couldnt plot time. Probably [absTimeIter.data] has wrong size. TODO! brb2022.05.17');
    end
    
    subplot(nRow, nCol, 2); hold on; box on; 
    scatter(iter, ifaccept); 
    xlabel('Iteration'); 
    title('Accepted or not'); 
    
    subplot(nRow, nCol, 3); hold on; box on; 
    scatter(iter, fail_chain); 
    xlabel('Iteration'); 
    title('Fail chain (at start of iteration)'); 
    ylabel('Fail chain'); 
    yyaxis right
    scatter(iter, fail_total); 
    ylabel('Fail total'); 
    
    subplot(nRow, nCol, 4); hold on; box on; 
    scatter(iter, loglik); 
    xlabel('iter'); 
    title('Log likelihood'); 
    prc10 = prctile(loglik, 10); 
    ylim([prc10, max(loglik)]); 
    ylabel('10th percentile to max'); 
    
    subplot(nRow, nCol, 5); hold on; box on; 
    scatter(iter, ptbnorm); 
    xlabel('Time'); 
    title('Norm of model change since last kernel reset');        
    
    exportgraphics(gcf, sprintf('%s/progress_%s.pdf',...
        par.res.resdir, par.res.chainID)); 
    
    if par.inv.verbose; %%% TODO can manually set this to true if you want to see it. 
        % Something about which perturbation happened. 
        % Slow to plot this and takes a lot of memory. Only do if verbose. 
        ptbs = string(ptb); 
        [all_ptb, first_occurance, ptb_num] = unique(ptbs); 

        figure(3002); set(gcf, 'pos', [719 711 1481 1346]); 
        clf; hold on; box on; grid on; 
        for i_allptb = [1:length(all_ptb)]; 
            scatter(find(ptb_num == i_allptb), i_allptb, 80, 'k', 'filled')
        end
        yticks([1:length(all_ptb)]); 
        yticklabels(all_ptb); 

        did_accept = ifaccept == 1; 
        scatter(find( did_accept), 0.5, 40, 'blue', 'filled')
        scatter(find(~did_accept), 0,   40, 'red' , 'filled')
        xlabel('Iteration'); 

        exportgraphics(gcf, sprintf('%s/progress_ptb_%s.pdf',...
            par.res.resdir, par.res.chainID)); 
    end
    
end   

end