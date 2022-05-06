function plot_progress_chain(absTimeIter, par, accept_info); 




if any(par.ii==[10 100 500 1000]) || (mod(par.ii, 1000) == 0) % Every thousand iterations or a few select iterations before that
    ifaccept           = [accept_info.ifaccept      ]; 
    loglik             = [accept_info.log_likelihood]; 
    ptbnorm            = [accept_info.ptbnorm       ]; 
    iter               = [accept_info.iter          ]; 
    
    nRow = 4; nCol = 1; 
    figure(3001); clf; hold on; set(gcf, 'color', 'white');
    
    subplot(nRow, nCol, 1); hold on; box on; 
    scatter([absTimeIter.data], iter); 
    xlabel('Time'); 
    title('Iteration'); 
    
    subplot(nRow, nCol, 2); hold on; box on; 
    scatter(iter, ifaccept); 
    xlabel('Iteration'); 
    title('Accepted or not'); 
    
    subplot(nRow, nCol, 3); hold on; box on; 
    scatter(iter, loglik); 
    xlabel('Time'); 
    title('Log likelihood'); 
    
    subplot(nRow, nCol, 4); hold on; box on; 
    scatter(iter, ptbnorm); 
    xlabel('Time'); 
    title('Norm of model change since last kernel reset'); 
    
    exportgraphics(gcf, sprintf('%s/progress_%s.pdf',...
        par.res.resdir, par.res.chainID)); 
    
end   

end