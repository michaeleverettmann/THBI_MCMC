function [h,hline,hpos,hneg] = plot_RF( ax,RF,tt,level,xshift,iffill,pcol,ncol )
%[h,hline,hpos,hneg] = plot_RF( ax,RF,tt,level,xshift,iffill )
%
% function to plot a single receiver function on a set of axes "ax". This
% will plot the receiver function vector "RF" in time "tt", shifted in the
% x-direction by "xshift". If "iffill" is true, it will fill in the regions
% with greater amplitude than "level" - red for positive, blue for
% negative. The handles to the plotted elements are returned, where
%  h = [hline,hpos,hneg].
%
% ** Must be in increasing time!!! Flip RF data if not!! **


if nargin < 4 || isempty(level)
    level = maxab(RF)/1000;
end
if nargin < 5 || isempty(xshift)
    xshift = 0;
end
if nargin < 6 || isempty(iffill)
    iffill = true;
end
if nargin < 7 || isempty(pcol)
    pcol = [1 0 0];
end
if nargin < 8 || isempty(ncol)
    ncol = [0 0 1];
end

RF = [0;RF(:);0]; % pad edges to avoid beginnning/end plot tie-line issues
tt = [tt(1);tt(:);tt(end)]; % pad edges to avoid beginnning/end plot tie-line issues

dt=tt(2)-tt(1); % edit by HEK May 2018
if dt < 0
    RF = flipud(RF);
    tt = flipud(tt);
end

% find and account for nans, if needed
rfnan = isnan(RF);
RFn = RF;
RF(rfnan) = 0;

RF = [0;RF(:);0]; % pad edges to avoid beginnning/end plot tie-line issues
tt = [tt(1);tt(:);tt(end)]; % pad edges to avoid beginnning/end plot tie-line issues

%% fill patches
axes(ax);hold on

if iffill
    RF_pos = RF;
    RF_pos(RF<level) = level;
    RF_neg = RF;
    RF_neg(RF>-level) = -level;
    % plot
    hold on
    hpos = patch(ax,RF_pos+xshift,tt,pcol,'linestyle','none','facealpha',0.6);
    hneg = patch(ax,RF_neg+xshift,tt,ncol,'linestyle','none','facealpha',0.6);
else
    hpos = [];
    hneg = [];
    fprintf('\nChose NO FILL\n');
end

hline = plot(ax,RF+xshift,tt,'-k','linewidth',1.5);

h = {hline,hpos,hneg};

end
