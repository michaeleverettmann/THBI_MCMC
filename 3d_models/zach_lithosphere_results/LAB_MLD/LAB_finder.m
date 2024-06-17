function [LAB,MLD] = LAB_finder(Vs,Z,Zmoh,smthpts,ifplot,minmlddv,minlabdv,mldmax)
% [LAB,MLD] = LAB_finder(Vs,Z,Zmoh,smthpts,ifplot,minmlddv,minlabdv,mldmax)
%LAB_FINDER Summary of this function goes here
%   Function to estimate LAB depth, using a variety of proxies. 
%   Along the way, calculate any MLDs
arguments
    Vs (:,1)
    Z (:,1)
    Zmoh (1,1) = Z(crossing(Vs,4.0))
    smthpts (1,1) = 5
    ifplot (1,1) = false;
    minmlddv (1,1) = 0.004;
    minlabdv (1,1) = 0.004;
    mldmax (1,1) = 120;
end
LAB = struct([]);
MLD = struct([]);

    smthVs = smooth(Z,Vs,smthpts);
    dvdz = gradient(smthVs);
    d2vdz2 = gradient(gradient(smthVs));

    % plotting while testing
    if ifplot
        figure(99); clf,hold on
        subplot(133), hold on, 
        plot(d2vdz2,Z,'linewidth',2)
        ylim([50 270]),set(gca,'ydir','reverse')
        xlim([-0.002 0.002]), xline(0)
        ylabel('depth, km'),xlabel('second derivative')
    
        subplot(132), hold on
        plot(dvdz,Z,'linewidth',2)
        xline(0)
        ylim([50 270]),set(gca,'ydir','reverse')
        xlim([-0.002 0.002]), 
        ylabel('depth, km'),xlabel('first derivative')
    
        subplot(131), hold on
        plot(smthVs,Z,'linewidth',2)
        plot(Vs,Z,'r','linewidth',2)
        ylim([50 270]),set(gca,'ydir','reverse')
        xlim([4.4 4.7])
        ylabel('depth, km'),zlabel('velocity (km/s)')
    end

    %% firstly, find all the significant max/min points
    [istat,zstat] = crossing(dvdz,Z);
    istat(zstat<Zmoh) = []; % must be in mantle
    zstat(zstat<Zmoh) = []; % must be in mantle
    % establish if maxima (+ive) or minima (-ive) 
    sstat = -sign(interp1(Z,d2vdz2,zstat));
    vstat = interp1(Z,smthVs,zstat); %get the actual velocities 
    % now loop through and group into clusters with small changes
    % goal is to eject max/min pts if only tiny velocity change, but need
    % to know which of several to eject
    if length(zstat)>4
%         keyboard
    end

    ikill = []; % min/max to kill
    maxivisit = 1;
    for imami = 2:length(zstat)
        if ismember(imami,ikill), continue; end
        if imami <= maxivisit, continue; end 
            iclust = [];
            maxivisit = max([maxivisit,imami]);
            dvi = abs(dv_by_v(vstat(imami),vstat(imami-1)));
        while dvi < 1e-3 % 0.1% threshold
            iclust = unique([iclust;imami;imami-1]);
            imami = imami+1; 
            maxivisit = max([maxivisit,imami]);
            if imami > length(zstat), break; end
            dvi = abs(dv_by_v(vstat(imami),vstat(imami-1)));
        end
        if length(iclust)>=2
            if sum(sstat(iclust)) < 0 % this is a wobble around a negative
                % choose the biggest negative
                [~,iclust_min] = min(vstat(iclust));
                iclust(iclust_min) = []; % remove biggest min from cluster
                ikill = unique([ikill;iclust]);
            elseif sum(sstat(iclust)) > 0 % this is a wobble around a positive
                % choose the biggest positive
                [~,iclust_max] = max(vstat(iclust));
                iclust(iclust_max) = []; % remove biggest max from cluster
                ikill = unique([ikill;iclust]);
            elseif sum(sstat(iclust)) == 0 % this is a kink in an otherwise monotonic trend
                ikill = unique([ikill;iclust]);
            end
        end % if on cluster type (min/max/kink)
    end % loop on indices
    % get rid of bad turning pts
    istat(ikill) = []; 
    zstat(ikill) = []; 
    sstat(ikill) = []; 
    vstat(ikill) = []; 

    % account for annoying case when there is a third stationary pair in
    % asthenosphere
    if length(istat)>4 % more than two pairs of turns
        vmas = vstat(sstat>0); % velocities of maxima
        if length(vmas)>2 % more than three maxima
            if dv_by_v(vmas(2),vmas(3)) > 0.02 %last maximum is actually super slow
                ikill = find(vstat == vmas(3)):length(istat);
                % get rid of bad turning pts
                istat(ikill) = []; 
                zstat(ikill) = []; 
                sstat(ikill) = []; 
                vstat(ikill) = []; 
            end
        end
    end

    % now have bounds on MLD and LAB depths
    % MLD has to be deeper than first max, shallower than last
    % LAB has to be deeper than last max
    zmld_min = min(zstat(sstat>0));
    zlab_min = max(zstat(sstat>0));
    zmld_max = zlab_min;

    %% method 1: stationary points of second derivative
    [invg,znvg] = crossing(d2vdz2,Z); % find maximum gradient points (turning)
    
    znvg(dvdz(invg)>0) = []; % have to be negative velocity gradients
    invg(dvdz(invg)>0) = []; % have to be negative velocity gradients
    
%     znvg(Z(invg)>maxlab) = []; % no nvg deeper than maxlab
%     invg(Z(invg)>maxlab) = []; % no nvg deeper than maxlab

    znvg(Z(invg)<Zmoh)  = []; % no nvg shallower than moho
    invg(Z(invg)<Zmoh) = []; % no nvg shallower than moho


    % SEEK MLD using stationary gradient point
    zmld = znvg(znvg > zmld_min & znvg < zmld_max);
    if isempty(zmld)
        MLD(1).z_stat = nan;
    else
        MLD(1).z_stat = zmld(1);
    end

    % SEEK LAB using stationary gradient point    
    invg(znvg < zlab_min) = [];% has to be deeper than the lab minimum
    znvg(znvg < zlab_min) = [];% has to be deeper than the lab minimum
    
    znvg = znvg(mindex(dvdz(invg))); % choose the one associated with steepest nvg
    if isempty(znvg)
        LAB(1).z_stat = nan;
    else
        LAB(1).z_stat = znvg(1);
    end

    
    %% method 2: velocity minimum
    % use previously found maximum gradient points (turning)     
    ztrn = zstat(sstat < 0);

    % SEEK MLD using minimum point
    zmld = ztrn(ztrn > zmld_min & ztrn < zmld_max);
    if isempty(zmld)
        MLD(1).z_vmi = nan;
    else
        MLD(1).z_vmi = zmld(1); % use shallowest if few
    end

    % SEEK LAB using minimum point
    zlab = ztrn(ztrn > zlab_min);
    zlab(zlab<LAB.z_stat) = []; % must be deeper than method 1 lab

    if isempty(zlab)
        LAB(1).z_vmi = nan;
    else
        LAB(1).z_vmi = zlab(1); % use shallowest if few
    end    
 
 
    %% method 3: halfway between max and min points
    % use previously found maximum gradient points (turning)     
    itrn = istat;
    ztrn = zstat;

    itrn(ztrn<Zmoh) = []; % must be in mantle
    ztrn(ztrn<Zmoh) = []; % must be in mantle

    % find pairs of max/min
    Npair = floor(length(itrn)/2);
    zpair = zeros(Npair,2);
    vpair = zeros(Npair,2);
    for ipair = 1:Npair
        zpair(ipair,:) = ztrn(ipair*2+[-1 0]);
        vpair(ipair,:) = smthVs(itrn(ipair*2+[-1 0]));
    end
    z_vmi50ma = mean(zpair,2);

    % SEEK MLD using halfway method
    zmld = z_vmi50ma(z_vmi50ma > zmld_min & z_vmi50ma < zmld_max);
    if isempty(zmld)
        MLD(1).z_vmi50ma = nan;
    else
        MLD(1).z_vmi50ma = zmld(1);
    end

    % SEEK LAB using halfway method
    zlab = z_vmi50ma(z_vmi50ma > zlab_min );
    if isempty(zlab)
        LAB.z_vmi50ma = nan;
    else
        LAB(1).z_vmi50ma = zlab(1);
    end
    
    %% method 4: 90% of *velocity* (not depth) between max and min points
    z90 = nan(Npair,1);
    for ipair = 1:Npair
        v90 = 0.1*vpair(ipair,1) + 0.9*vpair(ipair,2);
        izpair = (Z >= zpair(ipair,1)) & (Z <= zpair(ipair,2));
        if sum(izpair)==1
            z90(ipair) = mean(zpair(ipair,:)); % if literally touching line. Will get killed by later MLD test
        else
        [~,z90_provis] = crossing(smthVs(izpair),Z(izpair),v90);
        if ~isempty(z90_provis)
            z90(ipair) = z90_provis(1);
        end
        end
    end

    % SEEK MLD using v90 method
    zmld = z90(z90 > zmld_min & z90 < zmld_max);
    if isempty(zmld)
        MLD(1).z_vmi90ma = nan;
    else
        MLD(1).z_vmi90ma = zmld(1);
    end


    % SEEK LAB using v90 method
    zlab = z90(z90 > zlab_min );
    if isempty(zlab)
        LAB.z_vmi90ma=nan;
    else
        LAB(1).z_vmi90ma = zlab(1);
    end


%% test all the MLDs, make sure appreciably different velocity than max above of them

    MLDtypes = fieldnames(MLD);
    for imldt = 1:length(MLDtypes)
        if isnan(MLD.(MLDtypes{imldt})), continue; end
        %now find max values above this supposed mld
        vmld = interp1(Z,smthVs,MLD.(MLDtypes{imldt}));
        vmax_zlt_mld = vstat(zstat < MLD.(MLDtypes{imldt}));
%         vmax_zgt_mld = vstat(zstat > MLD.(MLDtypes{imldt}));
        dv_zlt_mld = (vmax_zlt_mld(1)-vmld)/vmax_zlt_mld(1); % fractional difference to shallower max
%         dv_zgt_mld = (vmax_zgt_mld-vmld)/vmld; % fractional difference to deeper max
%         dv_mld = max([dv_zgt_mld,dv_zlt_mld]);
    % has to be appreciably slower than the maximum ABOVE it
        dv_mld = dv_zlt_mld; 
        % kill mld if not 0.3% slower (bearing in mind this is "stat"
        % method
        if dv_mld < minmlddv
            MLD.(MLDtypes{imldt}) = nan;
        else
            MLD(1).(regexprep(MLDtypes{imldt},'z_','v_')) = vmld;
            MLD(1).(regexprep(MLDtypes{imldt},'z_','dv_')) = dv_mld;
        end
    end

    %% test all the LABS, make sure appreciably different velocity than max above them
    % if okay, record V and dV, for use later
    LABtypes = fieldnames(LAB);
    for ilabt = 1:length(LABtypes)
        if isnan(LAB.(LABtypes{ilabt})), continue; end
        %now find max values above this supposed lab
        vlab = interp1(Z,smthVs,LAB.(LABtypes{ilabt}));
        vmax_zlt_lab = vstat(zstat < LAB.(LABtypes{ilabt}) & sstat > 0);
        % has to be appreciably slower than the maximum closest ABOVE it
        dv_lab = (vmax_zlt_lab(end)-vlab)/vlab; % fractional difference to shallower max
%         dv_lab = dv_zlt_lab; 
        % kill lab if not 0.3% slower (bearing in mind this is "stat"
        % method
        if dv_lab < minlabdv
            LAB.(LABtypes{ilabt}) = nan;
        else
            LAB(1).(regexprep(LABtypes{ilabt},'z_','v_')) = vlab;
            LAB(1).(regexprep(LABtypes{ilabt},'z_','dv_')) = dv_lab;
        end
    end

    %% if mld is too deep and no LAB, conclude MLD is actually LAB
     if MLD.z_stat > mldmax  || MLD.z_vmi90ma > mldmax 
%         if ~isnan(LAB.z_stat)
%             keyboard
%         end
        LAB = MLD; 
        MLD = struct('z_stat',nan,'z_vmi',nan,'z_vmi50ma',nan,'z_vmi90ma',nan);
    end




    if ifplot
        for labtype = {'z_stat','z_vmi','z_vmi90ma','z_vmi50ma'}
            if ~isnan(LAB.(labtype{1}))
            subplot(131),hold on
            yline(LAB.(labtype{1}),'c','linewidth',2)
            end
        end
        for mldtype = {'z_stat','z_vmi','z_vmi90ma','z_vmi50ma'}
            if ~isnan(MLD.(mldtype{1}))
            subplot(131),hold on
            yline(MLD.(mldtype{1}),'m','linewidth',2)
            end
        end
    end

end

function dv = dv_by_v(v1,v2)
    dv = 2*(v1 - v2)/(v1 + v2);
end
