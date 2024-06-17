function [LAB,MLD] = LAB_finder(Vs,Z,Zmoh,smthpts,ifplot,maxlab)
%LAB_FINDER Summary of this function goes here
%   Function to estimate LAB depth, using a variety of proxies. 
%   Along the way, calculate any MLDs
arguments
    Vs (:,1)
    Z (:,1)
    Zmoh (1,1) = Z(crossing(Vs,4.0))
    smthpts (1,1) = 10
    ifplot (1,1) = false;
    maxlab (1,1) = 200;
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
        ylim([50 270]),set(gca,'ydir','reverse')
        xlim([4.4 4.7])
        ylabel('depth, km'),zlabel('velocity (km/s)')
    end


    % method 1: stationary points of first derivative
    [invg,znvg] = crossing(d2vdz2,Z); % find maximum gradient points (turning)
    
    znvg(dvdz(invg)>0) = []; % have to be negative velocity gradients
    invg(dvdz(invg)>0) = []; % have to be negative velocity gradients
    
    znvg(Z(invg)>maxlab) = []; % no nvg deeper than 250
    invg(Z(invg)>maxlab) = []; % no nvg deeper than 250

    znvg(Z(invg)<Zmoh) = []; % no nvg deeper than moho
    invg(Z(invg)<Zmoh) = []; % no nvg deeper than moho
    [~,imax] = max(smthVs(1:mindex(Z,maxlab)));
    zmax = Z(imax);

    if any(invg<imax) 
        MLD(1).z_stat = znvg(find(invg<imax,1,'first'));
    else
        MLD(1).z_stat = nan;
    end    

    znvg(invg<imax) = [];% has to be deeper than the Vs maximum
    invg(invg<imax) = [];% has to be deeper than the Vs maximum
    
    znvg = znvg(mindex(dvdz(invg))); % choose the one associated with steepest nvg
    invg = invg(mindex(dvdz(invg))); % choose the one associated with steepest nvg
    LAB(1).z_stat = znvg(find(invg>imax,1,'first'));

    
    % method 2: velocity minimum
    [itrn,ztrn] = crossing(dvdz,Z); % find maximum gradient points (turning)
    ztrn(Z(itrn)<LAB.z_stat) = []; % must be deeper than method 1 lab
    itrn(Z(itrn)<LAB.z_stat) = []; % must be deeper than method 1 lab
    LAB.z_vmi = ztrn(find(itrn>imax,1,'first'));

    % similarly method 2 for MLD
    [itrn,ztrn] = crossing(dvdz,Z); % find maximum gradient points (turning)
    ztrn(d2vdz2(itrn)<0) = []; % must be turning faster - minimum
    itrn(d2vdz2(itrn)<0) = []; % must be turning faster - minimum
%     itrn(ztrn>zmax) = []; % must be shallower than max point
%     ztrn(ztrn>zmax) = []; % must be shallower than max point
    MLD(1).z_vmi = ztrn(find(itrn<imax,1,'first'));


    % method 3: halfway between max and min points
    [itrn,ztrn] = crossing(dvdz,Z); % find maximum gradient points (turning)

    itrn(ztrn<Zmoh) = []; % must be in mantle
    ztrn(ztrn<Zmoh) = []; % must be in mantle

    ztrn(itrn < (imax-1)) = [];
    itrn(itrn < (imax-1)) = [];

    imin = itrn(end);
    if length(itrn)>2, itrn = itrn(1:2);ztrn = ztrn(1:2);end
    LAB.z_vmi50ma  = 0.5*diff(ztrn) + zmax;
    if isempty(LAB.z_vmi50ma),LAB.z_vmi50ma=nan;end
    
    % method 4: 90% between max and min points
    itrn = crossing(dvdz,Z); % find maximum gradient points (turning)
    itrn(Z(itrn)<Zmoh) = []; % must be in mantle
    itrn(itrn < (imax-1)) = [];
    V_mami = Vs(itrn);

    if length(V_mami)<2 % there is no minimum point - deeper than base of model
        LAB.z_vmi90ma = nan; 
    else
        [~,zzo] = crossing(Vs,Z,V_mami(2)+0.1*(max(V_mami)-min(V_mami)));
        zzo(zzo<=Z(imax)) = [];
        zzo(zzo>=Z(imin)) = [];
        LAB.z_vmi90ma = zzo(1);
    end

    if ifplot
        for labtype = {'z_stat','z_vmi','z_vmi90ma','z_vmi50ma'}
            subplot(131),hold on
            yline(LAB.(labtype{1}),'c','linewidth',2)
        end
        for mldtype = {'z_stat','z_vmi'}
            subplot(131),hold on
            yline(MLD.(mldtype{1}),'m','linewidth',2)
        end
    end

end

