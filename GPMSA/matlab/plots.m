function cipdssMortPlots(pout,pvec, plotnum, varargin)

writeOut = 0;
newFig=1;
for ii=1:2:length(varargin)
  switch varargin{ii}
  case 'newFig'
    newFig=varargin{ii+1};
  case 'compareCal'
    poutL=varargin{ii+1};
  case 'writeOutput';
    suffix = varargin{ii+1};
    writeOut=1;
  otherwise
    error('invalid extended argument passed to cosmoPlots');
  end
end

model=pout.model;
data=pout.data;
pvals=pout.pvals(pvec);
nreal=length(pvals);
simData=pout.simData;

pu=model.pu; pv=model.pv;
p=model.p; q=model.q;

doPlot(100)=0;
if exist('plotnum'); doPlot(plotnum)=1;
else
end

% Set up prediction grid
gridu=0:0.075:1; glen=length(gridu);  mlen=glen*glen;
[gridx,gridy]=meshgrid(gridu,gridu);



if doPlot(1);
    if newFig; figure(1); clf; end;
    % Collect the MCMC record of the betas
      bu=[pvals.betaU]';
    % Transform them into rhos
      ru = exp(-bu*0.25);  
%       pu=10;
      off=0;
    for ii=1:pu
        r=ru(:,(ii-1)*(p+q)+2:ii*(p+q));
%         r=ru(:,(ii+off-1)*p+1:(ii+off)*p);
        D2subplot(ii,pu,1,[.15 .15 .1 .02],[.01 0 .01 0]);
%         subplot(pu,3,ii+off);
        if ii~=pu
            boxplot(r, 'labels',{'',''});
        else
            boxplot(r,'labels',{'eta','alpha'})
%             xlabel('k');
        end
        a=axis; axis([a(1) a(2) 0 1]);
%         ylab = ['\rho_{w',num2str(ii+off),'k}'];
%         ylabel(ylab,'FontSize',14); 
%         if ii~=6
%             xticklabels('');
%         else
%             xticks(1:6);
%         end
%         if ii==6
%             xlabel('k');
%         end
        ylim([0 1.01]);
%         if(ii==1)
%             %axes('position',[0 0 1 1]);
%             D2subplot(1,pu,1,[.4 .15 .1 .02],[.01 0 0 0]);
%             set(gca,'Visible','off');
%             %text((.5+(0:4))/5,ones([1 q])+0.08,xlabs,'fontSize',14);
%         end
    end
    
%     D2subplot(1,pu,1,[.4 .12 .1 .02],[.01 0 0 0]);
%     set(gca,'Visible','off');
%     text(-.12,.5,'PC1','FontSize',12,'HorizontalAlignment','center','Rotation',90);
%     text(.09,1.6,'\theta_1','HorizontalAlignment','center','FontSize',12);
%     text(.09+0.15,1.6,'\theta_2','HorizontalAlignment','center','FontSize',12);
%     text(.09+0.3,1.6,'\theta_3','HorizontalAlignment','center','FontSize',12);
%     text(.09+0.45,1.6,'\theta_4','HorizontalAlignment','center','FontSize',12);
%     text(.09+0.6,1.6,'\theta_5','HorizontalAlignment','center','FontSize',12);
%     text(.09+0.75,1.6,'\theta_6','HorizontalAlignment','center','FontSize',12);
%     D2subplot(2,pu,1,[.39 .12 .1 .02],[.01 0 0 0]);
%     set(gca,'Visible','off');
%     text(-.12,-.25,'PC2','FontSize',12,'HorizontalAlignment','center','Rotation',90);
%     D2subplot(3,pu,1,[.3 .12 .1 .02],[.01 0 0 0]);
%     set(gca,'Visible','off');
%     text(-.12,-.35,'PC3','FontSize',12,'HorizontalAlignment','center','Rotation',90);
%     D2subplot(4,pu,1,[.24 .12 .1 .02],[.01 0 0 0]);
%     set(gca,'Visible','off');
%     text(-.12,-.45,'PC4','FontSize',12,'HorizontalAlignment','center','Rotation',90);
%     D2subplot(5,pu,1,[.2 .12 .1 .02],[.01 0 0 0]);
%     set(gca,'Visible','off');
%     text(-.12,-.55,'PC5','FontSize',12,'HorizontalAlignment','center','Rotation',90);
    
    %D2subplot(1,q,1,[0.4 0.1 0.1 0],[.01 0 0 0]);
    %set(gca,'Visible','off');
end

% bivariate calibration density plot

if doPlot(2)
    if newFig; figure(2); clf; end
    t=[pout.pvals(min(pvec):max(pvec)).theta]';
%     t=pout.simData.orig.designNative;
    [nreal numt] = size(t);
    thetamax = max(pout.simData.orig.designNative); 
    thetamin = min(pout.simData.orig.designNative); 
    
    delta = thetamax - thetamin;
    tnative = zeros(size(t));
    for ii=1:numt;
        tnative(:,ii) = t(:,ii).*delta(ii) + thetamin(ii);
    end
    
    for ii=1:numt; 
        for jj=1:numt;
            if(ii==jj)
                D2subplot((ii-1)*numt+ii,numt,numt,[.1 .1 .08 .16],[0 0 0 0]);
                xout = linspace(0,1,100); 
                xoutnative = xout.*delta(ii)+thetamin(ii);
                sdx = 1/15;
                dx = dens1d(t(:,ii),sdx,xout);
                plot(xoutnative,dx,'b'); set(gca,'YtickLabel','', 'XtickLabel','');
                xlim([thetamin(ii) thetamax(ii)]); ylim([0 6.5]);
            else
                D2subplot((ii-1)*numt+jj,numt,numt,[.1 .1 .08 .16],[0 0 0 0]);
                mins=[thetamin(jj) thetamin(ii)];
                maxes=[thetamax(jj) thetamax(ii)];

                hpd2d01(t(:,[jj ii]),.10,'b',mins,maxes); hold on;
%                 plot(pout.obsData.orig.thetatrueNative(jj),...
%                     pout.obsData.orig.thetatrueNative(ii),'.r',...
%                    'MarkerSize',20);
            end
            if(mod(ii,2)==0) & (jj==numt);
                set(gca,'YAxisLocation','right');
            elseif(mod(ii,2)==1) & (jj==1) & (ii > 1);
                set(gca,'YAxisLocation','left');
            elseif(ii==numt) & (mod(jj,2)==1);
                set(gca,'XAxisLocation','bottom','YtickLabel','');
            elseif(ii==1) & (mod(jj,2)==0);
                set(gca,'XAxisLocation','top','YtickLabel','');
            elseif(ii==numt) & (jj==numt);
                set(gca, 'XtickLabel','');
            else
                set(gca,'XtickLabel','','YtickLabel','','Clim',[0 .01]);  
            end
            set(gca,'Clim',[0 .01]);
                
        end
    end
    % now put in the labels
    xlabs = {'\eta','\alpha','sza'};
    D2subplot(1,1,1,[.1 .1 .08 .16],[0 0 0 0]);
    set(gca,'Visible','off');
    for ii=1:numt
        posx=(.4 + ii-1)/numt;
        posy=(.5 + numt-ii)/numt;
        text(posx,posy,xlabs{ii},'FontSize',14);
        %text(posx,posy,num2str(ii),'FontSize',14);
    end
    drawnow
end    
% figure(2); print -deps2c cosmopost.eps;

if doPlot(3)
    % Now the realizations over each theta
    if newFig; figure(3); end;
    clf; colormap('spring');
    gridu=0:0.1:1; glen=length(gridu);  mlen=glen*glen;
    [gridx,gridy]=meshgrid(gridu,gridu);
    npred = glen;
    AzEl=[45 10];    
    plabs = {'n','h','\sigma_8','\Omega_{CDM}','\Omega_b'};
    for ii=1:q
        plabs(ii) = {['\theta_',num2str(ii)]};
    end
    % Ranges
    % n       = 0.8 - 1.4
    % H       = 50 - 110
    % sigma_8 = 0.6 - 1.6
    % omega_m = 0 - 0.6
    % omega_b = 0.02 - 0.12
%     plims = [pout.simData.orig.colmin' pout.simData.orig.colmax'];
    plims = [min(simData.orig.designNative(1,:)) max(simData.orig.designNative(1,:));
                min(simData.orig.designNative(2,:)) max(simData.orig.designNative(2,:)) ];
    
    prange = diff(plims');
    for ii=2:(p+q)
        tdat=ones(npred,p+q)*0.5;
        tdat(:,ii)=gridu';
        xpred=tdat(:,1:p);
        theta=tdat(:,p+1:end);
        %pred=gPred(xpred,pvals,model,data,'wpred',theta);
        pred=gPredict(xpred,pvals,model,data,'mode','uvpred','theta',theta);
        pm=squeeze(mean(pred.u,1));
        
        r=(pout.simData.Ksim*pm)'*pout.simData.orig.ysd;
        v=r + repmat(pout.simData.orig.ymean(:)',npred,1);
        
        D2subplot(ii-1,2,3,[.1 .1 .02 .01],[0 0 0 0]);
%         mesh(pout.simData.orig.time,grid,r,'MeshStyle','row','LineWidth',1.5); view(AzEl);
%         mesh(repmat(pout.simData.orig.time,size(gridu)),...
%             repmat(gridu.*prange(ii-1)+plims(ii-1,1),size(pout.simData.orig.time)),r',...
%             'FaceLighting','gouraud','LineWidth',1.5); view(AzEl);

        plot3(repmat(pout.simData.orig.time,size(gridu)),...
            repmat(gridu.*prange(ii-1)+plims(ii-1,1),size(pout.simData.orig.time)),v','b'); view(AzEl); 
        if ii<=p; var=['x_' num2str(ii)];else; var=['\theta_' num2str(ii-p)];end
%         ylabel(plabs(ii-1)); 
        set(gca,'Fontsize',8,'XColor',[.3 .3 .3]);
        if(ii==2) || (ii==5); 
            zlabel('\Delta ln(n)','FontSize',12); 
            xlabel('time');
            set(gca, 'Ztick', [-4 -2 0 2 4], 'Xtick', [10 20 30 40 50]);
        else
            set(gca, 'Zticklabel', ',', 'Xticklabel', ','); 
        end;
%         if(ii~=2) || (ii~=4); 
%             set(gca, 'Zticklabel', ',', 'Xticklabel', ','); 
%         else
%             set(gca, 'ZTick', [-4 -2 0 2 4], 'XTick', [5 10 15 20]);
%         end;
%         if(ii > 2); 
%             set(gca,'Zticklabel',' | ','Ztick',[0 1],'Xtick',[-2 -1]);
%         end
%         if(ii == 2); 
%             set(gca,'Ztick',[0 1],'Xtick',[-2 -1]);
%         end
%         if(ii == 6); 
%             set(gca,'Ytick',[.02 .06 .08]);
%         end
        alpha(0.25)
        axis tight; 
        set(gca,'Xgrid','on','Ygrid','on','Zgrid','on','Zlim',[0 3]);
        %set(gca,'Xgrid','on','Ygrid','on','Zgrid','on');

        drawnow
    end
%     axisNorm('start');
%       for ii=2:p+q;subplot(3,3,ii);axisNorm('add');end
%     axisNorm('set');
%     drawnow
end

if doPlot(31)
    
    % Now the realizations over each theta
    
    % keep theta's at their posterior mean
    theta_mean = mean([pout.pvals.theta], 2);    
    
    if newFig; figure(31); end;
    clf; colormap('spring');
    gridu=0:0.2:1; glen=length(gridu);  mlen=glen*glen;
    [gridx,gridy]=meshgrid(gridu,gridu);
    npred = glen;
    for ii=1:q
        plabs(ii) = {['\theta_',num2str(ii)]};
    end
    plims = [pout.simData.orig.colmin' pout.simData.orig.colmax'];
    prange = diff(plims');
    col_mat = repmat(linspace(0.1, 0.95, 6)', 1, 3);
    for ii=2:(p+q)
        tdat=ones(npred,p+q)*0.5;
        tdat=[ones(npred, 1) repmat(theta_mean', [npred 1])];
        tdat(:,ii)=gridu';
        xpred=tdat(:,1:p);
        theta=tdat(:,p+1:end);
        %pred=gPred(xpred,pvals,model,data,'wpred',theta);
        pred=gPredict(xpred,pvals,model,data,'mode','uvpred','theta',theta);
        pm=squeeze(mean(pred.u,1));
        
        r=(pout.simData.Ksim*pm)'*pout.simData.orig.ysd;
        v=r + repmat(pout.simData.orig.ymean(:)',npred,1);
        
        ax = D2subplot(ii-1,2,3,[0 0.1 0 0.1],[0.05 0 0.05 0]);
%         ax = subplot_tight(2,3,ii-1,[-0.1 0.01]);
        pbaspect(ax, [3 2.5 1]);
        set(gca, 'ColorOrder', col_mat, 'NextPlot', 'replacechildren');
        plot(v');
        if(ii==2) || (ii==5)
            ylabel('log(cumulative infected)','FontSize',10); 
            set(gca, 'Ytick', [0 10], 'Xtick', [10 30 50]);
        else
            set(gca, 'Yticklabel', '', 'Xtick', [10 30 50]); 
        end;
        if (ii==3) || (ii==6)
            xlabel('weeks')
        else
            xlabel('');
        end
        alpha(0.25)
        axis([0 57 0 15])
        drawnow
    end
    D2subplot(1,1,1,[.1 .1 .08 .16],[0 0 0 0]);
    set(gca,'Visible','off');
    text(.17,1.02,'\theta_1', 'FontSize', 14);
    text(.53,1.02,'\theta_2', 'FontSize', 14);
    text(.89,1.02,'\theta_3', 'FontSize', 14);
    text(.17,0.41,'\theta_4', 'FontSize', 14);
    text(.53,0.41,'\theta_5', 'FontSize', 14);
    text(.89,0.41,'\alpha', 'FontSize', 14);
end
% figure(3); print -deps2c cosmosens.eps;

if doPlot(4)    
    if newFig; figure(4); end
    clf
    for ii=1:length(pout.obsData)
        D2subplot(0+ii,1,3,[.09 .04 .06 .01],[0.06 0.05 0 0]);
        pred=gPred(pout.obsData(ii).x(1),pvals,model,data,'uvpred');
        eta = pout.simData.Ksim*pred.u'.*pout.simData.orig.ysd;
        etabounds = prctile(eta,[5 95],2);
        yesd = sqrt(diag(pout.obsData(ii).Sigy)).*pout.simData.orig.ysd;
        for k=1:length(yesd) if(yesd(k) > 100.0) yesd(k) = NaN; end; end;

        plot(pout.obsData(ii).orig.time,pout.obsData(ii).orig.y,'bo-','MarkerSize',4);
        hold on;
%         errorbar(pout.obsData(ii).orig.time,pout.obsData(ii).orig.y,yesd,'b');
        axis([0 90 0 4]);grid on;
%         set(gca,'XtickLabel','');
%         if(ii ==1) 
%             set(gca,'Ytick',[-.1 0 .1]); 
%         else
%             set(gca,'YtickLabel','');
%         end
        meanmat = repmat(pout.simData.orig.ymean,[1 2]);
        plot(pout.simData.orig.time,etabounds+meanmat,'g','LineWidth',1);
        ylabel('EUVAC','FontSize',10);
        title('Calibrated simulator','Fontsize',10);
        
%         if(writeOut)
%             etafit = prctile(eta,[5 50 95],2);
%             fname = ['calfits' suffix];
%             save(fname,'etafit','-ASCII');
%         end

        
        % now delta
        D2subplot(1+ii,1,3,[.09 .04 .06 .01],[0.06 0.05 0 0]);
        pv = size(pout.obsData(ii).Dobs,2);
        deltaR = pout.simData.orig.Dsim*pred.v'.*pout.simData.orig.ysd;
        deltaRbounds = prctile(deltaR,[5 95],2);
        plot(pout.simData.orig.time,deltaRbounds,'c','LineWidth',0.5); hold on;
        axis([0 90 -1 1]);grid on;
%         set(gca,'YtickLabel','');
        xlabel('time','FontSize',12);
        title('Discrepancy','Fontsize',10);
%         axis([0 12.5 -.25 .22]);
%         set(gca,'XtickLabel','');
%         if(ii ==1) 
%             set(gca,'Ytick',[-.1 0 .1]); 
%         else
%             set(gca,'YtickLabel','');
%         end
        
        % now a discrepancy-adjusted prediction (right)
        yhat = deltaR+eta;
        D2subplot(2+ii,1,3,[.09 .04 .06 .01],[0.06 0.05 0 0]);
        yhatbounds = prctile(yhat,[5 95],2);
        plot(pout.obsData(ii).orig.time,pout.obsData(ii).orig.y,'bo-','MarkerSize',4);
        hold on;
%         errorbar(pout.obsData(ii).orig.time,pout.obsData(ii).orig.y,yesd,'b');        
        %plot(pout.simData.orig.time,pout.simData.orig.yho(:,17),'k--');
        axis([0 90 0 4]);grid on;
%         set(gca,'YtickLabel','');
        title('Calibrated prediction','Fontsize',10);
%         axis([0 12.5 -.22 .22]);
%         set(gca,'Xtick',[0 5 10]);
%         if(ii ==1) 
%             set(gca,'Ytick',[-.1 0 .1]); 
%         else
%             set(gca,'YtickLabel','');
%         end
        plot(pout.simData.orig.time,yhatbounds+meanmat,'k','LineWidth',1);
    end
    
    if newFig; figure(41); end
    clf
    for ii=1:length(pout.obsData)
        D2subplot(0+ii,1,1,[.09 .04 .06 .01],[0.06 0.05 0 0]);
        pred=gPred(pout.obsData(ii).x(1),pvals,model,data,'uvpred');
        eta = pout.simData.Ksim*pred.u'.*pout.simData.orig.ysd;
        etabounds = prctile(eta,[5 95],2);
        yesd = sqrt(diag(pout.obsData(ii).Sigy)).*pout.simData.orig.ysd;
        for k=1:length(yesd) if(yesd(k) > 100.0) yesd(k) = NaN; end; end;

        plot(pout.obsData(ii).orig.time,pout.obsData(ii).orig.y,'bo','MarkerSize',4);
        hold on;
%         errorbar(pout.obsData(ii).orig.time,pout.obsData(ii).orig.y,yesd,'b');
        axis([0 90 0 4]);grid on;
%         set(gca,'XtickLabel','');
%         if(ii ==1) 
%             set(gca,'Ytick',[-.1 0 .1]); 
%         else
%             set(gca,'YtickLabel','');
%         end
        meanmat = repmat(pout.simData.orig.ymean,[1 2]);
        plot(pout.simData.orig.time,etabounds+meanmat,'g','LineWidth',1);
        ylabel('absorption (in dB)','FontSize',10);
        title('Calibrated simulator','Fontsize',10);
    end
 
end

% a plot for the proposal
if doPlot(41)    
    if newFig; figure(41); end
    clf
    % read in full observation set
    %obsdat = textread('../Data/obs.txt');
    nreal = length(pvals);
    
    for ii=1:length(pout.obsData)
        D2subplot(0+ii,1,3,[.09 .04 .06 .01],[0 0.05 0.05 0]);
        pred=gPred(pout.obsData(ii).x(1),pvals,model,data,'uvpred');
        eta = pout.simData.Ksim*pred.u'.*pout.simData.orig.ysd;
        etar = eta + repmat(pout.simData.orig.ymean,[1 nreal]);
        etabounds = prctile(eta,[5 95],2);
        yesd = sqrt(diag(pout.obsData(ii).Sigy)).*pout.simData.orig.ysd;
        yplot = pout.obsData(ii).orig.y;
        for k=1:length(yesd) 
            if(yesd(k) > 100.0) 
                yesd(k) = NaN; 
                yplot(k) = NaN;
            end; 
        end;

        plot(pout.simData.orig.time,pout.simData.orig.y,'-','color',[.8 .8 .8]);
        hold on;
        plot(pout.obsData(ii).orig.time,yplot,'ko','MarkerSize',4);
        errorbar(pout.obsData(ii).orig.time,yplot,yesd,'k');
        %plot(obsdat(13:end,1),log(obsdat(13:end,2)),'--','color',[.2 .2 .2]);
        axis tight; grid on;
%         set(gca,'XtickLabel','');
%         if(ii ==1) 
%             set(gca,'Ytick',[-.1 0 .1]); 
%         else
%             set(gca,'YtickLabel','');
%         end
        meanmat = repmat(pout.simData.orig.ymean,[1 2]);
        plot(pout.simData.orig.time,etabounds+meanmat,'k','LineWidth',2,'Color','c');
        ylabel('log(cumulative infected)','FontSize',10);
        xlabel('weeks','FontSize',10);
        title('Calibrated prediction','Fontsize',10);
    end
    
    % now a plot of the realizations, and then their max time
    D2subplot(2,1,3,[.09 .04 .06 .01],[0 0.05 0.05 0]);
    plot(pout.simData.orig.time,etar,'-','color',[.2 .2 .2]);
    ylabel('ln Infected','FontSize',10);
    xlabel('time (weeks)','FontSize',10);
    title('Posterior realizations','Fontsize',10);
 
    % compute the maxes of differenced curves
    detar = [exp(etar(1,:)); diff(exp(etar))];
    [mm ii] = max(detar);
    D2subplot(3,1,3,[.09 .04 .06 .01],[0 0.05 0.05 0]);
    plot(pout.simData.orig.time(1:end),detar,'-','color',[.2 .2 .2]);
    hold on;
    plot(ii,mm,'.','MarkerSize',10);
    ylabel('Weekly Infected','FontSize',10);
    xlabel('time (weeks)','FontSize',10);
    title('Posterior Peak Realizations','Fontsize',10);
    axis([0 60 -0.5*(10^4) max(mm)+10000]);
%     set(gca, 'Ytick', [0e4 2e4 4e4 6e4]);
%     set(gca, 'YtickLabel', [0e4 2e4 4e4 6e4]);
    
    if(1)
        etafit = prctile(eta,[5 50 95],2);
        fname = 'realizations.txt';
        save(fname,'etar','-ASCII');
    end
    
end

% a modified doPlot(41)
if doPlot(42)    
    if newFig; figure(42); end
    clf
    % read in full observation set
    %obsdat = textread('../Data/obs.txt');
    nreal = length(pvals);
    
    for ii=1:length(pout.obsData)
        pred=gPred(pout.obsData(ii).x(1),pvals,model,data,'uvpred');
        eta = pout.simData.Ksim*pred.u'.*pout.simData.orig.ysd;
        etar = eta + repmat(pout.simData.orig.ymean,[1 nreal]);
        etabounds = prctile(etar,[5 95],2);
        yesd = sqrt(diag(pout.obsData(ii).Sigy)).*pout.simData.orig.ysd;
        yplot = pout.obsData(ii).orig.y;
        for k=1:length(yesd) 
            if(yesd(k) > 100.0) 
                yesd(k) = NaN; 
                yplot(k) = NaN;
            end; 
        end;

        
        %plot(obsdat(13:end,1),log(obsdat(13:end,2)),'--','color',[.2 .2 .2]);
        
%         set(gca,'XtickLabel','');
%         if(ii ==1) 
%             set(gca,'Ytick',[-.1 0 .1]); 
%         else
%             set(gca,'YtickLabel','');
%         end
        meanmat = repmat(pout.simData.orig.ymean,[1 2]);
    end
    
    % now a plot of the realizations, and then their max time
    D2subplot(1,1,3,[.13 .0 .06 .01],[0.05 0.1 0.05 0]);
    plot(pout.simData.orig.time,pout.simData.orig.y,'-','color',[.8 .8 .8]);
    hold on;
    plot(pout.simData.orig.time,etar,'-','color','c');
    plot(pout.simData.orig.time,etabounds,'--','color','k','Linewidth',1.5);
    plot(pout.obsData(ii).orig.time,yplot,'ko','MarkerSize',4);
%     errorbar(pout.obsData(ii).orig.time,yplot,yesd,'k');
    axis tight; grid on;
    ylabel('log(cumulative infected)','FontSize',10);
    xlabel('weeks','FontSize',10);
    %title('Posterior realizations','Fontsize',11);
 
    % compute the maxes of differenced curves
    detar = [exp(etar(1,:)); diff(exp(etar))];
    [mm ii] = max(detar);
    kk = datasample(1:nreal, 50, 'Replace', false);
    D2subplot(2,1,3,[.13 .0 .06 .01],[0.05 0.1 0.05 0]);
    plot(pout.simData.orig.time(1:end),detar(:,kk),'-','color','c');
    hold on;
    plot(ii(kk),mm(kk),'.','MarkerSize',10,'color','b');
    ylabel('Weekly Infected','FontSize',10);
    xlabel('weeks','FontSize',10);
    %title('Posterior Peak Realizations','Fontsize',11);
    axis([0 90 0 max(mm(kk))+100]); grid on;
    set(gca, 'YTick',[0 500 1000]);
%     set(gca, 'Ytick', [0e4 2e4 4e4 6e4]);
%     set(gca, 'YtickLabel', [0e4 2e4 4e4 6e4]);

    % plot a histogram of max calculated above
    D2subplot(3,1,3,[.13 .0 .06 .01],[0.05 0.1 0.05 0]);
    histogram(ii, 'BinWidth', 3, 'BinLimits', [0,50], 'Normalization','count', ...
        'FaceColor','b');
    ylabel('Count','FontSize',10);
    xlabel('weeks','FontSize',10);
    %title('Histogram of Peak Realizations','Fontsize',11);
    axis([0 57 0 1200]); grid on;
    set(gca, 'YTick',[0 400 800 1200]);
    
    D2subplot(1,1,1,[.0 .06 .0 .01],[0.05 0.06 0 0]);
    set(gca,'Visible','off');
    text(.01,0.96,'Posterior prediction', 'FontSize', 10);
    text(.42,0.96,'Posterior peak', 'FontSize', 10);
    text(.75,0.96,'Histogram of peak timing', 'FontSize', 10);
    
    if(1)
        etafit = prctile(eta,[5 50 95],2);
        fname = 'realizations.txt';
        save(fname,'etar','-ASCII');
    end
    
end

% a plot of the predicted w(theta) surface for each PC; use theta 1 and 2
% here as the varying grid.
if doPlot(5)    
    if newFig; figure(5); end
    % Set up prediction grid
    gridu=0:0.1:1; glen=length(gridu);  mlen=glen*glen;
    [gridx,gridy]=meshgrid(gridu,gridu);
    xxx = [.15 .15 .33 .56];
    npc = size(pout.simData.Ksim,2);
    clf
    %%%%%%%%
    npred=length(gridx(:));
    AzEl=[45 55];  
    % make the prediction...
    tdat=ones(npred,p+q)*0.5;
    tdat(:,2)=gridx(:); tdat(:,3)=gridy(:); %change column index for diff theta
    xpred=tdat(:,1:p);
    theta=tdat(:,p+1:end);
    %pred=gPred(xpred,pvals,model,data,'wpred',theta);
    pred=gPredict(xpred,pvals,model,data,'mode','uvpred','theta',theta);
    pm=squeeze(mean(pred.u,1));
    
    % now the plots...
    for ii=1:npc
        D2subplot(0+ii,2,3,[0 0 0 0],[.12 .1 .06 .01]);
        w1 = reshape(pm(ii,:),[11 11]);
        mesh(gridx,gridy,w1); hold on; axis([0 1 0 1 -3 3]);
        xlabel('\theta_3','fontSize',11);
        ylabel('\theta_1','fontSize',11);
        zlabel(['w_' num2str(ii) '(\theta)'],'fontSize',11);
        title(['PC ' num2str(ii)]);
    end
    % now make the decomposition plot
    % surfaces on the bottom row   
    figure(51); clf;
    for ii=1:2
        D2subplot(4+ii,2,3,[.1 0 0 0],[.1 .1 .06 .01]);
        w1 = reshape(pm(ii,:),[11 11]);
        mesh(gridx,gridy,w1); hold on; axis([0 1 0 1 -3 3]);
        xlabel('\alpha','fontSize',11);
        ylabel('\eta','fontSize',11);
        zlabel(['w_' num2str(ii) '(\eta, \alpha)'],'fontSize',11);
        %title(['PC ' num2str(ii)]);
    end
    % plot the mean on row 1 and then K1 and K2 on row 2
    ax = D2subplot(1,2,3,[0 0 0 0],[.15 .1 .1 .01]);
    pbaspect(ax, [3 2.5 1]); hold on;
     plot(pout.simData.orig.time,pout.simData.orig.y, 'Color',[.8,.8,.8]); hold on;
    plot(pout.simData.orig.time,pout.simData.orig.ymean, 'Color','b');
    xlabel('time', 'FontSize', 8);
    xticks([10 50 90]);
    ylabel('\phi_0', 'FontSize',11);
    axis([1 90 0 5]); grid on; box on;
    ax = D2subplot(2,2,3,[0 0 0 0],[.15 .1 .1 .01]);
    pbaspect(ax, [3 2.5 1]); hold on;
    plot(pout.simData.orig.time,pout.simData.Ksim(:,1)*pout.simData.orig.ysd); 
    ylabel('\phi_1', 'FontSize',11);
    xlabel('time', 'FontSize', 8);
    xticks([10 50 90]);
    axis([1 90 -0.5 1]);grid on; box on;
    ax = D2subplot(3,2,3,[0 0 0 0],[.15 .1 .1 .01]);
    pbaspect(ax, [3 2.5 1]); hold on;
    plot(pout.simData.orig.time,pout.simData.Ksim(:,2)*pout.simData.orig.ysd); 
    ylabel('\phi_2', 'FontSize',11);
    xlabel('time', 'FontSize', 8);
    xticks([10 50 90]);
    axis([1 90 -0.5 1]);grid on; box on;
    
    D2subplot(1,1,1,[0 0 0 0],[0 0 0 0]);
    set(gca,'Visible','off'); 
    text(0.17,0.95,'Mean');
    text(0.53,0.95,'PC1');
    text(0.85,0.95,'PC2');

end
simData = pout.simData;

%% now check the emulator to see how well it predicts holdouts
%% emulate the hold out sample
if doPlot(61)    
    [m pq] = size(pout.simData.x);
    ysimsho = simData.orig.yho;

 % read in the design
    designho01 = simData.orig.design01ho;
    nho = size(designho01,1);
    xho = repmat(.5,[nho 1]);
    xsamp = [xho designho01];
    
 % make predictions for each sample in pvals
    pred6=gPredict(xho,pvals,model,data,'mode','wpred','theta',designho01,'addResidVar',1);
    %pred6=gPred(xsamp,pvals,model,data,'wpred');
    % posterior draws of emulator
    wall = pred6.w;
    neta = size(ysimsho,1);
    % plot for each holdout
    [y iord] = sort(simData.orig.ymean);
    
    figure(63);clf;
    for kk=1:min([nho 56])
        etar = simData.Ksim*wall(:,:,kk)'.*simData.orig.ysd+...
              repmat(simData.orig.ymean,[1 nreal]);
        etaq = prctile(etar',[10 90]);
        D2subplot(kk,7,8,[.07 .1 .06 .01],[.0 0 0 0]);
        plot([ysimsho(:,kk) ysimsho(:,kk)]',etaq,'-k'); hold on;
        xlim([-.5 6.5]); ylim([-.5 6.5]); hold on; plot([-.5 6.5],[-.5 6.5],'y'); axis equal;
    end
    figure(64);clf;
    for kk=1:min([nho 56])
        etar = simData.Ksim*wall(:,:,kk)'.*simData.orig.ysd+...
              repmat(simData.orig.ymean,[1 nreal]);
        etaq = prctile(etar',[10 90]);
        D2subplot(kk,7,8,[.07 .1 .06 .01],[.0 0 0 0]);
        plot((exp([ysimsho(:,kk) ysimsho(:,kk)])-1)',exp(etaq)-1,'-k'); hold on;
        xlim([-.5 10^(6.5)-1]); ylim([-.5 10^(6.5)-1]); hold on; plot([-.5 10^6.5],[-.5 10^6.5],'y'); axis equal;
    end
    size(pvals)
    
    % make a plot of 3 different holdouts
    ihold = [3 7 39 95 120 199];
    figure(68);clf;
    for kk=1:length(ihold);
        etar = simData.Ksim*wall(:,:,ihold(kk))'.*simData.orig.ysd+...
              repmat(simData.orig.ymean,[1 nreal]);
        eta5 = prctile(etar',[50]);
        D2subplot(kk,2,3,[.17 .07 .06 .01],[.01 .01 .01 .01]);
        plot(simData.orig.time,eta5,'-k'); hold on;
        plot(simData.orig.time,ysimsho(:,ihold(kk)),':r');
        ylim([-1 6]);
    end

    
        
    
    % post mean predictions
    w6=squeeze(mean(pred6.w(:,:,:),1));
    etam6 = simData.Ksim*w6.*simData.orig.ysd+...
        repmat(simData.orig.ymean,[1 nho]);
    
    % compare emulations and sims
    res6 = ysimsho - etam6;
    res6Native = exp(ysimsho)-1 - exp(etam6)-1;
    rmse6 = sqrt(mean(res6.^2)');
    rmse6Native = sqrt(mean(res6Native.^2)');
    mean(rmse6)
    mean(rmse6Native)
    
    
    % now plot the post mean - ysimho residuals in fig 61
    figure(65);clf;
    D2subplot(1,2,2,[.17 .07 .06 .01],[.03 .05 .01 .01]);
    plot(simData.orig.time,res6);   
    D2subplot(2,2,2,[.17 .07 .06 .01],[.03 .05 .01 .01]);
    plot(simData.orig.time,res6Native);   
    
    
    h=rmse6;
    %return(eta6);

end

% now holdout predictions for all of the sims...
% we'll hold out two at a time...
   plabs = {'n','h','\sigma_8','\Omega_{CDM}','\Omega_b'};

if doPlot(7)    
    nhold = 2;
    m = size(pout.simData.x,1);
    holdoutfits = zeros(size(pout.simData.yStd));
    nrow = 6; ncol = 6;
    if newFig; figure(7); end
    for ii=1:(floor(m/nhold));
        iirow = floor((ii-1)/ncol)+1;
        iicol = rem(ii-1,ncol)+1;
        ihold = (ii-1)*nhold+(1:nhold);
        iok = ones([m 1]); iok(ihold) = 0; iok = (iok==1);
        obsData6 = pout.obsData;
        simData6 = pout.simData;
        simData6.x = simData6.x(iok,:);
        simData6.yStd = simData6.yStd(:,iok);
        p6 = setupModel(obsData6,simData6);
        model6 = p6.model;
        data6 = p6.data;
        xth6 = pout.simData.x(ihold,:);
        x6 = xth6(:,1:p);
        th6 = xth6(:,p+1:end);
        nsamp = length([pvals.lamOs]);
    
        pred6=gPred(x6,pvals,model6,data6,'wpred',th6);
        w6=squeeze(mean(pred6.w(:,:,:),1));
        eta6 = pout.simData.Ksim*w6.*pout.simData.orig.ysd+...
              repmat(pout.simData.orig.ymean,[1 nhold]);
        holdoutfits(:,ihold) = eta6;
        
        D2subplot(0+ii,10,10,[.09 .1 .06 .01],[0 0 0 0]);
        plot(pout.simData.orig.time,eta6,'g','LineWidth',1);
        hold on;
        plot(pout.simData.orig.time,pout.simData.orig.y(:,ihold),'.k',...
             'MarkerSize',4);
        set(gca,'Yticklabel',{'','',''},'Ytick',log([100 1e4 1e6]+1),'Xtick',[0 20 40],...
            'Xticklabel',{'','',''});
        axis([0 48 2 16]);
        if(iicol==1 & mod(iirow,2)==1)
            set(gca,'Yticklabel',{'10^2','10^4','10^6'},'Ytick',log([100 1e4 1e6]+1));
        end
        if(ii==13); ylabel('log_{10} Infected','FontSize',15); end
        if(iirow==6 & iicol==1) xlabel('weeks'); axis([0 48 2 16]); end
        if(iirow==5 & iicol==3) xlabel('weeks'); axis([0 48 2 16]); end
        if(iirow==5 & iicol==5) xlabel('weeks'); axis([0 48 2 16]); end
    end
    axes('Position',[0 0 1 1],'Visible','off');
    text(.52,.18,'simulations','FontSize',12);
    text(.52,.15,'response surface','FontSize',12,'Color',[0 1 0]);
    figure(17); clf;
    plot(pout.simData.orig.time,pout.simData.orig.y-holdoutfits,'.');
    axis([0 48 -5 5]); hold on;
    plot(pout.simData.orig.time,std(pout.simData.orig.y'-holdoutfits'),'c',...
        'linewidth',3);
    if(writeOut)
            fname = ['holdoutfits' suffix];
            save(fname,'holdoutfits','-ASCII');
    end

    % save holdoutfits holdoutfits -ASCII
end

% compare linear and non-linear posteriors 'compareCal' option
% cosmoPlots(pout,ip,8,'compareCal',pout03L)
if doPlot(8)
    if newFig; figure(8); clf; end
    t=[pout.pvals(min(pvec):max(pvec)).theta]';
    [nreal numt] = size(t);
    thetamax = pout.simData.orig.colmax;   %maximum values on native scale
    thetamin = pout.simData.orig.colmin;   %minimum values on native scale
    delta = thetamax - thetamin;
    tnative = zeros(size(t));
    for ii=1:numt;
        tnative(:,ii) = t(:,ii).*delta(ii) + thetamin(ii);
    end
    
    % now for linear piece
    tL=[poutL.pvals(min(pvec):end).theta]';
    [nrealL numt] = size(tL);
    tnativeL = zeros(size(tL));
    for ii=1:numt;
        tnativeL(:,ii) = tL(:,ii).*delta(ii) + thetamin(ii);
    end
    
    for ii=1:numt; 
        for jj=1:numt;
            if(ii==jj)
                D2subplot((ii-1)*numt+ii,numt,numt,[.08 .1 .08 .18],[0 0 0 0]);
                xout = linspace(0,1,100); 
                xoutnative = xout.*delta(ii)+thetamin(ii);
                sdx = 1/15;
                dx = dens1d(t(:,ii),sdx,xout);
                plot(xoutnative,dx,'b'); set(gca,'YtickLabel',''); hold on;
                xlim([thetamin(ii) thetamax(ii)]); ylim([0 6.5]);
                dxL = dens1d(tL(:,ii),sdx,xout);
                plot(xoutnative,dxL,'g'); set(gca,'YtickLabel',''); hold on;
            elseif(ii < jj)
                if(ii==1)|(jj==1)
                    ksd = .1;
                else
                    ksd = .08;
                end
                D2subplot((ii-1)*numt+jj,numt,numt,[.08 .1 .08 .18],[0 0 0 0]);
                mins=[thetamin(jj) thetamin(ii)];
                maxes=[thetamax(jj) thetamax(ii)];

                hpd2d01(t(:,[jj ii]),.09,'b',mins,maxes); hold on;
                hpd2d01(tL(:,[jj ii]),.09,'g',mins,maxes,0);
                plot(pout.obsData.orig.thetatrueNative(jj),...
                    pout.obsData.orig.thetatrueNative(ii),'.r',...
                   'MarkerSize',20);
            else
                D2subplot((ii-1)*numt+jj,numt,numt,[.08 .1 .08 .18],[0 0 0 0]);
                mins=[thetamin(jj) thetamin(ii)];
                maxes=[thetamax(jj) thetamax(ii)];

                hpd2d01(tL(:,[jj ii]),.09,'g',mins,maxes); hold on;
                hpd2d01(t(:,[jj ii]),.09,'b',mins,maxes,0);
                plot(pout.obsData.orig.thetatrueNative(jj),...
                    pout.obsData.orig.thetatrueNative(ii),'.r',...
                   'MarkerSize',20);
            end
            if(mod(ii,2)==0) & (jj==numt);
                set(gca,'YAxisLocation','right','XtickLabel','');
            elseif(mod(ii,2)==1) & (jj==1) & (ii == numt);
                set(gca,'YAxisLocation','left');
            elseif(mod(ii,2)==1) & (jj==1) & (ii > 1);
                set(gca,'YAxisLocation','left','XtickLabel','');
            elseif(ii==numt) & (mod(jj,2)==1);
                set(gca,'XAxisLocation','bottom','YtickLabel','');
            elseif(ii==1) & (mod(jj,2)==0);
                set(gca,'XAxisLocation','top','YtickLabel','');
            else
                set(gca,'XtickLabel','','YtickLabel','','Clim',[0 .01]);  
            end
            %set(gca,'XtickLabel','','YtickLabel','','Clim',[0 .01]);
                
        end
    end
    % now put in the labels
        xlabs = {'n','h','\sigma_8','\Omega_M','\Omega_B'};
    D2subplot(1,1,1,[.08 .1 .08 .18],[0 0 0 0]);
    set(gca,'Visible','off');
    for ii=1:numt
        posx=(.71 + ii-1)/numt;
        posy=(.75 + numt-ii)/numt;
        text(posx,posy,xlabs{ii},'FontSize',14);
    end
    drawnow
end    
% figure(8); print -deps2c cosmopost03NL.eps;

% show the non-linear posterior 'compareCal' option
% cosmoPlots(pout,ip,9,'compareCal')
if doPlot(9)
    if newFig; figure(9); clf; end
    t=[pout.pvals(min(pvec):max(pvec)).theta]';
    [nreal numt] = size(t);
    thetamax = pout.simData.orig.colmax;   %maximum values on native scale
    thetamin = pout.simData.orig.colmin;   %minimum values on native scale
    delta = thetamax - thetamin;
    tnative = zeros(size(t));
    for ii=1:numt;
        tnative(:,ii) = t(:,ii).*delta(ii) + thetamin(ii);
    end
        
    for ii=1:numt; 
        for jj=1:numt;
            if(ii==jj)
                D2subplot((ii-1)*numt+ii,numt,numt,[.08 .1 .08 .18],[0 0 0 0]);
                xout = linspace(0,1,100); 
                xoutnative = xout.*delta(ii)+thetamin(ii);
                sdx = 1/15;
                dx = dens1d(t(:,ii),sdx,xout);
                plot(xoutnative,dx,'b'); set(gca,'YtickLabel',''); hold on;
                xlim([thetamin(ii) thetamax(ii)]); ylim([0 6.5]);
            elseif(ii < jj)
                if(ii==1)|(jj==1)
                    ksd = .1;
                else
                    ksd = .08;
                end
                D2subplot((ii-1)*numt+jj,numt,numt,[.08 .1 .08 .18],[0 0 0 0]);
                mins=[thetamin(jj) thetamin(ii)];
                maxes=[thetamax(jj) thetamax(ii)];

                hpd2d01(t(:,[jj ii]),.09,'b',mins,maxes); hold on;
                plot(pout.obsData.orig.thetatrueNative(jj),...
                    pout.obsData.orig.thetatrueNative(ii),'.r',...
                   'MarkerSize',20);
            else
                D2subplot((ii-1)*numt+jj,numt,numt,[.08 .1 .08 .18],[0 0 0 0]);
                mins=[thetamin(jj) thetamin(ii)];
                maxes=[thetamax(jj) thetamax(ii)];

                hpd2d01(t(:,[jj ii]),.09,'b',mins,maxes); hold on;
                plot(pout.obsData.orig.thetatrueNative(jj),...
                    pout.obsData.orig.thetatrueNative(ii),'.r',...
                   'MarkerSize',20);
            end
            if(mod(ii,2)==0) & (jj==numt);
                set(gca,'YAxisLocation','right','XtickLabel','');
            elseif(mod(ii,2)==1) & (jj==1) & (ii == numt);
                set(gca,'YAxisLocation','left');
            elseif(mod(ii,2)==1) & (jj==1) & (ii > 1);
                set(gca,'YAxisLocation','left','XtickLabel','');
            elseif(ii==numt) & (mod(jj,2)==1);
                set(gca,'XAxisLocation','bottom','YtickLabel','');
            elseif(ii==1) & (mod(jj,2)==0);
                set(gca,'XAxisLocation','top','YtickLabel','');
            else
                set(gca,'XtickLabel','','YtickLabel','','Clim',[0 .01]);  
            end
            %set(gca,'XtickLabel','','YtickLabel','','Clim',[0 .01]);
                
        end
    end
    % now put in the labels
        xlabs = {'n','h','\sigma_8','\Omega_M','\Omega_B'};
    D2subplot(1,1,1,[.08 .1 .08 .18],[0 0 0 0]);
    set(gca,'Visible','off');
    for ii=1:numt
        posx=(.71 + ii-1)/numt;
        posy=(.75 + numt-ii)/numt;
        text(posx,posy,xlabs{ii},'FontSize',14);
    end
    drawnow
end    
% figure(9); print -deps2c cosmopost03N.eps;

%% now check the emulator to see how well it 64 holdouts
if doPlot(10)    
    if newFig; figure(10); end
    clf
    logysimho = textread('CosmoDat/logf.64.sims');
    logksimho = textread('CosmoDat/logk.64.sims');
    [m pq] = size(pout.simData.x);
    holddesign = textread('CosmoDat/native5x64');
    holddesign = to01(holddesign);
    [mh pqh] = size(holddesign);
    holddesign = [ones([mh 1])*.5 holddesign];
    xh = holddesign(:,1:p);
    th = holddesign(:,p+1:end);
    nsamp = length(pvec);
    
    
    ihold = [5 25 m]; % the sims being held out
    ihold = [5 6 7];
    iok = ones([m 1]); iok(ihold) = 0; iok = (iok==1);
    obsData6 = pout.obsData;
    simData6 = pout.simData;
    simData6.x = simData6.x(iok,:);
    simData6.yStd = simData6.yStd(:,iok);
    p6 = setupModel(obsData6,simData6);
    model6 = p6.model;
    data6 = p6.data;
    xth6 = pout.simData.x(ihold,:);
    x6 = xth6(:,1:p);
    th6 = xth6(:,p+1:end);
    nsamp = length([pvals.lamOs]);
    
    pred6=gPred(xh,pvals,model,data,'wpred',th);
 %   pred=gPred(x6,pvals,model,data,'wpred',th6);

    % now plot the predictions
    nrow=8; ncol=8;
    for ii=1:mh
        D2subplot(0+ii,8,8,[.09 .1 .06 .01],[0 0 0 0]);
        eta6 = pout.simData.Ksim*pred6.w(:,:,ii)'.*pout.simData.orig.ysd+...
              repmat(pout.simData.orig.ymean,[1 nsamp]);
        etabounds6 = prctile(eta6,[2.5 97.5],2);
        plot(pout.simData.orig.logk,etabounds6,'g','LineWidth',1);
        hold on;
        plot(logksimho,logysimho(:,ii),'.k',...
             'MarkerSize',2);
         
        iirow = floor((ii-1)/ncol)+1;
        iicol = rem(ii-1,ncol)+1;
 
        set(gca,'Yticklabel',' | | ','Ytick',[2 4 6],'Xtick',[-2 -1],...
            'Xticklabel',' | ');
        axis([-3 0 2 6]);
        if(iicol==1 & mod(iirow,2)==1)
            set(gca,'Yticklabel','2|4|6','Ytick',[2 4 6]);
        end
        %if(ii==13); ylabel('log spectrum','FontSize',15); end
        if(iirow==8 & iicol==1) xlabel('log k'); axis([-3 0 2 6]); set(gca,'Xticklabel','-3|-2|-1','Xtick',[-3 -2 -1]); end
        if(iirow==8 & iicol==3) xlabel('log k'); axis([-3 0 2 6]); set(gca,'Xticklabel','-3|-2|-1','Xtick',[-3 -2 -1]);end
        if(iirow==8 & iicol==5) xlabel('log k'); axis([-3 0 2 6]); set(gca,'Xticklabel','-3|-2|-1','Xtick',[-3 -2 -1]);end
        if(iirow==8 & iicol==7) xlabel('log k'); axis([-3 0 2 6]); set(gca,'Xticklabel','-3|-2|-1','Xtick',[-3 -2 -1]);end
         
             
        if(ii==25); ylabel('log spectrum','FontSize',14); end
        %xlabel('log k'); axis([-3 0 2 6]);
        %title(['holdout pred run ' num2str(ihold(ii))],'Fontsize',11);
    end

end

if doPlot(20)
    if newFig; figure(20); clf; end;
    % plot the sim and ksim
    ax = D2subplot(1,1,3,[.1 .1 .08 .16],[0.1 0 0.01 0]);
    pbaspect(ax, [3 2.5 1]); hold on;
    plot(simData.orig.time,simData.orig.y,'Color',[.8,.8,.8]);
    axis([0 57 0 15]); box on;
    xlabel('weeks')
    ylabel('log(cumulative infected)', 'FontSize', 10)
    ax = D2subplot(3,1,3,[.1 .1 .08 .16],[0.1 0 0.01 0]);
    pbaspect(ax, [3 2.5 1]);hold on;
    plot(simData.orig.time,simData.Ksim);
    xlabel('weeks')
    ylabel('', 'FontSize', 10)
    axis([0 57 -1 0.5]); box on;
end

return

% if ~exist('neddfigs'); mkdir('neddfigs'); end
% cd neddfigs
figure(10); print -depsc holdout64.eps
figure(1); print -deps2c rhobox.eps
figure(2); print -deps2c thetascatter.eps
figure(3); print -deps2c sens.eps
figure(4); print -deps2c preds.eps
% cd ..
% 
% return
