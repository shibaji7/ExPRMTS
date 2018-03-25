function params=ebola02eg(doPlot, loc);

if ~exist('doPlot'); doPlot=0; else; doPlot=1; end

ysim = textread(strcat('../R/sim_', loc,'.txt'));
[neta m] = size(ysim);
time = [1:neta]';

designNative = textread(strcat('../R/design_', loc,'.txt'));
dmat1 = repmat(designNative,[5 1]);
design01 = to01(designNative);

if(doPlot)
 figure(1); clf;
 plotmatrix(design01);
end
if(doPlot)
 figure(2); clf;
 plot(time,ysim); 
end


% read in the observations
obsdat = textread(strcat('../R/obs_', loc,'.txt'));
ydat = obsdat;
Sigy = 0.2 * eye(neta);

n=1;
yobs(1).y = ydat; yobs(1).time = time; 
yobs(1).Sigy = Sigy;

% Standardize the data
ysimmean=repmat(mean(ysim,2),size(time));
ysimmean=mean(ysim,2);
ysimsd=std(ysim(:));
ysimStd=(ysim-repmat(ysimmean,[1 m])) /ysimsd;
% interpolate to data grid and standardize experimental data
for ii=1:n
  yobs(ii).ymean=interp1(time,ysimmean,yobs(ii).time);
  yobs(ii).yStd=(yobs(ii).y-yobs(ii).ymean)/ysimsd;
end
if(doPlot)
 figure(3); clf;
 plot(time,ysimStd); hold on;
 plot(yobs.time,yobs.yStd,'ko')
end


[U,S,V]=svd(ysimStd,0);
lam=diag(S).^2/sum(diag(S).^2); lam=cumsum(lam);
pu = 3;
Ksim=U(:,1:pu)*S(1:pu,1:pu)./sqrt(m);
if(doPlot)
    figure(4);
    plot(time,Ksim);
%figure(2); print -deps2c cipdssMortpc.eps;
end

% interpolate K onto data grids
for ii=1:n
  yobs(ii).Kobs=zeros(length(yobs(ii).yStd),pu);
  for jj=1:pu
    yobs(ii).Kobs(:,jj)=interp1(time,Ksim(:,jj),yobs(ii).time);
  end
end

% D basis - build the discrepancy basis
% lay it out, and record decomposition on sim and data grids
% Kernel centers and widths
  Dgrid=-2:5:98; Dwidth=5; 
  pv=length(Dgrid);
% Compute the kernel function map, for each kernel
  Dsim=zeros(size(ysimStd,1),pv);
  for ii=1:n; yobs(ii).Dobs=zeros(length(yobs(ii).yStd),pv); end
  for jj=1:pv
    % first the obs
      for ii=1:n
        yobs(ii).Dobs(:,jj)=normpdf(yobs(ii).time,Dgrid(jj),Dwidth);
      end
    % now the sim
      Dsim(:,jj)=normpdf(time,Dgrid(jj),Dwidth);
  end
% normalize the D maps
  Dmax=max(max(Dsim*Dsim'));
  Dsim=Dsim/sqrt(Dmax);
  for ii=1:n; yobs(ii).Dobs=yobs(ii).Dobs/sqrt(Dmax); end


% record the data into the obsData and simData structures.
% First simData
% required fields
simData.x   = [repmat(.5, [m 1]),design01];
simData.x = [repmat(.5, [m 1]),design01(:,1)];
%simData.x   = [.5*ones(size(dmatQ)) design01];
simData.yStd=ysimStd;
simData.Ksim=Ksim;
% extra fields: original data and transform stuff
simData.orig.y=ysim;
simData.orig.ymean=ysimmean;
simData.orig.ysd=ysimsd;
simData.orig.Dsim=Dsim;
simData.orig.time=time;
simData.orig.designNative = designNative;

% obsData
% set the dummy x values to .5; we're treating all inputs as
% parameters to be calibrated.
x = [.5];
for ii=1:n
% required fields
obsData(ii).x   =x(ii);
obsData(ii).yStd=yobs(ii).yStd;
obsData(ii).Kobs=yobs(ii).Kobs;
obsData(ii).Dobs=yobs(ii).Dobs;
obsData(ii).Sigy=yobs(ii).Sigy./(ysimsd.^2); 
% extra fields
obsData(ii).orig.y=yobs(ii).y;
obsData(ii).orig.ymean=yobs(ii).ymean;
obsData(ii).orig.time =yobs(ii).time;
%         obsData(ii).orig.thetatrue = design01ho(itrue,:);
end

% pack up and leave
params.simData=simData;
params.obsData=obsData;









end