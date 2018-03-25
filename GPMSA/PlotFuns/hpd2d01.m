 %%% a function to give bivariate hpd contours
 function out = hpd2d01(x,ksd,col,minx,maxx,shade,ngrid,prob);
  % x is a n x 2 matrix; col1 is the first component; col 2 is the second.
   if ~exist('ksd'); ksd = .08; end;
   if ~(exist('ngrid')==1); ngrid=round(1.8/ksd); end;
   if ~exist('prob'); prob=[.9]; end;
   if ~exist('shade'); shade=1; end;
   if ~exist('col'); col='b'; end;
   if ~exist('minx'); minx=[0 0]; end;
   if ~exist('maxx'); maxx=[1 1]; end;
   gridvals = 0:(1/(ngrid-1)):1;
   [x1 x2] = meshgrid(gridvals,gridvals);
   x1v = x1(:); x2v = x2(:);
   ng = length(gridvals);
   [n m] = size(x);
   nim = length(x1v);
   dens = zeros(size(x1v));
   
   % now compute density for each image pixel;
   
   for i=1:nim
       f = normpdf(x(:,1),x1v(i),ksd).*normpdf(x(:,2),x2v(i),ksd);
       dens(i) = sum(f);
   end
  
   % normalize dens
   dens = dens/sum(dens);
   hlevels = zeros(size(prob));
   nprob = length(prob);
   for j=1:nprob
       hlevels(j) = fzero(@(x) getp(x,dens)-prob(j),[0 max(dens)]);
   end
   if nprob==1;
       hlevels = [hlevels 1.00000001*hlevels];
   end
%   plot(x(:,1),x(:,2),'y.'); hold on;
   %plot(x(:,1),x(:,2),'.'); 
   colormap(repmat([.9:-.02:.3]',[1 3]));
   x1im = x1.*(maxx(1)-minx(1))+minx(1);
   x2im = x2.*(maxx(2)-minx(2))+minx(2);
   densim = reshape(dens,size(x1));

   x1yr = x1im(:,ngrid:1);
   x2yr = x2im(:,ngrid:1);
   densimyr = densim(:,ngrid:1);

   linest = ':'; lwd=1.2;
   if(shade);
      h=imagesc(x1im(:),x2im(:),densim); hold on;
      set(gca,'YDir','normal');
      linest = '-'; lwd=1.5;
   end
   
   %imagesc(x1im(:),x2im(:),reshape(dens,size(x1))); hold on;
   %colormap('default');
   contour(reshape(x1im,size(x1)),reshape(x2im,size(x1)),...
           reshape(dens,size(x1)),hlevels,'LineWidth',lwd,'Color',col,...
           'LineStyle',linest); 
   contour(reshape(x1im,size(x1)),reshape(x2im,size(x1)),...
           densim,hlevels,'LineWidth',lwd,'Color',col,...
           'LineStyle',linest); 
   %hold off;
   
   % function to get probability of a given level h
%      function pout = getp(h,d);
%       %  d = dens; p = .9;
%       % function to find the density height corresponding
%       % to the p-hpd region
%       iabove = (d >= h);
%       pout = sum(d(iabove));
      
      
