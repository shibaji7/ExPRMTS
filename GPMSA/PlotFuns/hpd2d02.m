 %%% a function to give bivariate hpd contours
 function out = hpd2d02(x,mins,maxes,col,ksd,ngrid,prob,shade);
  % x is a n x 2 matrix; col1 is the first component; col 2 is the second.
   if ~exist('col'); col = 'b'; end;
   if ~exist('ksd'); ksd = .08; end;
   if ~exist('ngrid'); ngrid=round(1.8/ksd); end;
   if ~exist('prob'); prob=[.9]; end;
   if ~exist('shade'); shade=0; end;
   gridvals = 0:(1/(ngrid-1)):1;
   deltas=maxes-mins;
   grid1=gridvals*deltas(1)+mins(1);
   grid2=gridvals*deltas(2)+mins(2);
   [x1 x2] = meshgrid(grid1,grid2);
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
%    plot(x(:,1),x(:,2),'y.'); hold on;
   %plot(x(:,1),x(:,2),'.'); 
   colormap(repmat([.9:-.02:.3]',[1 3]));
   imagesc(x1(:),x2(:),reshape(dens,size(x1))); hold on;
   %colormap('default');
   contour(x1,x2,reshape(dens,size(x1)),hlevels,'LineWidth',1.0,'Color',col); 
   hold off;
   
   % function to get probability of a given level h
%      function pout = getp(h,d);
%       %  d = dens; p = .9;
%       % function to find the density height corresponding
%       % to the p-hpd region
%       iabove = (d >= h);
%       pout = sum(d(iabove));
      
      
