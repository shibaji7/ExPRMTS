function v = ovaltine(pout,varargin)
% Takes a pout and turns into the appropraite params.dat for use with the
% emu code.

% assume burn in is 500 unless specified in the call which is of the
% form ovaltine(pout,'burn',1000);
burn = 500;
outfile = 'params.dat';
for ii=1:2:length(varargin)
  switch varargin{ii}
  case 'burn'
    burn=varargin{ii+1};
  case 'outfile';
    outfile = varargin{ii+1};
    writeOut=1;
  otherwise
    error('invalid extended argument passed to cosmoPlots');
  end
end

% now remove the burnin
xx = pout.pvals(burn:end);
pout.pvals = xx;



% flength m neta p peta X(:,1)...X(:,p) Xmin Xrange W(1,:)...W(peta,:) 
%      K(:,1)...K(:,peta) beta(1,:)...beta(peta,:) lamws lamz mean sd

v = [pout.model.m, size(pout.simData.yStd,1), pout.model.p+pout.model.q, pout.model.pu];
% Add in the design
for i=1:v(3)
    v = [v, pout.simData.x(:,i)'];
end

% Add in the min and range for each X
v = [v, 0, pout.simData.orig.colmin, 1, pout.simData.orig.colmax-pout.simData.orig.colmin];

% Add in the W
for i=1:v(4)
    v = [v, pout.data.w(:,i)'];
end

% Add in K
for i=1:v(4)
    v = [v, pout.simData.Ksim(:,i)'];
end

% Add in the parameter estimates.  We'll be using the medians here.
v = [v, median([pout.pvals.betaU],2)'];

% lamws
v = [v, median([pout.pvals.lamWs],2)'];

% lamz
v = [v, median([pout.pvals.lamUz],2)'];

% mean
for i=1:4
    v = [v, pout.simData.orig.ymean(:,i)'];
end

v = [v, kron(pout.simData.orig.ysd, ones(1,1099))];

v = [length(v), v];
save(outfile,'v', '-ASCII');