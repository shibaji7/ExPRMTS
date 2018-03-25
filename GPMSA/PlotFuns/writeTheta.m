function theta = writeTheta(pout,varargin)
% Takes a pout and writes out the mcmc sample of thetas into the file
% theta.out, or whatever specified by 'outfile','theta.out'.  Also the call
% can specify which samples to write out using 'isample',1:100, or
% whatever.

% assume burn in is 20% unless specified in the call which is of the
% form writeTheta(pout,'isample',500:5:10000);

outfile = 'theta.out';
nmcmc = length(pout.pvals);
isample = round(linspace(.2*nmcmc,nmcmc,200));
for ii=1:2:length(varargin)
  switch varargin{ii}
  case 'isample'
    isample=varargin{ii+1};
  case 'outfile';
    outfile = varargin{ii+1};
    writeOut=1;
  otherwise
    error('invalid extended argument passed to writeTheta()');
  end
end

% now remove the burnin
theta01 = [pout.pvals(isample).theta]';
q = pout.model.q;
thmin = pout.simData.orig.colmin;
thrange = pout.simData.orig.colmax - thmin;
theta = 0*theta01;
for k=1:q
    theta(:,k) = theta01(:,k)*thrange(k) + thmin(k);
end


save(outfile,'theta', '-ASCII');
end
