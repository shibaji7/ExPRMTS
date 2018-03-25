function pvalsOut = medianPvals(pvals,varargin)
% Takes a pvals array and turns into a single set of parameter values by
% computing the median for each element in pvals.  This function will
% automatically chop off the first 500 values as burn-in.  The user can
% specify another burn-in value

% assume burn in is 500 unless specified in the call which is of the
% form medianPvals(pout,'burn',1000);
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

% now remove the burnin (burn:end)
pvals = pvals(burn:end);
pvalsOut = pvals(end);

% Compute point estimates for each parameter.
pvalsOut.theta = median([pvals.theta],2);
pvalsOut.betaV = median([pvals.betaV],2);
pvalsOut.betaU = median([pvals.betaU],2);
pvalsOut.lamVz = median([pvals.lamVz],2);
pvalsOut.lamUz = median([pvals.lamUz],2);
pvalsOut.lamWs = median([pvals.lamWs],2);
pvalsOut.lamWOs = median([pvals.lamWOs],2);
pvalsOut.lamOs = median([pvals.lamOs],2);

end