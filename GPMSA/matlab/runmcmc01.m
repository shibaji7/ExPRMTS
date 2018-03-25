addpath('/home/afadikar/Documents/shibaji/git/ExPRMTS/GPMSA/gpmsa_2015/matlab')
addpath('/home/afadikar/Documents/shibaji/git/ExPRMTS/GPMSA/matlab')
addpath('/home/afadikar/Documents/shibaji/git/ExPRMTS/GPMSA/PlotFuns')

ebdat = euvac_eg(1, 'ott');  
simData = ebdat.simData;
obsData = ebdat.obsData;
params = setupModel(obsData,simData);
params.priors.lamOs.bLower = .2;
params.priors.lamWs.params = [1 .0001];
params.model.lamVz=200;
params.model.lamOs=1.0;
params.priors.lamOs.params = [10 10];
params.priors.lamVz.params = [1 1.0000e-04];
%params.priors.lamVz.params = [10000000 1000];


nburn = 100; 
nlev = 9;
pout = stepsize(params,nburn,nlev);

pout=gpmmcmc(pout,10000);
save('pout_10000_ott.mat', 'pout');

showPvals(pout.pvals);

ip = round(linspace(900,10900,50));

% load previously saved pout object for plotting
load('pout_10000_ott.mat')

plots(pout, ip, 1);
plots(pout, ip, 2);
plots(pout, ip, 3);
plots(pout, ip, 5);
plots(pout, ip, 4);