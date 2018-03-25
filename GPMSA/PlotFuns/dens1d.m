function dens=dens1d(x,sd,xout);
%estimates a 1d density from the sample x using a normal
%kernel smooth of x.  The normal kernel has an sd of sd and
%the output is computed at xout
m = length(xout);
dens = zeros([m 1]);
for i=1:m
    dens(i) = mean(dnorm(x,xout(i),sd));
end

