function x01 = to01(x)
% normalizes the columns of x so that all of its values lie between
% 0 and 1.  (x - min(x))/(max(x)-min(x))
%
[n m] = size(x);

x01 = 0*x;
for k = 1:m
    x01(:,k) = (x(:,k) - min(x(:,k)))/(max(x(:,k))-min(x(:,k)));
end

