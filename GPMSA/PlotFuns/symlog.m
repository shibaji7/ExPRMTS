function y = symlog(x);
    % a symmetric log transform that returns sign(x)*log(|x|+1)
    y = sign(x).*log(abs(x)+1.0);
end