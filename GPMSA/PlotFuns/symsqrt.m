function y = symsqrt(x);
    % a symmetric sqrt transform that returns sign(x)*sqrt(|x|+.25)
    y = sign(x).*sqrt(abs(x)+.25);
end
