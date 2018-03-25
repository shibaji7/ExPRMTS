function stairplot(x,y,colorvec)
 xadj = mean(diff(x))/2;
 xm = [x'-xadj x'+xadj]';
 xx = xm(:);
 ym = [y' y']';
 yy = ym(:);
 plot(xx,yy,'color',colorvec);
end