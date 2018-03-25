     function pout = getp(h,d);
      %  d = dens; p = .9;
      % function to find the density height corresponding
      % to the p-hpd region
      iabove = (d >= h);
      pout = sum(d(iabove));
