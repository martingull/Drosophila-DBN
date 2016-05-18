function zz = ZZ(dd,nn,rr,vv,CC,XX);

zz = - nn*dd/2*log(pi) ...
     - dd/2*log(rr) ...
     - vv*sum(log(diag(cholupdate(CC,XX/sqrt(rr),'-')))) ...
     + sum(gammaln((vv-(0:dd-1))/2));