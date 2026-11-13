function v=usd_fourphase(n)
%USD_FOURPHASE Isotropic random unit vector (normalized complex Rademacher).
%   E[v*v']=I/n; every entry has modulus 1/sqrt(n).
roots=[1;1i;-1;-1i];
v=roots(randi(4,n,1))/sqrt(n);
end
