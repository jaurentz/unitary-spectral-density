function [rho,info]=usd_moment_density(U,theta,K,s,r,seed)
%USD_MOMENT_DENSITY Stochastic Fourier moments and truncated Poisson series.
%   U(X,'notransp') applies a unitary, U([],'size') reports dimensions.
%   Returns a Poisson-smoothed trace density; truncation may affect positivity.
%   E[z*z']=I for unnormalized four-phase probes; normalized probes are
%   used here so that E[v'*U^k*v]=trace(U^k)/n.
if nargin<6,seed=17;end
sz=U([],'size');n=sz(1);
assert(numel(sz)==2 && sz(2)==n && n>=1,'U must be square.');
validateattributes(K,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(s,{'numeric'},{'scalar','integer','positive'});
assert(isscalar(r) && r>0 && r<1,'Require 0<r<1.');
state=rng;restore=onCleanup(@()rng(state)); %#ok<NASGU>
rng(seed,'twister');
mu=zeros(K,1);
for j=1:s
    v=usd_fourphase(n);
    x=v;
    for k=1:K
        x=U(x,'notransp');
        mu(k)=mu(k)+v'*x/s;
    end
end
theta=theta(:);
rho=ones(size(theta))/(2*pi);
for k=1:K
    rho=rho+real((r^k*mu(k)/pi)*exp(-1i*k*theta));
end
info=struct('moments',mu,'matvecs',K*s,'probes',s, ...
    'degree',K,'radius',r);
end
