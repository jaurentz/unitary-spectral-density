function [rho,info]=usd_poisson_density(U,theta,s,r,seed,tol,restart,maxit)
%USD_POISSON_DENSITY Direct stochastic Poisson-resolvent spectral filtering.
%   rho_r(theta)=Re[tr((U+zI)/(U-zI))]/(2*pi*n), z=r*exp(1i*theta).
%   Uses only U(X,'notransp') through GMRES shifted solves, not
%   Fourier moments, quadrature, a matrix class, or an explicit U.
%   Small r is smoother; r close to one requires harder shifted solves.
%   When RESTART >= n, MATLAB interprets MAXIT as total iterations,
%   not as a number of restart cycles.
if nargin<5,seed=17;end
if nargin<6,tol=1e-9;end
sz=U([],'size');n=sz(1);
assert(numel(sz)==2 && sz(2)==n,'U must be square.');
validateattributes(s,{'numeric'},{'scalar','integer','positive'});
assert(isscalar(r) && r>0 && r<1,'Require 0<r<1.');
if nargin<7 || isempty(restart),restart=min(40,n);end
if nargin<8 || isempty(maxit)
    % MATLAB gmres counts MAXIT as *total* iterations for unrestarted
    % GMRES (restart == n), but as restart cycles otherwise.
    if restart>=n
        maxit=n;
    else
        maxit=ceil(n/restart)+5;
    end
end
state=rng;restore=onCleanup(@()rng(state)); %#ok<NASGU>
rng(seed,'twister');
probes=zeros(n,s);
for j=1:s,probes(:,j)=usd_fourphase(n);end
theta=theta(:);rho=zeros(size(theta));calls=0;
maxrel=0;
for k=1:numel(theta)
    z=r*exp(1i*theta(k));
    acc=0;
    for j=1:s
        v=probes(:,j);
        [y,flag,relres]=gmres(@(x)shifted(x,z),v, ...
                             restart,tol,maxit);
        if flag~=0
            error('usd_poisson_density:GMRES', ...
                'Shifted solve did not converge (theta index %d, flag %d).',k,flag);
        end
        maxrel=max(maxrel,relres);
        acc=acc+v'*y/s;
    end
    rho(k)=real(1+2*z*acc)/(2*pi);
end
info=struct('matvecs',calls,'probes',s,'radius',r, ...
    'max_relative_residual',maxrel,'solver_tolerance',tol);
    function y=shifted(x,z)
        calls=calls+1;
        y=U(x,'notransp')-z*x;
    end
end
