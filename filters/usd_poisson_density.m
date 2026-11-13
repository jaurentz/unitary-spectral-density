function [rho,info]=usd_poisson_density(U,theta,s,r,seed,tol,restart,maxit,opts)
%USD_POISSON_DENSITY Poisson-resolvent trace density by shifted solves.
%   rho_r(theta) = Re[trace((U+zI)/(U-zI))]/(2*pi*n),
%   z=r*exp(1i*theta), 0<r<1. Only the Anymatrix operator action
%   U(x,'notransp') and U([],'size') are required.
%
%   The legacy 8-positional-argument call defaults to MATLAB GMRES.
%   Ninth argument OPTS selects experimental projected alternatives:
%       OPTS.solver = 'gmres' (default) or 'sketch'
%       OPTS.sketch_order = Krylov subspace dimension (default min(n,32))
%       OPTS.sketch_dimension = projection rows (default 2*order)
%       OPTS.sketch_type = 'gaussian' (default) or 'countsketch'
%       OPTS.sketch_fallback = false (default): if true, retry with
%          full GMRES when the TRUE relative residual is too large
%       OPTS.sketch_residual_tolerance = 1e-6 (for fallback)
%       OPTS.probe_type / probe_vectors as in USD_PROBE_VECTORS.
%
%   The sketch method forms a full-reorthogonalized forward Krylov
%   basis once per probe and solves sketched residual least squares
%      min_c ||S*((U-zI)*Q*c-v)||_2
%   for every theta, reusing U*Q. This is a prototype of sketched
%   Krylov projection, NOT an optimized low-memory sGMRES: storage
%   remains O(n*order), and sketch error is separately diagnosed.
%   CountSketch uses sparse storage O(n+sketch_dimension).
%   Sketches are independent of the trace probes. No convergence
%   guarantee is implied without a residual check or fallback.
if nargin<5 || isempty(seed),seed=17;end
if nargin<6 || isempty(tol),tol=1e-9;end
if nargin<9 || isempty(opts),opts=struct;end
sz=U([],'size');n=sz(1);
assert(numel(sz)==2 && sz(2)==n && n>=1,'U must be square.');
validateattributes(s,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
assert(isscalar(r) && isreal(r) && r>0 && r<1, ...
    'Poisson radius must satisfy 0<r<1.');
if nargin<7 || isempty(restart),restart=min(40,n);end
if nargin<8 || isempty(maxit)
    if restart>=n
        maxit=n;
    else
        maxit=ceil(n/restart)+5;
    end
end
if ~isfield(opts,'solver'),opts.solver='gmres';end
solver=lower(char(opts.solver));
if strcmp(solver,'sketched'),solver='sketch';end
assert(any(strcmp(solver,{'gmres','sketch'})), ...
    'solver must be gmres or sketch.');
opts.seed=seed;
[V,kind]=usd_probe_vectors(n,s,opts);
p=size(V,2);
theta=theta(:);
samples=zeros(numel(theta),p);
residuals=zeros(numel(theta),p);
calls=0;fallbacks=0;
orders=zeros(p,1);
if strcmp(solver,'sketch')
    if ~isfield(opts,'sketch_order') || isempty(opts.sketch_order)
        opts.sketch_order=min(n,32);
    end
    m=opts.sketch_order;
    validateattributes(m,{'numeric'}, ...
        {'scalar','integer','positive','<=',n});
    if ~isfield(opts,'sketch_dimension') || isempty(opts.sketch_dimension)
        opts.sketch_dimension=2*m;
    end
    q=opts.sketch_dimension;
    validateattributes(q,{'numeric'}, ...
        {'scalar','integer','>=',m,'finite'});
    if ~isfield(opts,'sketch_type'),opts.sketch_type='gaussian';end
    if ~isfield(opts,'sketch_fallback'),opts.sketch_fallback=false;end
    if ~isfield(opts,'sketch_residual_tolerance')
        opts.sketch_residual_tolerance=1e-6;
    end
    validateattributes(opts.sketch_residual_tolerance,{'numeric'}, ...
        {'scalar','real','positive','finite'});
    % Independent reproducible sketch stream, without altering caller RNG.
    state=rng;restore=onCleanup(@()rng(state)); %#ok<NASGU>
    rng(seed+7919,'twister');
else
    m=0;q=0;
end
for j=1:p
    % Some matrix-free operators require dense RHS (e.g. FFT kernels).
    % Keep basis probes sparse in storage; materialize one column only.
    v=full(V(:,j));
    if strcmp(solver,'sketch')
        [Q,AQ]=krylov_data(v,m);
        orders(j)=size(Q,2);
        S=make_sketch(q,opts.sketch_type);
        SQ=S*Q;SA=S*AQ;sv=S*v;
    end
    for k=1:numel(theta)
        z=r*exp(1i*theta(k));
        if strcmp(solver,'gmres')
            [y,relres]=full_solve(v,z,k);
            scalar=v'*y;
        else
            T=SA-z*SQ;
            coeff=T\sv;
            residual=norm(AQ*coeff-z*(Q*coeff)-v)/norm(v);
            if ~isfinite(residual),residual=inf;end
            if opts.sketch_fallback && ...
                    residual>opts.sketch_residual_tolerance
                [y,relres]=full_solve(v,z,k);
                scalar=v'*y;
                fallbacks=fallbacks+1;
            else
                scalar=(v'*Q)*coeff;
                relres=residual;
            end
        end
        residuals(k,j)=relres;
        samples(k,j)=real(1+2*z*scalar)/(2*pi);
    end
end
rho=mean(samples,2);
if p>1 && ~strcmp(kind,'basis')
    stderr=std(samples,0,2)/sqrt(p);
else
    stderr=zeros(size(rho));
end
info=struct('matvecs',calls,'probes',p,'radius',r, ...
    'max_relative_residual',max(residuals(:)), ...
    'relative_residuals',residuals,'solver_tolerance',tol, ...
    'solver',solver,'probe_type',kind,'standard_error',stderr, ...
    'sketch_order',m,'effective_orders',orders, ...
    'sketch_dimension',q,'fallback_solves',fallbacks);
if strcmp(solver,'sketch')
    info.sketch_type=lower(char(opts.sketch_type));
else
    info.sketch_type='none';
end

    function [y,relres]=full_solve(v,z,k)
        [y,flag,relres]=gmres(@(x)shifted(x,z),v, ...
            restart,tol,maxit);
        if flag~=0
            error('usd_poisson_density:GMRES', ...
                'Shifted solve failed at theta index %d (flag %d).',k,flag);
        end
    end

    function w=shifted(x,z)
        calls=calls+1;
        w=U(x,'notransp')-z*x;
    end

    function [Q,AQ]=krylov_data(v,order)
        Q=complex(zeros(n,order));
        AQ=complex(zeros(n,order));
        Q(:,1)=v;
        used=order;
        for t=1:order
            w=U(Q(:,t),'notransp');calls=calls+1;
            AQ(:,t)=w;
            if t==order,break;end
            h=w;
            for pass=1:2
                h=h-Q(:,1:t)*(Q(:,1:t)'*h);
            end
            b=norm(h);
            if b<=1e-12*max(1,norm(w))
                used=t;
                break
            end
            Q(:,t+1)=h/b;
        end
        Q=Q(:,1:used);
        AQ=AQ(:,1:used);
    end

    function S=make_sketch(nrows,type)
        switch lower(char(type))
            case 'gaussian'
                S=(randn(nrows,n)+1i*randn(nrows,n))/sqrt(2*nrows);
            case 'countsketch'
                rows=randi(nrows,n,1);
                phase=exp(0.5i*pi*randi(4,n,1));
                S=sparse(rows,(1:n)',phase,nrows,n);
            otherwise
                error('usd_poisson_density:SketchType', ...
                    'Unknown sketch_type: %s.',char(type));
        end
    end
end
