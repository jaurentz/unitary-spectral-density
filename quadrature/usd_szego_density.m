function [rho,info]=usd_szego_density(U,theta,m,s,r,seed,opts)
%USD_SZEGO_DENSITY Averaged random-start Gauss-Szego quadrature.
%   [RHO,INFO]=USD_SZEGO_DENSITY(U,THETA,M,S,R,SEED,OPTS).
%   The m-node rule matches each starting vector's unitary moments
%   through degree m-1. Its NONUNIFORM spectral weights are essential.
%   Averaging independent isotropic starts estimates the trace measure.
%
%   OPTS.probe_type: 'fourphase' (default), 'rademacher',
%      'complex_gaussian', 'basis'; OPTS.probe_vectors: custom starts.
%   OPTS.smoothing: 'poisson' (default), or a method supported by
%      USD_WINDOW_WEIGHTS. For finite windows, OPTS.degree defaults
%      to m-1, ensuring exact moment transfer from the quadrature.
%   Returned STANDARD_ERROR is the pointwise Monte Carlo SE, not
%   the deterministic Arnoldi/quadrature truncation error.
if nargin<6 || isempty(seed),seed=17;end
if nargin<7 || isempty(opts),opts=struct;end
sz=U([],'size');n=sz(1);
assert(numel(sz)==2 && sz(2)==n,'U must be square.');
validateattributes(s,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
validateattributes(m,{'numeric'}, ...
    {'scalar','integer','positive','<=',n});
if ~isfield(opts,'smoothing'),opts.smoothing='poisson';end
if ~isfield(opts,'degree') || isempty(opts.degree)
    opts.degree=m-1;
end
validateattributes(opts.degree,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
if opts.degree>m-1 && ~strcmpi(char(opts.smoothing),'poisson')
    warning('usd_szego_density:Degree', ...
        'Requested smoothing degree exceeds guaranteed quadrature exactness.');
end
opts.seed=seed;
[V,kind]=usd_probe_vectors(n,s,opts);
p=size(V,2);
theta=theta(:);
samples=zeros(numel(theta),p);
counts=zeros(p,1);unitarity=zeros(p,1);sizes=zeros(p,1);
for j=1:p
    [nodes,weights,rule]=usd_szego_rule(U,V(:,j),m);
    samples(:,j)=usd_circular_density(theta,nodes,weights, ...
        r,opts.smoothing,opts.degree,opts);
    counts(j)=rule.matvecs;
    unitarity(j)=rule.unitarity_residual;
    sizes(j)=rule.nodes;
end
rho=mean(samples,2);
if p>1 && ~strcmp(kind,'basis')
    stderr=std(samples,0,2)/sqrt(p);
else
    stderr=zeros(size(rho));
end
info=struct('matvecs',sum(counts),'probes',p, ...
    'quadrature_order',m,'effective_orders',sizes, ...
    'max_unitarity_residual',max(unitarity),'radius',r, ...
    'probe_type',kind,'smoothing',lower(char(opts.smoothing)), ...
    'degree',opts.degree,'standard_error',stderr);
end
