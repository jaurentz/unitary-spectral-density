function [rho,info]=usd_moment_density(U,theta,K,s,r,seed,opts)
%USD_MOMENT_DENSITY Trace Fourier moments with selectable circle windows.
%   [RHO,INFO]=USD_MOMENT_DENSITY(U,THETA,K,S,R,SEED,OPTS)
%   only calls U(X,'notransp') and U([],'size').
%
%   OPTS.smoothing: 'poisson' (default), 'jackson', 'fejer',
%     'dirichlet', 'lanczos', 'heat', 'poisson-jackson',
%     'poisson-fejer'. HEAT requires OPTS.bandwidth (radians).
%   OPTS.trace_mode: 'stochastic' (default) or 'exact'.
%     EXACT averages all canonical basis vectors, regardless of S;
%     it costs n*K matrix-vector products and is for small matrices.
%   OPTS.probe_type: 'fourphase' (default), 'rademacher',
%     'complex_gaussian', 'basis'; OPTS.probe_vectors: custom vectors.
%   For Jackson/Fejer/Dirichlet/Lanczos/Heat, the positional R is
%   unused unless a POISSON-* composite window is selected.
%
%   Each random probe is unit norm and isotropic, so its moments
%   estimate trace(U^k)/n without bias. Results at different angles
%   reuse probes, giving correlated (not independent) point errors.
if nargin<6 || isempty(seed),seed=17;end
if nargin<7 || isempty(opts),opts=struct;end
sz=U([],'size');n=sz(1);
assert(numel(sz)==2 && sz(2)==n && n>=1,'U must be square.');
validateattributes(K,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
validateattributes(s,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
if ~isfield(opts,'smoothing'),opts.smoothing='poisson';end
if ~isfield(opts,'trace_mode'),opts.trace_mode='stochastic';end
switch lower(char(opts.trace_mode))
    case 'exact'
        assert(~isfield(opts,'probe_vectors') || isempty(opts.probe_vectors), ...
            'Exact trace mode cannot use custom probes.');
        opts.probe_type='basis';
    case 'stochastic'
        % User may also select the complete basis via probe_type.
    otherwise
        error('usd_moment_density:TraceMode', ...
            'trace_mode must be stochastic or exact.');
end
opts.seed=seed;
[V,probeKind]=usd_probe_vectors(n,s,opts);
p=size(V,2);
moments=zeros(K,p);
for j=1:p
    v=V(:,j);
    x=v;
    for k=1:K
        x=U(x,'notransp');
        moments(k,j)=v'*x;
    end
end
g=usd_window_weights(K,r,opts.smoothing,opts);
theta=theta(:);
samples=ones(numel(theta),p)/(2*pi);
for k=1:K
    samples=samples+(g(k)/pi)*real( ...
        exp(-1i*k*theta)*moments(k,:));
end
rho=mean(samples,2);
if p>1 && ~strcmp(probeKind,'basis')
    stderr=std(samples,0,2)/sqrt(p);
else
    stderr=zeros(size(rho));
end
info=struct('moments',mean(moments,2), ...
    'probe_moments',moments,'matvecs',K*p,'probes',p, ...
    'degree',K,'radius',r,'smoothing',lower(char(opts.smoothing)), ...
    'weights',g,'probe_type',probeKind, ...
    'trace_mode',lower(char(opts.trace_mode)), ...
    'standard_error',stderr);
end
