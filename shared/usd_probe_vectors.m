function [V,kind]=usd_probe_vectors(n,s,opts)
%USD_PROBE_VECTORS Unit-norm vectors for trace/quadrature experiments.
%   Options: probe_type = 'fourphase' (default), 'rademacher',
%   'complex_gaussian', or 'basis'. probe_vectors overrides these.
%   Random probes have E[v*v']=I/n. The complete canonical basis
%   gives an exact trace average. Custom probes need not be isotropic.
if nargin<3 || isempty(opts),opts=struct;end
if isfield(opts,'probe_vectors') && ~isempty(opts.probe_vectors)
    V=opts.probe_vectors;
    assert(isnumeric(V) && size(V,1)==n && size(V,2)>0, ...
        'Custom probes must have n rows and at least one column.');
    lengths=sqrt(sum(abs(V).^2,1));
    assert(all(isfinite(lengths) & lengths>0), ...
        'Custom probes must be finite and nonzero.');
    V=bsxfun(@rdivide,V,lengths);
    kind='custom';
    return
end
kind='fourphase';
if isfield(opts,'probe_type') && ~isempty(opts.probe_type)
    kind=lower(char(opts.probe_type));
end
if strcmp(kind,'basis')
    V=speye(n);
    return
end
validateattributes(s,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
seed=17;
if isfield(opts,'seed') && ~isempty(opts.seed),seed=opts.seed;end
state=rng;cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
rng(seed,'twister');
V=complex(zeros(n,s));
for j=1:s
    switch kind
        case 'fourphase'
            V(:,j)=usd_fourphase(n);
        case 'rademacher'
            V(:,j)=(2*randi(2,n,1)-3)/sqrt(n);
        case {'complex_gaussian','gaussian'}
            z=randn(n,1)+1i*randn(n,1);
            V(:,j)=z/norm(z);
        otherwise
            error('usd_probe_vectors:Type', ...
                'Unknown probe_type: %s.',kind);
    end
end
end
