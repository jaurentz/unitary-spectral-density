function g=usd_window_weights(K,r,method,opts)
%USD_WINDOW_WEIGHTS Fourier multipliers for circular spectral smoothing.
%   g(k) multiplies normalized trace moment tr(U^k)/n, k=1,...,K.
%   POISSON is r^k; FEJER, JACKSON and DIRICHLET do not use r.
%   POISSON-JACKSON and POISSON-FEJER explicitly compose windows.
%   JACKSON is the normalized SQUARE of the circle Fejer kernel,
%   NOT the Hermitian/Chebyshev KPM Jackson coefficients.
if nargin<4 || isempty(opts),opts=struct;end
if nargin<3 || isempty(method),method='poisson';end
method=lower(char(method));
validateattributes(K,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
if strncmp(method,'poisson-',8)
    base=usd_window_weights(K,r,method(9:end),opts);
    assert(isscalar(r) && isreal(r) && r>0 && r<1, ...
        'Poisson radius must satisfy 0<r<1.');
    g=(r.^(1:K)').*base;
    return
end
k=(1:K)';
switch method
    case 'poisson'
        assert(isscalar(r) && isreal(r) && r>0 && r<1, ...
            'Poisson radius must satisfy 0<r<1.');
        g=r.^k;
    case 'fejer'
        g=1-k/(K+1);
    case 'jackson'
        q=floor(K/2)+1;
        a=q-abs(-(q-1):(q-1));
        c=conv(a,a);
        mid=2*q-1;
        g=zeros(K,1);
        last=min(K,2*(q-1));
        g(1:last)=c(mid+(1:last)).'/c(mid);
    case {'none','dirichlet'}
        g=ones(K,1);
    case 'lanczos'
        t=k/(K+1);
        g=sin(pi*t)./(pi*t);
    case {'heat','wrapped_gaussian'}
        assert(isfield(opts,'bandwidth') && ...
            isscalar(opts.bandwidth) && opts.bandwidth>0, ...
            'HEAT requires a positive opts.bandwidth (radians).');
        g=exp(-0.5*(opts.bandwidth*k).^2);
    otherwise
        error('usd_window_weights:Method', ...
            'Unknown smoothing method: %s.',method);
end
end
