function rho=usd_circular_density(theta,nodes,weights,r,method,K,opts)
%USD_CIRCULAR_DENSITY Smoothed weighted discrete measure on the circle.
%   The multipliers for non-Poisson windows truncate at degree K.
if nargin<7 || isempty(opts),opts=struct;end
if nargin<6 || isempty(K),K=0;end
if nargin<5 || isempty(method),method='poisson';end
if strcmpi(char(method),'poisson')
    rho=usd_poisson_kernel(theta,nodes,weights,r);
    return
end
theta=theta(:);nodes=nodes(:);weights=weights(:);
assert(numel(nodes)==numel(weights),'Node/weight mismatch.');
g=usd_window_weights(K,r,method,opts);
rho=ones(size(theta))*real(sum(weights))/(2*pi);
for k=1:K
    mk=sum(weights.*exp(1i*k*nodes));
    rho=rho+(g(k)/pi)*real(mk*exp(-1i*k*theta));
end
end
