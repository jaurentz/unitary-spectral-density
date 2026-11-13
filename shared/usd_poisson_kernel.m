function rho=usd_poisson_kernel(theta,nodes,weights,r)
%USD_POISSON_KERNEL Weighted circular Poisson density (unit integral).
%   THETA and NODES are angles in radians. WEIGHTS sum to one.
theta=theta(:);nodes=nodes(:);weights=weights(:);
if numel(nodes)~=numel(weights)
    error('usd_poisson_kernel:Weights','Node/weight lengths differ.');
end
if ~(isscalar(r) && isreal(r) && r>0 && r<1)
    error('usd_poisson_kernel:Radius','Radius must satisfy 0<r<1.');
end
d=theta-nodes.';
p=(1-r^2)./(1+r^2-2*r*cos(d));
rho=real(p*weights)/(2*pi);
end
