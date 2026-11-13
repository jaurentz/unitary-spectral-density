function [rho,info]=usd_szego_density(U,theta,m,s,r,seed)
%USD_SZEGO_DENSITY Stochastic isometric-Arnoldi/Gauss-Szego density.
%   Random vector spectral measures are averaged to recover the trace
%   spectral measure in expectation. Keep quadrature weights; do not
%   use the equal-weight spectrum of the reduced unitary matrices.
if nargin<6,seed=17;end
sz=U([],'size');n=sz(1);
assert(numel(sz)==2 && sz(2)==n,'U must be square.');
validateattributes(s,{'numeric'},{'scalar','integer','positive'});
state=rng;restore=onCleanup(@()rng(state)); %#ok<NASGU>
rng(seed,'twister');
theta=theta(:);rho=zeros(size(theta));
counts=zeros(s,1);unitarity=zeros(s,1);sizes=zeros(s,1);
for j=1:s
    v=usd_fourphase(n);
    [nodes,weights,rule]=usd_szego_rule(U,v,m);
    rho=rho+usd_poisson_kernel(theta,nodes,weights,r)/s;
    counts(j)=rule.matvecs;
    unitarity(j)=rule.unitarity_residual;
    sizes(j)=rule.nodes;
end
info=struct('matvecs',sum(counts),'probes',s, ...
    'quadrature_order',m,'effective_orders',sizes, ...
    'max_unitarity_residual',max(unitarity),'radius',r);
end
