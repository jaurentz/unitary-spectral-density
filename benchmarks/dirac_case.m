function p=dirac_case(nsites,mass,r,ngrid)
%DIRAC_CASE Dirac automaton using the Anymatrix spectral reference.
if nargin<1,nsites=64;end
if nargin<2,mass=.45;end
if nargin<3,r=.82;end
if nargin<4,ngrid=97;end
theta=linspace(-pi,pi,ngrid).';
S=dirac_automaton(nsites,mass,'spectral');
phases=angle(S.finite.eigenvalues());
reference=usd_poisson_kernel(theta,phases, ...
    ones(S.dimension,1)/S.dimension,r);
limiting=zeros(size(theta));
if isa(S.limit.density,'function_handle')
    limiting=S.limit.density(theta);
end
p=struct('U',dirac_automaton(nsites,mass,'operator'), ...
    'spectral',S,'sites',nsites,'dimension',S.dimension, ...
    'mass',mass,'theta',theta,'radius',r, ...
    'exact_phases',phases,'exact_smoothed',reference, ...
    'limiting_density',limiting, ...
    'band_edges',[asin(mass),pi-asin(mass)]);
end
