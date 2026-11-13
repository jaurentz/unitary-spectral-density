function p=dirac_case(nsites,mass,r,ngrid)
%DIRAC_CASE Two-arc Dirac automaton and exact finite Poisson reference.
%   The model-specific construction is isolated from estimator routines.
if nargin<1,nsites=64;end
if nargin<2,mass=.45;end
if nargin<3,r=.82;end
if nargin<4,ngrid=97;end
theta=linspace(-pi,pi,ngrid).';
nu=sqrt((1-mass)*(1+mass));
k=(0:nsites-1).';
omega=acos(nu*cos(2*pi*k/nsites));
phases=[omega;-omega];
reference=usd_poisson_kernel(theta,phases, ...
                            ones(2*nsites,1)/(2*nsites),r);
limiting=zeros(size(theta));
inside=abs(cos(theta))<nu;
limiting(inside)=abs(sin(theta(inside)))./ ...
    (2*pi*sqrt(nu^2-cos(theta(inside)).^2));
p=struct('U',dirac_automaton(nsites,mass,'operator'), ...
    'sites',nsites,'dimension',2*nsites,'mass',mass, ...
    'theta',theta,'radius',r,'exact_phases',phases, ...
    'exact_smoothed',reference,'limiting_density',limiting, ...
    'band_edges',[asin(mass),pi-asin(mass)]);
end
