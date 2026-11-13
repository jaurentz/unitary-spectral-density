function tests=test_experimental_options
%TEST_EXPERIMENTAL_OPTIONS Smoothing, stochastic probes and sketched solves.
tests=functiontests(localfunctions);
end

function test_circle_windows_positive_and_normalized(tc)
theta=linspace(-pi,pi,2001).';
nodes=[-.7;.3;2.4];
weights=[.2;.5;.3];
for method={'fejer','jackson','poisson','poisson-jackson', ...
        'poisson-fejer','heat'}
    opts=struct('bandwidth',.18);
    rho=usd_circular_density(theta,nodes,weights,.65, ...
        method{1},24,opts);
    if ~strcmp(method{1},'heat')
        % Finite heat-series truncation need not preserve positivity.
        verifyGreaterThanOrEqual(tc,min(rho),-2e-12);
    end
    verifyLessThan(tc,abs(trapz(theta,rho)-1),1e-11);
end
g=usd_window_weights(19,.6,'jackson',struct);
verifyEqual(tc,numel(g),19);
verifyLessThan(tc,max(abs(g)),1+1e-12);
verifyEqual(tc,g(end),0,'AbsTol',1e-14); % degree <= 18
end

function test_exact_moments_match_circle_convolution(tc)
n=7;
phi=[-.6;.1;.5;1.7;2.8;-2.3;-.95];
U=@(x,f)diagapply(x,f,exp(1i*phi));
theta=linspace(-pi,pi,85).';
opts=struct('trace_mode','exact','smoothing','jackson');
[rho,info]=usd_moment_density(U,theta,16,2,.73,99,opts);
exact=usd_circular_density(theta,phi,ones(n,1)/n,.73, ...
    'jackson',16,opts);
verifyLessThan(tc,max(abs(rho-exact)),1e-12);
verifyEqual(tc,info.matvecs,n*16);
verifyEqual(tc,info.probes,n);
verifyEqual(tc,info.probe_type,'basis');
verifyLessThan(tc,max(info.standard_error),1e-15);
end

function test_random_moments_reproducible_and_rng_preserved(tc)
n=80;
phases=linspace(-1,2.5,n).';
U=@(x,f)diagapply(x,f,exp(1i*phases));
theta=linspace(-pi,pi,40).';
before=rng;
opts=struct('smoothing','poisson-jackson', ...
    'probe_type','complex_gaussian');
[a,ia]=usd_moment_density(U,theta,12,5,.7,44,opts);
verifyEqual(tc,rng,before);
[b,ib]=usd_moment_density(U,theta,12,5,.7,44,opts);
verifyEqual(tc,a,b,'AbsTol',0);
verifyEqual(tc,ia.probe_moments,ib.probe_moments,'AbsTol',0);
verifyEqual(tc,ia.matvecs,60);
verifyTrue(tc,all(isfinite(ia.standard_error)));
end

function test_large_sparse_stochastic_trace(tc)
n=350;
lambda=exp(1i*linspace(-2.5,2.5,n).');
U=@(x,f)diagapply(x,f,lambda);
opts=struct('probe_type','rademacher','smoothing','fejer');
[rho,info]=usd_moment_density(U,[-1;0;1],9,3,.7,33,opts);
verifyEqual(tc,info.matvecs,27);
verifyEqual(tc,info.probes,3);
verifyTrue(tc,all(isfinite(rho)));
end

function test_random_start_quadrature_repeats(tc)
n=14;
A=diag(exp(1i*(0:n-1)'/n))*fft(eye(n))/sqrt(n);
U=@(x,f)denseapply(x,f,A);
theta=linspace(-pi,pi,41).';
opts=struct('probe_type','complex_gaussian','smoothing','fejer');
[rho,info]=usd_szego_density(U,theta,8,5,.6,57,opts);
[other,again]=usd_szego_density(U,theta,8,5,.6,57,opts);
verifyEqual(tc,rho,other,'AbsTol',0);
verifyEqual(tc,info.standard_error,again.standard_error,'AbsTol',0);
verifyEqual(tc,info.probes,5);
verifyLessThan(tc,info.max_unitarity_residual,1e-9);
verifyTrue(tc,all(isfinite(info.standard_error)));
end

function test_basis_quadrature_recovers_trace(tc)
n=8;phi=[-.7;-.3;.1;.55;1.2;2.3;2.8;-2.1];
U=@(x,f)diagapply(x,f,exp(1i*phi));
theta=linspace(-pi,pi,53).';
opts=struct('probe_type','basis','smoothing','jackson');
[rho,info]=usd_szego_density(U,theta,n,1,.6,5,opts);
exact=usd_circular_density(theta,phi,ones(n,1)/n,.6, ...
    'jackson',n-1,opts);
verifyLessThan(tc,max(abs(rho-exact)),1e-10);
verifyEqual(tc,info.probes,n);
end

function test_gaussian_sketched_exact_full_space(tc)
n=10;
phi=linspace(-2.4,2.1,n).';
U=@(x,f)diagapply(x,f,exp(1i*phi));
theta=linspace(-pi,pi,19).';
opts=struct('solver','sketch','sketch_order',n, ...
    'sketch_dimension',3*n,'sketch_type','gaussian');
[rho,info]=usd_poisson_density(U,theta,2,.6,41, ...
    1e-12,n,n,opts);
exact=usd_poisson_kernel(theta,phi,ones(n,1)/n,.6);
verifyLessThan(tc,max(abs(rho-exact)),1e-9);
verifyLessThan(tc,info.max_relative_residual,1e-9);
verifyEqual(tc,info.matvecs,2*n);
verifyEqual(tc,info.fallback_solves,0);
verifyEqual(tc,info.solver,'sketch');
end

function test_countsketch_fallback_to_gmres(tc)
n=9;
phi=linspace(-2.7,2.2,n).';
U=@(x,f)diagapply(x,f,exp(1i*phi));
theta=[-.9;0;.8];
opts=struct('solver','sketch','sketch_order',2, ...
    'sketch_dimension',12,'sketch_type','countsketch', ...
    'sketch_fallback',true,'sketch_residual_tolerance',1e-10);
[rho,info]=usd_poisson_density(U,theta,1,.6,13, ...
    1e-11,n,n,opts);
exact=usd_poisson_kernel(theta,phi,ones(n,1)/n,.6);
verifyGreaterThan(tc,info.fallback_solves,0);
verifyLessThan(tc,info.max_relative_residual,1e-8);
verifyLessThan(tc,max(abs(rho-exact)),1e-8);
end

function out=diagapply(x,flag,lambda)
n=numel(lambda);
switch lower(char(flag))
    case 'size',out=[n n];
    case 'notransp',out=lambda.*x;
    case {'transp','ctransp'},out=conj(lambda).*x;
    otherwise,error('test_experimental_options:Flag','Unknown flag');
end
end

function out=denseapply(x,flag,A)
switch lower(char(flag))
    case 'size',out=size(A);
    case 'notransp',out=A*x;
    case {'transp','ctransp'},out=A'*x;
    otherwise,error('test_experimental_options:Flag','Unknown flag');
end
end
