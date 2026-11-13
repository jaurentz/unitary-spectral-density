function tests=test_matrix_library_integration
%TEST_MATRIX_LIBRARY_INTEGRATION Exercise Anymatrix operators via all methods.
%   These tests run twice in CI: pinned and latest anymatrix-unitary.
tests=functiontests(localfunctions);
end

function test_operator_contract_on_distinct_families(tc)
fixtures=matrix_fixtures;
for j=1:numel(fixtures)
    f=fixtures{j};
    A=f.matrix();
    U=f.operator();
    n=size(A,1);
    verifyEqual(tc,U([],'size'),[n n],f.name);
    X=reshape(sin(.17*(1:3*n)),n,3) + ...
      1i*reshape(cos(.23*(1:3*n)),n,3);
    tol=3e-11*max(1,norm(X,'fro'));
    verifyLessThan(tc,norm(U(X,'notransp')-A*X,'fro'),tol,f.name);
    verifyLessThan(tc,norm(U(X,'transp')-A'*X,'fro'),tol,f.name);
    verifyLessThan(tc,norm(U(X,'ctransp')-A'*X,'fro'),tol,f.name);
    verifyLessThan(tc,norm(U(U(X,'notransp'),'transp')-X,'fro'),tol,f.name);
    verifyLessThan(tc,norm(full(A'*A-speye(n)),'fro'),3e-11*n,f.name);
end
end

function test_exact_trace_reconstruction_on_matrix_families(tc)
% The exact reference uses eigenvalues of a SMALL explicit matrix only in
% this test, never inside the generic spectral density estimator.
fixtures=matrix_fixtures;
theta=linspace(-pi,pi,19).';
r=.58; K=60;
opts=struct('trace_mode','exact','smoothing','poisson');
for j=1:numel(fixtures)
    f=fixtures{j};
    A=f.matrix();
    U=f.operator();
    n=size(A,1);
    phases=angle(eig(full(A)));
    exact=usd_poisson_kernel(theta,phases,ones(n,1)/n,r);
    [rho,info]=usd_moment_density(U,theta,K,1,r,19,opts);
    verifyLessThan(tc,max(abs(rho-exact)),1e-10,f.name);
    verifyEqual(tc,info.matvecs,n*K,f.name);
end
end

function test_weighted_szego_trace_on_cmv(tc)
% Using a complete canonical basis eliminates stochastic trace error.
n=10;
alpha=[.2*exp(1i*(1:n-1)'*.37);exp(1i*.5)];
U=finite_cmv(n,alpha,'operator');
A=finite_cmv(n,alpha);
theta=linspace(-pi,pi,19).';
r=.6;
exact=usd_poisson_kernel(theta,angle(eig(full(A))),ones(n,1)/n,r);
opts=struct('probe_type','basis');
[rho,info]=usd_szego_density(U,theta,n,1,r,13,opts);
verifyLessThan(tc,max(abs(rho-exact)),1e-9);
verifyLessThan(tc,info.max_unitarity_residual,1e-9);
verifyEqual(tc,info.probes,n);
end

function test_resolvent_trace_on_prescribed_arcs(tc)
n=10;
phi=[linspace(-2.6,-1.5,n/2),linspace(.4,1.4,n/2)].';
U=phase_circulant(n,phi,'operator');
theta=linspace(-pi,pi,11).';
r=.55;
exact=usd_poisson_kernel(theta,phi,ones(n,1)/n,r);
opts=struct('probe_type','basis');
[rho,info]=usd_poisson_density(U,theta,1,r,13,1e-11,n,n,opts);
verifyLessThan(tc,max(abs(rho-exact)),1e-9);
verifyGreaterThan(tc,info.matvecs,0);
end

function fixtures=matrix_fixtures
% Distinct matrix representations exercise FFT, sparse quantum dynamics,
% and CMV core products. They are independent upstream regression targets.
n=12;
phases=[linspace(-2.6,-1.5,n/2),linspace(.4,1.5,n/2)].';
alpha=[.2*exp(1i*(1:n-1)'*.37);exp(1i*.5)];
fixtures={ ...
    struct('name','Dirac automaton', ...
           'matrix',@()dirac_automaton(6,.45), ...
           'operator',@()dirac_automaton(6,.45,'operator')), ...
    struct('name','Prescribed-phase circulant', ...
           'matrix',@()phase_circulant(n,phases), ...
           'operator',@()phase_circulant(n,phases,'operator')), ...
    struct('name','Finite CMV', ...
           'matrix',@()finite_cmv(n,alpha), ...
           'operator',@()finite_cmv(n,alpha,'operator'))};
end
