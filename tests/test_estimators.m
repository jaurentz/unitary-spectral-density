function tests=test_estimators
%TEST_ESTIMATORS Regression tests for all three portable estimators.
tests=functiontests(localfunctions);
end

function test_moments_exact_diagonal(tc)
[U,theta,exact,r]=diagonal_fixture;
[rho,info]=usd_moment_density(U,theta,80,2,r,23);
verifyLessThan(tc,max(abs(rho-exact)),1e-9);
verifyEqual(tc,info.matvecs,160);
end

function test_szego_exact_diagonal(tc)
[U,theta,exact,r]=diagonal_fixture;
[rho,info]=usd_szego_density(U,theta,12,1,r,23);
verifyLessThan(tc,max(abs(rho-exact)),1e-8);
verifyLessThan(tc,info.max_unitarity_residual,1e-10);
end

function test_resolvent_exact_diagonal(tc)
[U,theta,exact,r]=diagonal_fixture;
[rho,info]=usd_poisson_density(U,theta,2,r,23,1e-12,12,3);
verifyLessThan(tc,max(abs(rho-exact)),1e-8);
verifyGreaterThan(tc,info.matvecs,0);
end

function test_szego_moment_matching(tc)
U=dirac_automaton(16,.5,'operator');
v=usd_fourphase(32);
[nodes,weights,info]=usd_szego_rule(U,v,13);
verifyLessThan(tc,abs(sum(weights)-1),1e-11);
verifyLessThan(tc,info.unitarity_residual,1e-10);
x=v;
for k=1:12
    x=U(x,'notransp');
    actual=v'*x;
    quad=sum(weights.*exp(1i*k*nodes));
    verifyLessThan(tc,abs(actual-quad),3e-9);
end
end

function test_all_three_dirac_drivers(tc)
p=dirac_case(16,.5,.7,33);
a=run_dirac_moments(p,false);
b=run_dirac_quadrature(p,false);
c=run_dirac_filters(p,false);
for d={a,b,c}
    q=d{1};
    verifyTrue(tc,all(isfinite(q.rho)));
    verifyLessThan(tc,q.L1_error,.65);
end
verifyGreaterThanOrEqual(tc,min(b.rho),-1e-12);
verifyGreaterThanOrEqual(tc,min(c.rho),-1e-9);
end

function [U,theta,exact,r]=diagonal_fixture
n=12;phases=2*pi*(0:n-1)'/n;
lambda=exp(1i*phases);
U=@(X,flag)diagonal_apply(X,flag,lambda);
theta=linspace(-pi,pi,41).';
r=.65;
exact=usd_poisson_kernel(theta,phases,ones(n,1)/n,r);
end

function Y=diagonal_apply(X,flag,lambda)
n=numel(lambda);
switch lower(char(flag))
    case 'size',Y=[n n];
    case 'notransp',Y=lambda.*X;
    case {'transp','ctransp'},Y=conj(lambda).*X;
    otherwise,error('test_estimators:Flag','Unknown action.');
end
end
