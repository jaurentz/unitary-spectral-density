%[text] # Spectral density of a Dirac quantum automaton
%[text] A periodic Dirac quantum automaton is a unitary matrix with eigenvalues on two arcs of the circle. We compare three methods of estimating its eigenphase density with a known finite-dimensional reference. The full mathematics and additional options are in [mathematics.md](../mathematics.md).
%%
%[text] ## 1. The Dirac operator and its spectrum
%[text] Anymatrix-unitary supplies both a matrix-vector operator and analytic spectral information. We use its exact finite eigenvalues to plot the Poisson-smoothed density.
addpath(fileparts(fileparts(which('dirac_spectral_density'))))
setup
sites = 48; mass = 0.45; r = 0.8;
theta = linspace(-pi,pi,121).';
U = dirac_automaton(sites,mass,'operator');
S = dirac_automaton(sites,mass,'spectral');
lambda = S.finite.eigenvalues();
rhoTrue = usd_poisson_kernel(theta,angle(lambda), ...
    ones(numel(lambda),1)/numel(lambda),r);
plot(theta,rhoTrue,'k','LineWidth',1.7)
xlabel('Eigenphase'), ylabel('Density'), grid on
title('Exact finite spectrum (Poisson smoothed)')
%%
%[text] ## 2. Fourier moments
%[text] Randomized traces of powers of the operator give Fourier moments. A truncated Poisson series turns these moments into a density. Jackson and Fejér are other available smoothing options.
K = 40; probes = 6; seed = 17;
[rhoMoments,momentsInfo] = usd_moment_density(U,theta,K,probes,r,seed);
plot(theta,rhoTrue,'k',theta,rhoMoments,'r--','LineWidth',1.5)
xlabel('Eigenphase'), ylabel('Density'), grid on
legend('Exact','Fourier moments','Location','best')
title('Moment reconstruction')
%%
%[text] ## 3. Arnoldi and Gauss-Szegő quadrature
%[text] A small unitary Hessenberg matrix supplies quadrature nodes and weights for each randomized starting vector. We average their Poisson densities. The modest order deliberately leaves a visible quadrature error.
m = 12;
[rhoArnoldi,arnoldiInfo] = usd_szego_density(U,theta,m,probes,r,seed);
plot(theta,rhoTrue,'k',theta,rhoArnoldi,'b--','LineWidth',1.5)
xlabel('Eigenphase'), ylabel('Density'), grid on
legend('Exact','Gauss-Szego','Location','best')
title('Random-start unitary Arnoldi')
%%
%[text] ## 4. A Poisson resolvent
%[text] The same density can be estimated using shifted linear systems. For small matrices, the default solver uses MATLAB backslash; larger problems use GMRES, with optional randomized projected solves.
[rhoKernel,kernelInfo] = usd_poisson_density(U,theta,probes,r,seed);
plot(theta,rhoTrue,'k',theta,rhoKernel,'m--','LineWidth',1.5)
xlabel('Eigenphase'), ylabel('Density'), grid on
legend('Exact','Poisson resolvent','Location','best')
title(['Resolvent with ' kernelInfo.solver])
%%
%[text] ## 5. Comparison
%[text] The three estimators target the same Poisson-smoothed finite spectral measure. Their matrix-vector counts differ, so this is a demonstration rather than a comparison at equal work.
plot(theta,rhoTrue,'k',theta,rhoMoments,'r--', ...
     theta,rhoArnoldi,'b-.',theta,rhoKernel,'m:','LineWidth',1.5)
xlabel('Eigenphase'), ylabel('Density'), grid on
legend('Exact','Moments','Arnoldi','Resolvent','Location','best')
title('Three methods for the Dirac automaton')
Method = {'Moments';'Arnoldi';'Resolvent'};
L1Error = [trapz(theta,abs(rhoMoments-rhoTrue)); ...
           trapz(theta,abs(rhoArnoldi-rhoTrue)); ...
           trapz(theta,abs(rhoKernel-rhoTrue))];
Matvecs = [momentsInfo.matvecs;arnoldiInfo.matvecs;kernelInfo.matvecs];
comparison = table(Method,L1Error,Matvecs)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
