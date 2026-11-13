%% 1. Dirac quantum cellular automaton: two spectral arcs
% DIRAC SPECTRAL DENSITY -- a five-section, plain-text Live Editor example.
%
% The one-dimensional periodic Dirac automaton acts on 2*n states:
%
%   U = [nu*P, -1i*m*I; -1i*m*I, nu*P'],  nu = sqrt(1-m^2),
%
% where P is a cyclic shift. Its EXACT finite eigenphases are
%
%   +/-acos(nu*cos(2*pi*k/n)),  k=0,...,n-1.
%
% For 0<m<1, the eigenphases lie in the two arcs
%   [asin(m),pi-asin(m)] and [-pi+asin(m),-asin(m)].
%
% The finite eigenphase-counting measure is atomic, not a smooth density.
% Accordingly, all three estimators below target the SAME Poisson-smoothed
% empirical measure, with smoothing radius r. The analytic limiting
% (unsmoothed) density is plotted here as additional context.
%
% This is an ordinary, version-control-friendly .m script. MATLAB's Live
% Editor recognizes each %% section and can run the sections in sequence.
% Set showPlots=false before running to suppress figures (e.g. in CI).
%
% Experiment settings may be overridden by assigning exampleSites,
% exampleGridPoints, exampleProbes, exampleRadius, exampleMomentDegree,
% and exampleArnoldiOrder before running this script.

exampleRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(exampleRoot);
setup;

if ~exist('showPlots','var') || isempty(showPlots), showPlots = true; end
if ~exist('exampleSites','var') || isempty(exampleSites), exampleSites = 48; end
if ~exist('exampleGridPoints','var') || isempty(exampleGridPoints)
    exampleGridPoints = 97;
end
if ~exist('exampleProbes','var') || isempty(exampleProbes), exampleProbes = 6; end
if ~exist('exampleRadius','var') || isempty(exampleRadius), exampleRadius = .80; end
if ~exist('exampleMomentDegree','var') || isempty(exampleMomentDegree)
    exampleMomentDegree = 70;
end
if ~exist('exampleArnoldiOrder','var') || isempty(exampleArnoldiOrder)
    exampleArnoldiOrder = min(28,2*exampleSites);
end

mass = .45;
seed = 17;
problem = dirac_case(exampleSites,mass,exampleRadius,exampleGridPoints);
U = problem.U;                 % Matrix-free only: U(X,'notransp'), U([],'size')
theta = problem.theta;        % Eigenphase grid in [-pi,pi]
rhoTrue = problem.exact_smoothed; % EXACT finite, Poisson-smoothed trace density
rhoLimit = problem.limiting_density; % Analytic LIMITING unsmoothed density
r = problem.radius;
bandEdge = problem.band_edges(1);

if showPlots
    figure('Name','Dirac automaton: spectral arcs');
    hExact = plot(theta,rhoTrue,'LineWidth',1.8); hold on
    hLimit = plot(theta,rhoLimit,'--','LineWidth',1.2);
    hEigen = plot(problem.exact_phases,zeros(size(problem.exact_phases)), ...
        '.','MarkerSize',7);
    xline([-pi+bandEdge,-bandEdge,bandEdge,pi-bandEdge],':');
    hold off
    xlim([-pi pi]); xlabel('Eigenphase \theta (radians)');
    ylabel('Spectral density');
    title(sprintf('Dirac automaton: 2n=%d, mass=%.2f, r=%.2f', ...
        problem.dimension,mass,r));
    legend([hExact hLimit hEigen], ...
        {'Exact smoothed finite measure','Limiting density', ...
         'Exact finite eigenphases'},'Location','best');
    grid on
end

%% 2. Trace/moment methods: stochastic Fourier reconstruction
% With normalized isotropic probes v, we estimate
%
%   m_k = tr(U^k)/(2*n)  from  v'*U^k*v.
%
% The estimator reconstructs a truncated Poisson series:
%
%   rho(theta) = [1 + 2*Re(sum_{k=1}^K r^k*m_k*exp(-1i*k*theta))]/(2*pi).
%
% No explicit Dirac matrix or eigenvalues are passed to the estimator.

[rhoMoments,momentsInfo] = usd_moment_density( ...
    U,theta,exampleMomentDegree,exampleProbes,r,seed);
errMoments = trapz(theta,abs(rhoMoments-rhoTrue));
fprintf('Moments:    L1 error = %.4g; matvecs = %d\n', ...
    errMoments,momentsInfo.matvecs);

if showPlots
    figure('Name','Dirac density: trace moments');
    plot(theta,rhoTrue,'k-','LineWidth',1.8); hold on
    plot(theta,rhoMoments,'--','LineWidth',1.5); hold off
    xlim([-pi pi]); xlabel('Eigenphase \theta (radians)');
    ylabel('Poisson-smoothed trace density');
    title('Stochastic Fourier moments');
    legend('Exact smoothed finite density','Moment estimate','Location','best');
    grid on
end

%% 3. Isometric Arnoldi and Gauss-Szego quadrature
% Full reorthogonalization builds a unitary Hessenberg quadrature matrix.
% Its NODES and (essential) spectral WEIGHTS approximate each starting
% vector's spectral measure. Averaging over random starting vectors gives
% a trace-density estimator. The reduced matrix eigenvalues must NOT be
% given equal weights.
%
% This initial implementation uses Arnoldi plus a unitary completion; it
% is not yet a specialized short-recurrence isometric Arnoldi algorithm.

[rhoArnoldi,arnoldiInfo] = usd_szego_density( ...
    U,theta,exampleArnoldiOrder,exampleProbes,r,seed);
errArnoldi = trapz(theta,abs(rhoArnoldi-rhoTrue));
fprintf('Arnoldi:    L1 error = %.4g; matvecs = %d\n', ...
    errArnoldi,arnoldiInfo.matvecs);

if showPlots
    figure('Name','Dirac density: Gauss-Szego quadrature');
    plot(theta,rhoTrue,'k-','LineWidth',1.8); hold on
    plot(theta,rhoArnoldi,'--','LineWidth',1.5); hold off
    xlim([-pi pi]); xlabel('Eigenphase \theta (radians)');
    ylabel('Poisson-smoothed trace density');
    title('Isometric Arnoldi / Gauss-Szego quadrature');
    legend('Exact smoothed finite density','Arnoldi estimate', ...
        'Location','best');
    grid on
end

%% 4. Kernel methods: direct Poisson-resolvent filtering
% A third numerical strategy evaluates the same smoothing kernel through
% shifted linear solves, rather than through moments or quadrature:
%
%   z = r*exp(1i*theta)
%   rho_r(theta) = Re tr[(U+z*I)*(U-z*I)^{-1}]/(2*pi*2*n).
%
% Here stochastic trace estimation and MATLAB's GMRES use ONLY the
% matrix-free action U(x,'notransp'). This route can be considerably more
% expensive when many evaluation angles or sharp kernels are requested.

[rhoKernel,kernelInfo] = usd_poisson_density( ...
    U,theta,exampleProbes,r,seed);
errKernel = trapz(theta,abs(rhoKernel-rhoTrue));
fprintf('Kernel:     L1 error = %.4g; matvecs = %d\n', ...
    errKernel,kernelInfo.matvecs);

if showPlots
    figure('Name','Dirac density: Poisson-resolvent kernel');
    plot(theta,rhoTrue,'k-','LineWidth',1.8); hold on
    plot(theta,rhoKernel,'--','LineWidth',1.5); hold off
    xlim([-pi pi]); xlabel('Eigenphase \theta (radians)');
    ylabel('Poisson-smoothed trace density');
    title('Direct kernel / shifted-resolvent estimate');
    legend('Exact smoothed finite density','Kernel estimate', ...
        'Location','best');
    grid on
end

%% 5. Compare all three techniques with the true smoothed density
% All four curves use the SAME smoothing radius and phase grid.
% The exact reference comes from known finite Dirac eigenphases.
% Different matvec counts mean this is a demonstration, not a comparison
% at equal computational budgets. In particular, the direct kernel method
% solves a separate linear system at each point of the phase grid.

method = {'Fourier moments';'Arnoldi / Gauss-Szego';'Poisson resolvent'};
L1Error = [errMoments;errArnoldi;errKernel];
Matvecs = [momentsInfo.matvecs;arnoldiInfo.matvecs;kernelInfo.matvecs];
comparison = table(method,L1Error,Matvecs);
disp(comparison)

% Named arrays remain in the workspace for further Live Editor exploration:
% rhoTrue, rhoMoments, rhoArnoldi, rhoKernel, theta, comparison.
if showPlots
    figure('Name','Dirac density: three estimator comparison');
    plot(theta,rhoTrue,'k-','LineWidth',2); hold on
    plot(theta,rhoMoments,'--','LineWidth',1.3);
    plot(theta,rhoArnoldi,'-.','LineWidth',1.3);
    plot(theta,rhoKernel,':','LineWidth',1.8);
    hold off
    xlim([-pi pi]); xlabel('Eigenphase \theta (radians)');
    ylabel('Poisson-smoothed trace density');
    title(sprintf('Three methods vs exact Dirac density (r=%.2f)',r));
    legend('Exact smoothed finite density','Moments', ...
        'Arnoldi / Gauss-Szego','Direct kernel','Location','best');
    grid on
end
