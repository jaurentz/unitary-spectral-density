%%
%[text] ## 1. Dirac quantum cellular automaton: two spectral arcs
% DIRAC SPECTRAL DENSITY -- a plain-text, five-section MATLAB Live script.
%
% The periodic one-dimensional Dirac quantum cellular automaton is
%
%   U = [nu*P, -1i*m*I; -1i*m*I, nu*P'],    nu = sqrt(1-m^2).
%
% Here P is the cyclic shift on n sites, U has order N=2*n, and the
% EXACT finite eigenphases are
%
%   +/-acos(nu*cos(2*pi*k/n)),   k=0,...,n-1.
%
% For 0<m<1 they occupy the two arcs
%   +/-[asin(m), pi-asin(m)].
%
% The finite empirical trace measure is ATOMIC. All comparisons below use
% analytically known finite eigenphases convolved with the EXACT SAME
% circle kernel as the numerical estimator. In particular, Jackson,
% Fejer, and heat smoothed references are NOT Poisson references.
% The limiting unsmoothed density is shown only as additional context.
%
% Each numerical estimator receives only U(X,'notransp'), U([],'size').
% R2025a+ recognizes the plain-text Live Script appendix and headings.
% The figures consist of five windows, with multiple panels in 2--5.
%
% Set any of these in the workspace to override defaults BEFORE running:
% showPlots, exampleSites, exampleGridPoints, exampleProbes, exampleRadius,
% exampleMomentDegree, exampleArnoldiOrder, exampleStudyDegree,
% exampleStudyProbes, exampleSketchOrder, exampleSketchDimension,
% exampleKernelStudyGridPoints, exampleHeatBandwidth.
%
% The full sweep demonstrates all smoothing choices, exact/randomized
% traces, custom and basis probes, randomized Szego starting vectors,
% direct backslash, GMRES, Gaussian/CountSketch and GMRES fallback.
% For the mathematics and corresponding Hermitian references see
% ../mathematics.md.
exampleRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(exampleRoot);
setup;
if ~exist('showPlots','var') || isempty(showPlots),showPlots=true;end
if ~exist('exampleSites','var') || isempty(exampleSites),exampleSites=48;end
if ~exist('exampleGridPoints','var') || isempty(exampleGridPoints)
    exampleGridPoints=97;
end
if ~exist('exampleProbes','var') || isempty(exampleProbes),exampleProbes=6;end
if ~exist('exampleRadius','var') || isempty(exampleRadius),exampleRadius=.80;end
if ~exist('exampleMomentDegree','var') || isempty(exampleMomentDegree)
    exampleMomentDegree=70;
end
if ~exist('exampleArnoldiOrder','var') || isempty(exampleArnoldiOrder)
    exampleArnoldiOrder=min(28,2*exampleSites);
end
if ~exist('exampleStudyDegree','var') || isempty(exampleStudyDegree)
    % Matched finite-window studies use K <= m-1.
    exampleStudyDegree=min([exampleMomentDegree, ...
        exampleArnoldiOrder-1,24]);
end
if ~exist('exampleStudyProbes','var') || isempty(exampleStudyProbes)
    exampleStudyProbes=min(exampleProbes,4);
end
if ~exist('exampleSketchOrder','var') || isempty(exampleSketchOrder)
    exampleSketchOrder=min(exampleArnoldiOrder,24);
end
if ~exist('exampleSketchDimension','var') || isempty(exampleSketchDimension)
    exampleSketchDimension=2*exampleSketchOrder;
end
if ~exist('exampleKernelStudyGridPoints','var') || ...
        isempty(exampleKernelStudyGridPoints)
    exampleKernelStudyGridPoints=25;
end
if ~exist('exampleHeatBandwidth','var') || isempty(exampleHeatBandwidth)
    exampleHeatBandwidth=.18; % radians (wrapped heat/Gaussian kernel)
end
mass=.45;
seed=17;
problem=dirac_case(exampleSites,mass,exampleRadius,exampleGridPoints);
U=problem.U;
theta=problem.theta;
r=problem.radius;
N=problem.dimension;
rhoTrue=problem.exact_smoothed;       % exact finite Poisson convolution
rhoLimit=problem.limiting_density;   % unsmoothed infinite-size limit
bandEdge=problem.band_edges(1);
finiteWeights=ones(N,1)/N;
% A separate, smaller study budget keeps the comprehensive option sweep
% suitable for interactive use without changing the primary demonstration.
studyDegree=max(0,min([exampleStudyDegree,exampleMomentDegree, ...
    exampleArnoldiOrder-1]));
studyProbes=max(1,exampleStudyProbes);
studyOrder=min(exampleArnoldiOrder,N);
sketchOrder=max(1,min(exampleSketchOrder,N));
sketchDimension=max(sketchOrder,exampleSketchDimension);
kernelGridSize=max(3,min(exampleKernelStudyGridPoints,numel(theta)));
kernelTheta=linspace(-pi,pi,kernelGridSize).';
kernelReference=usd_poisson_kernel(kernelTheta, ...
    problem.exact_phases,finiteWeights,r);
fprintf('\nDirac: N=%d, mass=%.2f, Poisson r=%.2f\n',N,mass,r);
fprintf('Main: K=%d, Arnoldi m=%d, s=%d. Sweeps: K=%d, s=%d.\n', ...
    exampleMomentDegree,studyOrder,exampleProbes,studyDegree,studyProbes);
if showPlots
    figure('Name','1. Dirac automaton: spectral arcs');
    hExact=plot(theta,rhoTrue,'LineWidth',1.8);hold on
    hLimit=plot(theta,rhoLimit,'--','LineWidth',1.2);
    hEigen=plot(problem.exact_phases,zeros(size(problem.exact_phases)), ...
        '.','MarkerSize',7);
    xline([-pi+bandEdge,-bandEdge,bandEdge,pi-bandEdge],':');
    hold off
    xlim([-pi pi]);xlabel('Eigenphase \theta (radians)');
    ylabel('Spectral density');
    title(sprintf('Dirac QCA: N=%d, mass=%.2f, r=%.2f',N,mass,r));
    legend([hExact hLimit hEigen], ...
        {'Exact finite Poisson smoothing','Unsmoothed limiting law', ...
        'Exact finite eigenphases'},'Location','best');
    grid on
end
%%
%[text] ## 2. Trace/moment methods: every circle window and probe option
% Each isotropic unit vector v gives E[v'*U^k*v]=tr(U^k)/N.
% The reconstruction is
%
%   rho(theta) = [1+2*Re(sum_{k=1}^K g_k*m_k*exp(-1i*k*theta))]/(2*pi).
%
% Demonstrate all implemented smoothing labels. Jackson here means the
% positive squared-Fejer CIRCLE kernel; it does NOT use the real-line
% Chebyshev KPM Jackson coefficients. The heat Fourier SERIES is truncated
% and therefore may take slightly negative values despite its positive
% infinite kernel. Dirichlet/none are equivalent aliases; heat and
% wrapped_gaussian are equivalent aliases.
%
% Separate trace/probe experiments show stochastic fourphase, Rademacher,
% normalized complex Gaussian, all canonical basis vectors, custom
% prescribed vectors and trace_mode='exact'. For custom vectors an
% unbiased trace estimator is NOT guaranteed.
[rhoMoments,momentsInfo]=usd_moment_density( ...
    U,theta,exampleMomentDegree,exampleProbes,r,seed);
errMoments=trapz(theta,abs(rhoMoments-rhoTrue));
fprintf('Moments (Poisson): L1 %.4g, %d matvecs\n', ...
    errMoments,momentsInfo.matvecs);
windowNames={'poisson','jackson','fejer','dirichlet','none', ...
    'lanczos','heat','wrapped_gaussian', ...
    'poisson-fejer','poisson-jackson'};
numWindows=numel(windowNames);
windowCurves=zeros(numel(theta),numWindows);
windowTruth=zeros(numel(theta),numWindows);
windowError=zeros(numWindows,1);
windowMatvecs=zeros(numWindows,1);
windowMeanSE=zeros(numWindows,1);
windowInfo=cell(numWindows,1);
for j=1:numWindows
    opts=struct('smoothing',windowNames{j}, ...
        'probe_type','fourphase','bandwidth',exampleHeatBandwidth);
    [windowCurves(:,j),windowInfo{j}]=usd_moment_density( ...
        U,theta,studyDegree,studyProbes,r,seed,opts);
    % Full Poisson reference for the infinite kernel, finite Fourier
    % reference for the other smoothing windows.
    windowTruth(:,j)=usd_circular_density( ...
        theta,problem.exact_phases,finiteWeights,r, ...
        windowNames{j},studyDegree,opts);
    windowError(j)=trapz(theta, ...
        abs(windowCurves(:,j)-windowTruth(:,j)));
    windowMatvecs(j)=windowInfo{j}.matvecs;
    windowMeanSE(j)=mean(windowInfo{j}.standard_error);
end
windowSummary=table(string(windowNames(:)),windowError,windowMatvecs, ...
    windowMeanSE,'VariableNames', ...
    {'Window','L1Error','Matvecs','MeanSE'});
probeChoices={'fourphase','rademacher','complex_gaussian', ...
    'basis','custom','exact'};
numProbeChoices=numel(probeChoices);
probeCurves=zeros(numel(theta),numProbeChoices);
probeError=zeros(numProbeChoices,1);
probeMatvecs=zeros(numProbeChoices,1);
probeMeanSE=zeros(numProbeChoices,1);
probeInfo=cell(numProbeChoices,1);
jacksonTruth=usd_circular_density(theta,problem.exact_phases, ...
    finiteWeights,r,'jackson',studyDegree,struct);
for j=1:numProbeChoices
    opts=struct('smoothing','jackson');
    if strcmp(probeChoices{j},'exact')
        opts.trace_mode='exact'; % all N basis vectors, N*K matvecs
    elseif strcmp(probeChoices{j},'custom')
        opts.probe_vectors=full(speye(N,min(N,2)));
    else
        opts.probe_type=probeChoices{j};
    end
    [probeCurves(:,j),probeInfo{j}]=usd_moment_density( ...
        U,theta,studyDegree,studyProbes,r,seed,opts);
    probeError(j)=trapz(theta,abs(probeCurves(:,j)-jacksonTruth));
    probeMatvecs(j)=probeInfo{j}.matvecs;
    probeMeanSE(j)=mean(probeInfo{j}.standard_error);
end
probeSummary=table(string(probeChoices(:)),probeError, ...
    probeMatvecs,probeMeanSE,'VariableNames', ...
    {'ProbeMode','L1Error','Matvecs','MeanSE'});
fprintf('Moment window sweep (%d choices):\n',numWindows);
disp(windowSummary);
fprintf('Moment probe/trace sweep (%d choices, Jackson):\n', ...
    numProbeChoices);
disp(probeSummary);
if showPlots
    figure('Name','2. Moments: smoothing and trace options');
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile;
    plot(theta,rhoTrue,'k-','LineWidth',1.6);hold on
    plot(theta,rhoMoments,'--','LineWidth',1.3);
    hold off;grid on;xlim([-pi pi]);
    title('Poisson trace estimate (full K)');
    legend('Exact Poisson','Four-phase trace','Location','best');
    xlabel('\theta');ylabel('Density');
    nexttile;
    chosenWindows=[2 3 4 8]; % Jackson, Fejer, Dirichlet, heat alias
    plot(theta,windowTruth(:,chosenWindows),'LineWidth',1.2);
    grid on;xlim([-pi pi]);
    title('Different targets: exact finite window densities');
    legend(windowNames(chosenWindows),'Location','best');
    xlabel('\theta');ylabel('Density');
    nexttile;
    plot(theta,jacksonTruth,'k-','LineWidth',1.6);hold on
    plot(theta,probeCurves(:,1),'--','LineWidth',1.2);
    plot(theta,probeCurves(:,end),':','LineWidth',1.6);
    hold off;grid on;xlim([-pi pi]);
    title('Jackson: stochastic vs exact trace');
    legend('Exact Jackson','Four-phase','Exact basis','Location','best');
    xlabel('\theta');ylabel('Density');
    nexttile;
    bar(windowError);
    set(gca,'XTick',1:numWindows,'XTickLabel',windowNames);
    xtickangle(45);grid on;
    title('Window-specific L1 error');
    ylabel('L1 error');
end
%%
%[text] ## 3. Isometric Arnoldi / Gauss-Szego: repeated random starts
% The m-node unitary Hessenberg quadrature rule from each start v has
% nodes theta_j and nonuniform spectral weights |W(1,j)|^2. A weighted
% quadrature (never equally weighted reduced eigenvalues) matches
% v'*U^k*v for |k|<=m-1 in exact arithmetic.
%
% The Poisson quadrature is approximating an INFINITE Fourier series.
% The JACKSON quadrature with degree K<=m-1 is, for the same starting
% vectors, equivalent to stochastic Fourier moments up to roundoff.
% Repeat/average across randomized four-phase, Rademacher, complex
% Gaussian, all basis and custom starting vectors. Standard errors
% measure probe variability only, not deterministic quadrature error.
[rhoArnoldi,arnoldiInfo]=usd_szego_density( ...
    U,theta,studyOrder,exampleProbes,r,seed);
errArnoldi=trapz(theta,abs(rhoArnoldi-rhoTrue));
fprintf('Szego (Poisson): L1 %.4g, %d matvecs\n', ...
    errArnoldi,arnoldiInfo.matvecs);
startChoices={'fourphase','rademacher','complex_gaussian', ...
    'basis','custom'};
numStartChoices=numel(startChoices);
quadratureCurves=zeros(numel(theta),numStartChoices);
quadratureError=zeros(numStartChoices,1);
quadratureMatvecs=zeros(numStartChoices,1);
quadratureMeanSE=zeros(numStartChoices,1);
quadratureInfo=cell(numStartChoices,1);
for j=1:numStartChoices
    opts=struct('smoothing','jackson','degree',studyDegree);
    if strcmp(startChoices{j},'custom')
        opts.probe_vectors=full(speye(N,min(N,2)));
    else
        opts.probe_type=startChoices{j};
    end
    [quadratureCurves(:,j),quadratureInfo{j}]=usd_szego_density( ...
        U,theta,studyOrder,studyProbes,r,seed,opts);
    quadratureError(j)=trapz(theta, ...
        abs(quadratureCurves(:,j)-jacksonTruth));
    quadratureMatvecs(j)=quadratureInfo{j}.matvecs;
    quadratureMeanSE(j)=mean(quadratureInfo{j}.standard_error);
end
rhoArnoldiJackson=quadratureCurves(:,1);
momentQuadratureDiscrepancy= ...
    max(abs(probeCurves(:,1)-rhoArnoldiJackson));
quadratureSummary=table(string(startChoices(:)),quadratureError, ...
    quadratureMatvecs,quadratureMeanSE, ...
    'VariableNames',{'StartType','L1Error','Matvecs','MeanSE'});
fprintf('Arnoldi start-vector sweep (%d choices, Jackson):\n', ...
    numStartChoices);
disp(quadratureSummary);
fprintf('Same-start Jackson moments vs Szego (max difference): %.3g\n', ...
    momentQuadratureDiscrepancy);
% With the SAME probe vectors, Jackson-smoothed moment reconstruction and
% a Szego rule of order m >= K+1 should OVERLAP (moment exactness). Rather
% than plot those identical curves again, visualize quadrature's defining
% weighted spectral nodes and its distinct order-dependent Poisson error.
quadratureOrders=unique(max(1,round(studyOrder*[.20 .50 1])));
quadratureOrderCurves=zeros(numel(theta),numel(quadratureOrders));
quadratureOrderErrors=zeros(numel(quadratureOrders),1);
quadratureOrderMatvecs=zeros(numel(quadratureOrders),1);
for idx=1:numel(quadratureOrders)
    [quadratureOrderCurves(:,idx),orderInfo]=usd_szego_density( ...
        U,theta,quadratureOrders(idx),studyProbes,r,seed);
    quadratureOrderErrors(idx)=trapz(theta, ...
        abs(quadratureOrderCurves(:,idx)-rhoTrue));
    quadratureOrderMatvecs(idx)=orderInfo.matvecs;
end
quadProbe=usd_probe_vectors(N,1,struct( ...
    'probe_type','fourphase','seed',seed));
[quadratureNodes,quadratureWeights,quadratureRuleInfo]= ...
    usd_szego_rule(U,quadProbe,studyOrder);
fprintf('Poisson Szego order study (order / L1 / matvecs):\n');
disp(table(quadratureOrders(:),quadratureOrderErrors, ...
    quadratureOrderMatvecs,'VariableNames', ...
    {'Order','PoissonL1','Matvecs'}));
if showPlots
    figure('Name','3. Szego: nodes, weights, and convergence');
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile;
    stem(quadratureNodes,quadratureWeights,'filled','LineWidth',1.1);
    grid on;xlim([-pi pi]);
    title(sprintf('One %d-node Szego rule (nonuniform weights)', ...
        numel(quadratureNodes)));
    xlabel('Quadrature phase \theta');ylabel('Spectral weight');
    nexttile;
    plot(theta,rhoTrue,'k-','LineWidth',1.8);hold on
    plot(theta,quadratureOrderCurves,'LineWidth',1.2);hold off
    grid on;xlim([-pi pi]);
    title('Poisson quadrature: change Arnoldi order');
    legend([{'Exact Poisson'}, ...
        arrayfun(@(m)sprintf('order %d',m),quadratureOrders, ...
        'UniformOutput',false)],'Location','best');
    xlabel('\theta');ylabel('Density');
    nexttile;
    bar(quadratureError);grid on
    set(gca,'XTick',1:numStartChoices,'XTickLabel',startChoices);
    xtickangle(35);
    title('Random starting vectors, fixed Jackson window');
    ylabel('Jackson L1 error');
    nexttile;
    semilogy(quadratureOrders,quadratureOrderErrors,'o-', ...
        'LineWidth',1.5);grid on
    title('Poisson quadrature L1 error vs order');
    xlabel('Arnoldi order m');ylabel('Poisson L1 error');
end
%%
%[text] ## 4. Kernel methods: direct backslash, GMRES and sketched Krylov
% Direct Poisson filtering uses the Caratheodory/Herglotz transform:
%
%   z=r*exp(1i*theta),
%   rho_r(theta)=Re[1+(2*z/N)*tr((U-z*I)^(-1))]/(2*pi).
%
% For small N the default is direct backslash on the explicit matrix
% assembled once via N operator actions. GMRES is retained for larger
% matrices or explicit selection. Randomized sketches project a
% shared forward Krylov basis Q with Gaussian or sparse CountSketch S:
%
%   min_c ||S*((U-z*I)*Q*c-v)||_2.
%
% A true (UNSKETCHED) residual is reported for each shift. The option
% sketch_fallback=true retries when that residual is too large.
% Experiments deliberately show BOTH unchecked and fallback versions;
% an unchecked sketch may have substantial systematic error or negative
% apparent density. Sketched forward Krylov is an experimental baseline,
% not an optimized low-memory sGMRES implementation.
%
% To control the interactive cost, GMRES/sketches use a coarser angle grid;
% their plotted/reference densities are compared on THAT same grid.
[rhoKernel,kernelInfo]=usd_poisson_density( ...
    U,theta,exampleProbes,r,seed);
errKernel=trapz(theta,abs(rhoKernel-rhoTrue));
fprintf('Main resolvent (%s): L1 %.4g, %d matvecs, residual %.3g\n', ...
    kernelInfo.solver,errKernel,kernelInfo.matvecs, ...
    kernelInfo.max_relative_residual);
% Explicitly exercise matrix-free GMRES even when the small-problem
% default is direct backslash. Both use the same seeded trace probes.
[rhoGMRES,gmresInfo]=usd_poisson_density( ...
    U,kernelTheta,studyProbes,r,seed,1e-9,min(40,N),[], ...
    struct('solver','gmres'));
gmresError=trapz(kernelTheta,abs(rhoGMRES-kernelReference));
sketchNames={'gaussian','countsketch', ...
    'gaussian + fallback','countsketch + fallback'};
sketchTypes={'gaussian','countsketch','gaussian','countsketch'};
sketchFallback=[false false true true];
numSketches=numel(sketchNames);
sketchCurves=zeros(numel(kernelTheta),numSketches);
sketchError=zeros(numSketches,1);
sketchMatvecs=zeros(numSketches,1);
sketchMaxResidual=zeros(numSketches,1);
sketchFallbackSolves=zeros(numSketches,1);
sketchMeanSE=zeros(numSketches,1);
sketchInfo=cell(numSketches,1);
for j=1:numSketches
    opts=struct('solver','sketch','sketch_type',sketchTypes{j}, ...
        'sketch_order',sketchOrder, ...
        'sketch_dimension',sketchDimension, ...
        'sketch_fallback',sketchFallback(j), ...
        'sketch_residual_tolerance',1e-5);
    [sketchCurves(:,j),sketchInfo{j}]=usd_poisson_density( ...
        U,kernelTheta,studyProbes,r,seed, ...
        1e-9,min(40,N),[],opts);
    sketchError(j)=trapz(kernelTheta, ...
        abs(sketchCurves(:,j)-kernelReference));
    sketchMatvecs(j)=sketchInfo{j}.matvecs;
    sketchMaxResidual(j)=sketchInfo{j}.max_relative_residual;
    sketchFallbackSolves(j)=sketchInfo{j}.fallback_solves;
    sketchMeanSE(j)=mean(sketchInfo{j}.standard_error);
end
sketchSummary=table(string(sketchNames(:)),sketchError,sketchMatvecs, ...
    sketchMaxResidual,sketchFallbackSolves,sketchMeanSE, ...
    'VariableNames',{'Solver','L1Error','Matvecs', ...
    'MaxTrueResidual','FallbackSolves','MeanSE'});
fprintf('Randomized projected shifted-system sweep:\n');
disp(sketchSummary);
if showPlots
    figure('Name','4. Resolvent: backslash, GMRES and sketches');
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile;
    plot(theta,rhoTrue,'k-','LineWidth',1.6);hold on
    plot(theta,rhoKernel,'--','LineWidth',1.3);hold off
    grid on;xlim([-pi pi]);
    title(sprintf('%s baseline on full phase grid',kernelInfo.solver));
    legend('Exact Poisson',kernelInfo.solver,'Location','best');
    xlabel('\theta');ylabel('Density');
    nexttile;
    plot(kernelTheta,kernelReference,'k-','LineWidth',1.6);hold on
    plot(kernelTheta,rhoGMRES,'--','LineWidth',1.1);
    plot(kernelTheta,sketchCurves(:,2),':','LineWidth',1.3);
    plot(kernelTheta,sketchCurves(:,4),'-.','LineWidth',1.3);
    hold off;grid on;xlim([-pi pi]);
    title('Gaussian, CountSketch, fallback');
    legend('Exact Poisson','GMRES','CountSketch', ...
        'CountSketch + fallback','Location','best');
    xlabel('\theta');ylabel('Density');
    nexttile;
    semilogy(kernelTheta,max(sketchInfo{1}.relative_residuals,[],2), ...
        'LineWidth',1.2);hold on
    semilogy(kernelTheta,max(sketchInfo{2}.relative_residuals,[],2), ...
        '--','LineWidth',1.2);
    semilogy(kernelTheta,max(sketchInfo{4}.relative_residuals,[],2), ...
        ':','LineWidth',1.6);
    hold off;grid on;xlim([-pi pi]);
    title('True solve residual by phase');
    legend('Gaussian','CountSketch','CountSketch + fallback', ...
        'Location','best');
    xlabel('\theta');ylabel('Relative residual');
    nexttile;
    scatter(sketchMatvecs,sketchError,55,'filled');hold on
    for j=1:numSketches
        text(sketchMatvecs(j),sketchError(j),[' ' sketchNames{j}]);
    end
    hold off;grid on;
    title('Sketched resolvent: error vs matvecs');
    xlabel('U matvecs');ylabel('Poisson L1 error');
end
%%
%[text] ## 5. Compare methods, options, errors, work, and target measures
% The three baseline curves below ALL target the same exact FINITE
% Poisson-smoothed trace measure. Their different matvec budgets must
% NOT be interpreted as a fair complexity comparison.
%
% The optionSummary table also compares all window, starting-vector and
% sketch variants. Each row carries its target: a finite Poisson measure
% versus the corresponding finite WINDOW measure. These heterogeneous
% targets must not be ranked as though they were the same quantity.
% All computed arrays, diagnostics, and tables remain in the workspace.
method={'Fourier moments';'Arnoldi / Gauss-Szego';'Poisson resolvent'};
L1Error=[errMoments;errArnoldi;errKernel];
Matvecs=[momentsInfo.matvecs;arnoldiInfo.matvecs;kernelInfo.matvecs];
comparison=table(method,L1Error,Matvecs);
fprintf('\nThree main Poisson comparisons:\n');
disp(comparison);
Option=[ "Moment window: " + string(windowNames(:)); ...
    "Moment trace: " + string(probeChoices(:)); ...
    "Szego start: " + string(startChoices(:)); ...
    "Resolvent: " + string(kernelInfo.solver); ...
    "Resolvent: GMRES"; ...
    "Resolvent sketch: " + string(sketchNames(:)) ];
Target=[repmat("Window-specific finite",numWindows,1); ...
    repmat("Jackson finite",numProbeChoices+numStartChoices,1); ...
    repmat("Poisson finite",numSketches+2,1)];
OptionL1=[windowError;probeError;quadratureError; ...
    errKernel;gmresError;sketchError];
OptionMatvecs=[windowMatvecs;probeMatvecs;quadratureMatvecs; ...
    kernelInfo.matvecs;gmresInfo.matvecs;sketchMatvecs];
OptionMeanSE=[windowMeanSE;probeMeanSE;quadratureMeanSE; ...
    mean(kernelInfo.standard_error); ...
    mean(gmresInfo.standard_error);sketchMeanSE];
OptionResidual=[nan(numWindows+numProbeChoices+numStartChoices,1); ...
    kernelInfo.max_relative_residual;gmresInfo.max_relative_residual; ...
    sketchMaxResidual];
optionComparison=table(Option,Target,OptionL1,OptionMatvecs, ...
    OptionMeanSE,OptionResidual);
fprintf('All %d experimental choices (targets identified):\n', ...
    height(optionComparison));
disp(optionComparison);
if showPlots
    figure('Name','5. Dirac density: all method families');
    tiledlayout(2,2,'TileSpacing','compact');
    nexttile([1 2]);
    plot(theta,rhoTrue,'k-','LineWidth',2);hold on
    plot(theta,rhoMoments,'--','LineWidth',1.3);
    plot(theta,rhoArnoldi,'-.','LineWidth',1.3);
    plot(theta,rhoKernel,':','LineWidth',1.8);hold off
    grid on;xlim([-pi pi]);xlabel('Eigenphase \theta (radians)');
    ylabel('Poisson-smoothed finite trace density');
    title(sprintf('Three estimators vs exact Dirac measure (r=%.2f)',r));
    legend('Exact Poisson','Fourier moments', ...
        'Szego quadrature','Direct resolvent','Location','best');
    nexttile;
    bar(L1Error);grid on
    set(gca,'XTick',1:3,'XTickLabel', ...
        {'Moments','Szego','Direct'});
    title('Poisson L1 errors');
    ylabel('L1 error');
    nexttile;
    bar(Matvecs);grid on
    set(gca,'XTick',1:3,'XTickLabel', ...
        {'Moments','Szego','GMRES'});
    title('U matvec budget (not equal)');
    ylabel('Matvecs');
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
