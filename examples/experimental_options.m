%% 1. Dirac matrix-free benchmark and exact finite phase reference
% EXPERIMENTAL_OPTIONS Additional plain-text experiments for MATLAB.
% This intentionally does NOT modify the existing five-section Live script.
% Every estimator receives only U(X,'notransp') and U([],'size').
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root);setup;
if ~exist('showPlots','var'),showPlots=true;end
p=dirac_case(24,.45,.72,91);
U=p.U; theta=p.theta;r=p.radius;
n=p.dimension;s=6;K=22;m=K+1;seed=19;

%% 2. Circle Jackson damping and exact trace reference
% Compare the *same Jackson kernel*, not the Poisson reference.
weights=ones(n,1)/n;
exactJackson=usd_circular_density(theta,p.exact_phases, ...
    weights,r,'jackson',K,struct);
momOpts=struct('smoothing','jackson','probe_type','fourphase');
[rhoMom,im]=usd_moment_density(U,theta,K,s,r,seed,momOpts);
exactOpts=struct('smoothing','jackson','trace_mode','exact');
[rhoBasis,ie]=usd_moment_density(U,theta,K,1,r,seed,exactOpts);
fprintf('Jackson moments: random L1 %.3g; exact-basis L1 %.3g\n', ...
    trapz(theta,abs(rhoMom-exactJackson)), ...
    trapz(theta,abs(rhoBasis-exactJackson)));
if showPlots
    figure('Name','Circle Jackson smoothing');
    plot(theta,exactJackson,'k','LineWidth',2);hold on;
    plot(theta,rhoMom,'--');plot(theta,rhoBasis,':');hold off;
    legend('Exact finite Jackson','Random trace','Exact basis');
    xlabel('Phase');ylabel('Density');grid on
end

%% 3. Random-start Szego quadrature with matched Jackson degree
szOpts=struct('smoothing','jackson','degree',K, ...
    'probe_type','fourphase');
[rhoSzego,iq]=usd_szego_density(U,theta,m,s,r,seed,szOpts);
fprintf('Jackson Szego: L1 %.3g; %d matvecs\n', ...
    trapz(theta,abs(rhoSzego-exactJackson)),iq.matvecs);

%% 4. Direct Poisson resolvent and sketched projected solves
% Kernel methods target Poisson, not Jackson. Compare each to the
% correct reference and inspect true (unsketched) residuals.
poissonTruth=p.exact_smoothed;
[rhoGMRES,ig]=usd_poisson_density(U,theta,s,r,seed, ...
    1e-9,min(40,n),[]);
skOpts=struct('solver','sketch','sketch_order',min(36,n), ...
    'sketch_dimension',120,'sketch_type','countsketch', ...
    'sketch_fallback',true,'sketch_residual_tolerance',1e-5);
[rhoSketch,is]=usd_poisson_density(U,theta,s,r,seed, ...
    1e-9,min(40,n),[],skOpts);
fprintf('GMRES:  %d matvecs, true residual %.3g\n', ...
    ig.matvecs,ig.max_relative_residual);
fprintf('Sketch: %d matvecs, true residual %.3g, fallbacks %d\n', ...
    is.matvecs,is.max_relative_residual,is.fallback_solves);
if showPlots
    figure('Name','Resolvent sketch experiment');
    plot(theta,poissonTruth,'k','LineWidth',2);hold on;
    plot(theta,rhoGMRES,'--');plot(theta,rhoSketch,':');hold off;
    legend('Exact finite Poisson','GMRES','Sketched projection');
    xlabel('Phase');ylabel('Density');grid on
end

%% 5. Matvec budgets and errors against matched targets
method={'Random Jackson moments';'Exact Jackson moments'; ...
    'Szego Jackson';'Poisson GMRES';'Poisson sketch'};
L1=[trapz(theta,abs(rhoMom-exactJackson)); ...
    trapz(theta,abs(rhoBasis-exactJackson)); ...
    trapz(theta,abs(rhoSzego-exactJackson)); ...
    trapz(theta,abs(rhoGMRES-poissonTruth)); ...
    trapz(theta,abs(rhoSketch-poissonTruth))];
Matvecs=[im.matvecs;ie.matvecs;iq.matvecs;ig.matvecs;is.matvecs];
experimentComparison=table(method,L1,Matvecs);
disp(experimentComparison);
