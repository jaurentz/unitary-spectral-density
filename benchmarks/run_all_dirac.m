function out=run_all_dirac(showPlot)
%RUN_ALL_DIRAC Compare all three families with the same Dirac operator.
if nargin<1,showPlot=true;end
setup;
p=dirac_case;
out=struct('problem',rmfield(p,'U'), ...
    'moments',run_dirac_moments(p,false), ...
    'quadrature',run_dirac_quadrature(p,false), ...
    'filters',run_dirac_filters(p,false));
fprintf('Dirac order %d, mass %.3f, Poisson radius %.3f\n', ...
    p.dimension,p.mass,p.radius);
fprintf('%-14s  %-12s  %s\n','Estimator','L1 error','Matvec actions');
names={'moments','quadrature','filters'};
for k=1:numel(names)
    a=out.(names{k});
    fprintf('%-14s  %12.4g  %d\n',names{k},a.L1_error,a.info.matvecs);
end
if showPlot
    figure;plot(p.theta,p.exact_smoothed,'k-','LineWidth',2);hold on;
    plot(p.theta,out.moments.rho,'--','LineWidth',1.2);
    plot(p.theta,out.quadrature.rho,'-.','LineWidth',1.2);
    plot(p.theta,out.filters.rho,':','LineWidth',1.5);
    xlabel('Eigenphase (rad)');ylabel('Smoothed trace density');
    title(sprintf('Dirac automaton: n=%d, m=%g, r=%g', ...
        p.dimension,p.mass,p.radius));
    legend('Exact finite','Moments','Gauss-Szego','Resolvent');
    xlim([-pi pi]);grid on;
end
end
