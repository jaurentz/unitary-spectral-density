function out=run_dirac_quadrature(p,showPlot)
%RUN_DIRAC_QUADRATURE Isometric Arnoldi/Gauss-Szego Dirac experiment.
if nargin<1 || isempty(p),p=dirac_case;end
if nargin<2,showPlot=true;end
order=min(32,p.dimension);probes=8;seed=17;
[rho,info]=usd_szego_density(p.U,p.theta,order,probes,p.radius,seed);
out=struct('rho',rho,'info',info, ...
    'L1_error',trapz(p.theta,abs(rho-p.exact_smoothed)));
if showPlot
    figure;plot(p.theta,p.exact_smoothed,'k-',p.theta,rho,'--','LineWidth',1.4);
    xlabel('Eigenphase (rad)');ylabel('Poisson density');
    title('Dirac: isometric Arnoldi / Gauss-Szego');
    legend('Exact smoothed','Szego quadrature');
end
end
