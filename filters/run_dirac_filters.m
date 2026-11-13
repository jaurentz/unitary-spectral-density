function out=run_dirac_filters(p,showPlot)
%RUN_DIRAC_FILTERS Direct resolvent/Poisson Dirac experiment (GMRES).
if nargin<1 || isempty(p),p=dirac_case;end
if nargin<2,showPlot=true;end
probes=8;seed=17;
[rho,info]=usd_poisson_density(p.U,p.theta,probes,p.radius,seed);
out=struct('rho',rho,'info',info, ...
    'L1_error',trapz(p.theta,abs(rho-p.exact_smoothed)));
if showPlot
    figure;plot(p.theta,p.exact_smoothed,'k-',p.theta,rho,'--','LineWidth',1.4);
    xlabel('Eigenphase (rad)');ylabel('Poisson density');
    title('Dirac: direct Poisson resolvent');legend('Exact smoothed','Resolvent');
end
end
