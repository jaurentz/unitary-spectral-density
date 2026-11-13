function out=run_dirac_moments(p,showPlot)
%RUN_DIRAC_MOMENTS Fourier/trace moment experiment on Dirac bands.
if nargin<1 || isempty(p),p=dirac_case;end
if nargin<2,showPlot=true;end
degree=80;probes=8;seed=17;
[rho,info]=usd_moment_density(p.U,p.theta,degree,probes,p.radius,seed);
out=struct('rho',rho,'info',info, ...
    'L1_error',trapz(p.theta,abs(rho-p.exact_smoothed)));
if showPlot
    figure;plot(p.theta,p.exact_smoothed,'k-',p.theta,rho,'--','LineWidth',1.4);
    xlabel('Eigenphase (rad)');ylabel('Poisson density');
    title('Dirac: stochastic trace moments');legend('Exact smoothed','Moments');
end
end
