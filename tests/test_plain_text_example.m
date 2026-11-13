function tests=test_plain_text_example
%TEST_PLAIN_TEXT_EXAMPLE Compact five-section Dirac plain-text Live Script.
tests=functiontests(localfunctions);
end

function test_five_sections_and_spectral_reference(tc)
root=fileparts(fileparts(mfilename('fullpath')));
script=fullfile(root,'examples','dirac_spectral_density.m');
source=fileread(script);
verifyEqual(tc,numel(regexp(source,'(?m)^%%$','match')),5);
verifyEqual(tc,numel(regexp(source,'(?m)^%\[text\] ## [1-5]\. ','match')),5);
verifyEqual(tc,numel(strfind(source,'%[appendix]{"version":"1.0"}')),1);
verifyTrue(tc,contains(source,'dirac_automaton(sites,mass,''spectral'')'));
verifyTrue(tc,contains(source,'dirac_automaton(sites,mass,''operator'')'));
verifyTrue(tc,contains(source,'S.finite.eigenvalues()'));
verifyFalse(tc,contains(source,'dirac_case('));
verifyFalse(tc,contains(source,'windowSummary'));
verifyLessThan(tc,numel(splitlines(string(source))),100);
end

function test_run_example(tc)
root=fileparts(fileparts(mfilename('fullpath')));
script=fullfile(root,'examples','dirac_spectral_density.m');
prior=get(groot,'DefaultFigureVisible');
cleanup=onCleanup(@()restore_visibility(prior)); %#ok<NASGU>
set(groot,'DefaultFigureVisible','off');
close all force
run(script);
verifyEqual(tc,height(comparison),3);
verifyEqual(tc,numel(rhoTrue),121);
verifyEqual(tc,S.dimension,96);
verifyEqual(tc,S.finite.kind,'exact');
verifyTrue(tc,all(isfinite([rhoTrue rhoMoments rhoArnoldi rhoKernel]),'all'));
verifyGreaterThan(tc,min(Matvecs),0);
verifyLessThan(tc,max(L1Error),.65);
verifyEqual(tc,kernelInfo.solver,'backslash');
end

function test_benchmark_uses_spectral_descriptor(tc)
p=dirac_case(16,.45,.7,41);
verifyEqual(tc,p.spectral.family,'dirac_automaton');
verifyEqual(tc,p.dimension,32);
verifyEqual(tc,numel(p.exact_phases),p.dimension);
verifyEqual(tc,p.exact_smoothed,usd_poisson_kernel( ...
    p.theta,p.exact_phases,ones(32,1)/32,p.radius),'AbsTol',1e-13);
end

function restore_visibility(prior)
close all force
set(groot,'DefaultFigureVisible',prior);
end
