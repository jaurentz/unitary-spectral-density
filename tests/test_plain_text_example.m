function tests = test_plain_text_example
%TEST_PLAIN_TEXT_EXAMPLE Run the sectioned Dirac .m example without figures.
tests = functiontests(localfunctions);
end

function test_all_five_sections(tc)
root = fileparts(fileparts(mfilename('fullpath')));
script = fullfile(root,'examples','dirac_spectral_density.m');
source = fileread(script);
sections = regexp(source,'(?m)^%% [1-5]\. ','match');
verifyEqual(tc,numel(sections),5);

% Smaller settings make the CI check quick; the Live Editor defaults
% remain higher-resolution for interactive use.
showPlots = false; %#ok<NASGU>
exampleSites = 12; %#ok<NASGU>
exampleGridPoints = 33; %#ok<NASGU>
exampleProbes = 2; %#ok<NASGU>
exampleRadius = .70; %#ok<NASGU>
exampleMomentDegree = 32; %#ok<NASGU>
exampleArnoldiOrder = 12; %#ok<NASGU>
run(script);

verifyEqual(tc,height(comparison),3);
verifyEqual(tc,numel(rhoTrue),exampleGridPoints);
verifyTrue(tc,all(isfinite([rhoTrue rhoMoments rhoArnoldi rhoKernel]),'all'));
verifyGreaterThan(tc,min([momentsInfo.matvecs, ...
    arnoldiInfo.matvecs,kernelInfo.matvecs]),0);
verifyLessThan(tc,max(comparison.L1Error),.7);

% Every implemented window, trace mode, Arnoldi start, and sketch option
% must be exercised in the same plain-text five-section example.
verifyEqual(tc,numel(windowNames),10);
verifyEqual(tc,numel(probeChoices),6);
verifyEqual(tc,numel(startChoices),5);
verifyEqual(tc,numel(sketchNames),4);
verifyEqual(tc,height(optionComparison),26);
verifyEqual(tc,height(windowSummary),10);
verifyEqual(tc,height(probeSummary),6);
verifyEqual(tc,height(quadratureSummary),5);
verifyEqual(tc,height(sketchSummary),4);
verifyTrue(tc,all(isfinite(optionComparison.OptionL1)));
verifyGreaterThanOrEqual(tc,min(windowMatvecs),1);
verifyGreaterThan(tc,probeInfo{end}.matvecs, ...
    probeInfo{1}.matvecs);
verifyLessThan(tc,momentQuadratureDiscrepancy,2e-8);
verifyEqual(tc,sketchInfo{1}.solver,'sketch');
verifyEqual(tc,sketchInfo{2}.sketch_type,'countsketch');
verifyGreaterThanOrEqual(tc,sketchInfo{4}.fallback_solves,0);
verifyEqual(tc,numel(kernelReference),numel(kernelTheta));
end

function test_all_five_figures(tc)
root = fileparts(fileparts(mfilename('fullpath')));
script = fullfile(root,'examples','dirac_spectral_density.m');
prior = get(groot,'DefaultFigureVisible');
cleanup = onCleanup(@()restore_visibility(prior)); %#ok<NASGU>
set(groot,'DefaultFigureVisible','off');
close all force

showPlots = true; %#ok<NASGU>
exampleSites = 8; %#ok<NASGU>
exampleGridPoints = 21; %#ok<NASGU>
exampleProbes = 1; %#ok<NASGU>
exampleRadius = .6; %#ok<NASGU>
exampleMomentDegree = 25; %#ok<NASGU>
exampleArnoldiOrder = 8; %#ok<NASGU>
run(script);

figures = findall(groot,'Type','figure');
verifyEqual(tc,numel(figures),5);
verifyTrue(tc,all(isfinite([rhoTrue rhoMoments rhoArnoldi rhoKernel]),'all'));
verifyEqual(tc,height(optionComparison),26);
end

function restore_visibility(value)
close all force
set(groot,'DefaultFigureVisible',value);
end
