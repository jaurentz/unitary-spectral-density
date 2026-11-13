function root=setup
%SETUP Add the three portable estimators and pinned Anymatrix dependency.
root=fileparts(mfilename('fullpath'));
dep=fullfile(root,'external','anymatrix-unitary');
if exist(fullfile(dep,'dirac_automaton.m'),'file')~=2
    error('usd:Dependency', ...
        'Missing Git submodule. Run git submodule update --init --recursive.');
end
addpath(dep,fullfile(root,'shared'),fullfile(root,'moments'), ...
    fullfile(root,'quadrature'),fullfile(root,'filters'), ...
    fullfile(root,'benchmarks'));
end
