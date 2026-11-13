# Unitary Spectral Density

Portable MATLAB experiments for approximating the **normalized empirical
eigenphase measure** of a unitary matrix from matrix-vector products.
The pinned public `anymatrix-unitary` library supplies test matrices,
not the numerical estimators.

## First benchmark: two-arc Dirac quantum cellular automaton

`dirac_automaton(n,m,'operator')` has order `2*n` and eigenphases
`+/-acos(sqrt(1-m^2)*cos(2*pi*k/n))`. For `0<m<1` these lie in
`+/-[asin(m),pi-asin(m)]`. The exact *finite* eigenphases are used to
validate the Poisson-smoothed trace density. The limiting two-band density
is a separate reference, not the finite empirical measure.

## Three independent computational strategies

| Folder | Estimator | Main routine |
|---|---|---|
| `moments/` | Fourier trace moments and selectable circle smoothing windows | `usd_moment_density` |
| `quadrature/` | Random-start isometric Arnoldi and unitary Szegő quadrature | `usd_szego_density` |
| `filters/` | Poisson-resolvent traces via GMRES or sketched Krylov solves | `usd_poisson_density` |

Estimators use only `U(X,'notransp')`, `U([],'size')` and MATLAB built-ins.
They do **not** import Dirac-specific spectral formulas or Anymatrix functions.
All three default to seeded four-phase trace probes; randomized probe
distributions and canonical-basis references are now configurable.
A common Poisson radius gives comparable smoothing but **not equal
matvec budgets**. Each routine reports its matvec count.

`quadrature/` matches moments of each random vector spectral measure.
The **Gauss–Szegő quadrature weights must be retained**; the equal-weight
eigenvalues of the reduced matrix are not a trace-density estimator.
`filters/` solves `(U-r*exp(1i*theta)*I)y=v` with GMRES at each angle.

## Experimental controls

All functions retain their original positional arguments; pass a MATLAB
options struct as the final argument. Consult
[mathematics.md](mathematics.md) for equations, error sources, literature
references, and distinctions between established unitary methods and
experimental Hermitian-to-unitary adaptations.

**Moment smoothing:** choose Poisson (default), Fejer, Jackson,
Dirichlet, Lanczos sigma, heat/wrapped Gaussian, or explicit
Poisson-Fejer / Poisson-Jackson compositions. Jackson uses the
*circle squared-Fejer kernel*, not the real-line Chebyshev KPM damping
coefficients. Non-Poisson windows ignore the positional radius r.

~~~matlab
opts = struct('smoothing','jackson', ...
              'trace_mode','stochastic','probe_type','fourphase');
[rho,info] = usd_moment_density(U,theta,40,20,.8,17,opts);

opts.trace_mode = 'exact';  % complete canonical basis, n*K matvecs
[rhoExact,exactInfo] = usd_moment_density(U,theta,40,1,.8,17,opts);
~~~

**Arnoldi / Gauss-Szego:** repeatedly sample independent starting
vectors, form nonuniform-weight quadrature measures, and average.
Select the same finite smoothing kernel and degree as moments
(with degree <= order-1) for a mathematically matched experiment.

~~~matlab
opts = struct('smoothing','jackson','degree',30, ...
              'probe_type','complex_gaussian');
[rho,info] = usd_szego_density(U,theta,31,20,.8,17,opts);
% info.standard_error and info.effective_orders are per-run diagnostics.
~~~

**Poisson resolvents:** use classic GMRES (default) or a randomized
sketch of Krylov projected shifted systems. Gaussian and sparse
CountSketch embeddings are supported. The sketch is independently
seeded, reused over phase shifts, and checked using TRUE residuals.
Optional fallback to GMRES is deliberately separate from the
projected method, so approximate-sketch experiments remain possible.

~~~matlab
opts = struct('solver','sketch','sketch_order',40, ...
              'sketch_dimension',120,'sketch_type','countsketch', ...
              'sketch_fallback',true, ...
              'sketch_residual_tolerance',1e-5);
[rho,info] = usd_poisson_density(U,theta,8,.8,17,1e-9,40,[],opts);
disp([info.matvecs,info.max_relative_residual,info.fallback_solves]);
~~~

Available trace-probe choices are fourphase (default), rademacher,
complex_gaussian, and basis; a custom n-by-s matrix can be passed
using opts.probe_vectors. Random vectors are normalized and isotropic.
A complete basis gives an exact trace, while arbitrary custom starts
are not necessarily unbiased.

See [examples/experimental_options.m](examples/experimental_options.m)
for a plain-text comparative script with error diagnostics. The existing
five-section Dirac Live Editor script remains backward compatible.

## Sectioned MATLAB Live Editor example

Open [`examples/dirac_spectral_density.m`](examples/dirac_spectral_density.m)
in MATLAB's Live Editor or ordinary Editor. It is a **plain-text `.m` file**
with five `%%` executable sections: Dirac introduction and reference density,
Fourier moments, Gauss–Szegő/Arnoldi, Poisson-resolvent filtering, and a
four-curve comparison. Each method has its own plot; the exact reference
is the **Poisson-smoothed finite trace measure**, not the singular limiting
density.

Run all five sections in sequence, or run the file as a script. Figures are
enabled by default. Set `showPlots=false` beforehand for batch execution.
The example accepts optional `exampleSites`, `exampleGridPoints`,
`exampleProbes`, `exampleRadius`, `exampleMomentDegree`, and
`exampleArnoldiOrder` variables.

The CI suite includes a reduced-size smoke test of all five sections.

## Install and run

```bash
git clone --recurse-submodules https://github.com/jaurentz/unitary-spectral-density.git
cd unitary-spectral-density
```

From MATLAB R2020b or later (no additional toolboxes):

```matlab
setup;
out=run_all_dirac(true);       % all three, plots and errors
a=run_dirac_moments([],false);  % individual experiment
b=run_dirac_quadrature([],false);
c=run_dirac_filters([],false);
results=runtests('tests');
assertSuccess(results);
```

Individual drivers accept an optional `p=dirac_case(...)` parameter and
`showPlot` flag; e.g. `run_dirac_filters(dirac_case(16,.5,.7,33),false)`.
Shared probe generation and circular Poisson kernels live in `shared/`.
Benchmark-only formulas live in `benchmarks/`.
No figures or large experimental outputs are committed by default.

## Git submodule

The `external/anymatrix-unitary` Gitlink pins a tested public commit.
To advance deliberately:

```bash
git submodule update --init --recursive
git submodule update --remote external/anymatrix-unitary
git add external/anymatrix-unitary
git commit -m "Update unitary matrix dependency"
```

## Continuous integration and parallel development

The MATLAB workflow tests **both** the recorded Anymatrix Git submodule
commit and a separately resolved latest `anymatrix-unitary/main` commit
on pushes, PRs, manual runs and a daily schedule. A cross-repository
dispatch can additionally trigger integration tests whenever upstream
MATLAB CI passes. After both modes pass, the scheduled/dispatched run
can propose an updated submodule pin by pull request.

See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) for full details,
required one-time GitHub permissions, and the local two-clone workflow.

## GitHub Actions / MATLAB regression

`.github/workflows/matlab-tests.yml` uses MathWorks' `setup-matlab@v3`
and `run-command@v3`, checking out the pinned submodule and running
the full estimator suite on push, pull request, and manual dispatch.

The repository is public, so MATLAB Actions can run on GitHub-hosted runners
without an additional MATLAB licensing token. The CI workflow runs on each
push to `main`, on pull requests, and by manual dispatch.

Source: https://github.com/matlab-actions/setup-matlab

## Limitations

The quadrature prototype uses twice-reorthogonalized Arnoldi and dense
unitary closure, not a specialized short-recurrence isometric Arnoldi
algorithm. The sketched solver also builds and stores a full Krylov
basis; it is not a low-memory or optimized sketched GMRES implementation.
For a general unitary spectrum, forward Krylov spaces may approximate
interior-disk resolvents poorly. Check actual residuals and enable
fallback when appropriate. Truncated Poisson and heat reconstructions
can be slightly negative even though their untruncated kernels are positive.

MATLAB tests are configured; inspect the latest Actions run for results.
