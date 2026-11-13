# Unitary Spectral Density (private research repository)

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
| `moments/` | Stochastic Fourier traces, truncated Poisson expansion | `usd_moment_density` |
| `quadrature/` | Isometric Arnoldi and unitary Szegő quadrature | `usd_szego_density` |
| `filters/` | Direct Poisson-resolvent traces via shifted GMRES solves | `usd_poisson_density` |

Estimators use only `U(X,'notransp')`, `U([],'size')` and MATLAB built-ins.
They do **not** import Dirac-specific spectral formulas or Anymatrix functions.
All three use independent seeded four-phase random vectors and a common
Poisson radius `0<r<1`. This gives comparable smoothing but **not equal
matvec budgets**. Each routine reports its matvec count.

`quadrature/` matches moments of each random vector spectral measure.
The **Gauss–Szegő quadrature weights must be retained**; the equal-weight
eigenvalues of the reduced matrix are not a trace-density estimator.
`filters/` solves `(U-r*exp(1i*theta)*I)y=v` with GMRES at each angle.

## Install and run

```bash
git clone --recurse-submodules https://github.com/jaurentz/unitary-spectral-density.git
cd unitary-spectral-density
```

From MATLAB R2020b or later (no additional toolboxes):

```matlab
setup;
out=run_all_dirac(true);       % all three, plots and errors
a=run_dirac_moments(false);   % individual experiment
b=run_dirac_quadrature(false);
c=run_dirac_filters(false);
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

## GitHub Actions and MATLAB licensing

`.github/workflows/matlab-tests.yml` uses MathWorks' `setup-matlab@v3`
and `run-command@v3`, checking out the pinned submodule and running
the full estimator suite on push, pull request, and manual dispatch.

**MATLAB Actions do not automatically license MATLAB for a private repository.**
Create the GitHub Actions repository secret `MLM_LICENSE_TOKEN`, containing
a valid MATLAB batch licensing token, before expecting the workflow to pass.
The workflow deliberately fails with a clear message when the secret is absent.
An appropriately licensed self-hosted runner is an alternative.

Source: https://github.com/matlab-actions/setup-matlab

## Limitations

The initial quadrature prototype uses twice-reorthogonalized Arnoldi
and a dense order-`m` unitary closure; it is not yet a specialized
short-recurrence isometric Arnoldi implementation. The direct-filter
method can require many shifted GMRES iterations as `r` approaches one.
Moment truncation can produce a slightly negative reconstruction even
though the untruncated Poisson density is positive.

MATLAB tests are configured but not claimed to have passed until a
licensed CI run has completed.
