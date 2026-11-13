# Parallel development and dependency testing

The numerical matrix library [`jaurentz/anymatrix-unitary`](https://github.com/jaurentz/anymatrix-unitary)
and [`jaurentz/unitary-spectral-density`](https://github.com/jaurentz/unitary-spectral-density)
are **separate, actively developed public repositories**.

## Two MATLAB test modes

[`.github/workflows/matlab-tests.yml`](../.github/workflows/matlab-tests.yml)
runs the **same complete MATLAB regression suite** in two independent jobs:

- **`anymatrix-pinned`**: the `external/anymatrix-unitary` Git submodule
  at the SHA recorded in the spectral-density commit. This supports
  reproducible benchmarks and catches regressions in density code.
- **`anymatrix-latest`**: a specific upstream SHA resolved once at the
  start of the workflow. This catches changes to matrix generators,
  operator handles, and spectral references before advancing the submodule.

Both MATLAB jobs print the numerical research commit and matrix-library
commit in the GitHub Actions run summary. `fail-fast: false` ensures a
failure in one dependency mode does not hide the other result.

The workflow runs on spectral-density pushes/PRs, manually, daily at
**03:23 UTC**, and on cross-repository `repository_dispatch` with event
type `anymatrix-main-passed`. On dispatch, the supplied upstream SHA is
tested, rather than silently substituting an unrelated later commit.

The additional `test_matrix_library_integration.m` test validates FFT
circulants, sparse Dirac quantum automata, and structured CMV operators,
including their forward/adjoint matrix-free actions, trace moments,
weighted Szegő quadrature, and Poisson-resolvent reconstruction against
finite-dimensional reference data.

## Trigger after a successful Anymatrix CI run

The upstream `anymatrix-unitary` workflow has a separate job that
dispatches an event to spectral density only after its **entire standalone
MATLAB test job succeeds on a push to `main`**. A regular
`GITHUB_TOKEN` cannot dispatch to a different repository; this needs a
separately authorized token.

1. Create a **fine-grained GitHub personal access token** limited to the
   `jaurentz/unitary-spectral-density` repository, with **Contents:
   Read and write** permission (and the default Metadata: Read).
2. Add it to the **upstream** repository's Actions secrets as
   `SPECTRAL_DENSITY_DISPATCH_TOKEN`:
   [`anymatrix-unitary` secrets](https://github.com/jaurentz/anymatrix-unitary/settings/secrets/actions).
3. The upstream job will then issue GitHub's `repository_dispatch` request
   with the successful source SHA. Until configured, dispatch reports an
   explicit notice and skips, while the daily downstream fallback still runs.

Do **not** commit a token into either repository. Keep the token scoped to
the single downstream repository; there is no need for administration,
other repository access, or general GitHub API privileges.

## Automatically propose a tested submodule advance

On successful **scheduled or dispatched** runs, the downstream workflow
attempts to open/update a pull request on branch
`deps/anymatrix-unitary-main` using the **exact upstream SHA that was
tested**. It never silently advances the submodule on `main`.

To enable this PR creation, open
[Actions settings](https://github.com/jaurentz/unitary-spectral-density/settings/actions)
and enable **Allow GitHub Actions to create and approve pull requests**.
The PR job requests `contents: write` and `pull-requests: write`;
the two MATLAB test jobs retain only read access. If GitHub blocks
automatic PR creation, the job emits a warning but does not invalidate
the successful numerical tests.

**GitHub limitation:** pushes and PRs created using the default
`GITHUB_TOKEN` ordinarily do not trigger new push/PR workflows.
The PR includes the originating successful two-mode MATLAB CI run.
For a new CI run on the PR itself, update its branch through a user or
dedicated GitHub App token, or manually rerun the workflow as needed.

## Local development

Work in two independent clones:

```bash
git clone https://github.com/jaurentz/anymatrix-unitary.git
git clone --recurse-submodules https://github.com/jaurentz/unitary-spectral-density.git
```

Change and test generators in the standalone library first. Then in
the research repository, test the latest library manually:

```bash
git submodule update --remote external/anymatrix-unitary
# Run MATLAB tests: setup; results=runtests('tests'); assertSuccess(results);
```

After validation, pin the chosen upstream commit in a reviewed commit/PR:

```bash
git add external/anymatrix-unitary
git commit -m "Update tested anymatrix-unitary"
git push
```

Do not modify the submodule working copy as the main library development
workspace; it is typically at a detached HEAD.
