# Mathematics and literature: unitary spectral-density estimators

Experimental MATLAB implementations, November 2026. This document distinguishes
established unit-circle theory, Hermitian analogues, and unproved combinations.

## 1. Spectral measure and normalization

For a unitary n-by-n matrix U with eigenvalues exp(i phi_j), the finite
empirical phase measure is mu_U = (1/n) sum_j delta_{phi_j}. Its Fourier
moments are m_k = tr(U^k)/n. All densities here are measured against
dtheta on the circle (integral over one period equals 1). A finite
measure is atomic, even when a smooth limiting density exists.

The full Poisson regularization at 0<r<1 is

    P_r(t) = (1-r^2)/(1+r^2-2*r*cos(t))
    rho_r(theta) = (1/(2*pi*n)) sum_j P_r(theta-phi_j).

For Fourier multipliers g_k and a finite degree K, the estimator targets

    rho_g(theta) = [1+2*Re(sum_{k=1}^K g_k*m_k*exp(-i*k*theta))]/(2*pi).

Poisson multipliers g_k=r^k have an INFINITE series. Truncating this
series may cause negative values. Distinguish this empirical trace measure
from the spectral measure of an OPUC distinguished vector and from the
infinite-size or ensemble-averaged spectral distribution. See the
[companion matrix library mathematics](https://github.com/jaurentz/anymatrix-unitary/blob/main/docs/MATHEMATICS.md)
and [spectral reference interface](https://github.com/jaurentz/anymatrix-unitary/blob/main/docs/SPECTRAL_REFERENCE.md).

## 2. Trace moments and smoothing on the circle

If v is a random unit vector with E[v*v']=I/n, then

    E[v'*U^k*v] = tr(U^k)/n.

The implementation supports independent four-phase, real Rademacher,
or normalized complex-Gaussian probes. The same probe is reused for
powers k and for every evaluation angle. A complete canonical basis
instead gives exact moments at cost n*K matvecs. Arbitrary custom
vectors are not necessarily isotropic and need not yield an unbiased
trace estimator. Pointwise sample standard errors are not uniform
high-probability error guarantees.

The smoothing methods are:

| Name | Fourier multiplier g_k | Notes |
|---|---|---|
| poisson | r^k | Full kernel positive; finite truncation need not be |
| fejer | 1-k/(K+1) | Positive trigonometric kernel |
| jackson | Normalized squared-Fejer Fourier coefficient | Positive trigonometric kernel |
| none / dirichlet | 1 | Gibbs oscillations; no positivity |
| lanczos | sin(pi*k/(K+1))/(pi*k/(K+1)) | Positivity not guaranteed |
| heat / wrapped_gaussian | exp(-sigma^2*k^2/2) | Full wrapped heat kernel positive |
| poisson-fejer | r^k times Fejer multiplier | Positive finite kernel |
| poisson-jackson | r^k times Jackson multiplier | Positive finite kernel |

For the **circle** Jackson kernel, let

    F_q(t) = (1/q) * abs(sum_{j=0}^{q-1} exp(i*j*t))^2
    J_q(t) = F_q(t)^2 / [(1/2pi) integral F_q(s)^2 ds].

Its coefficients are the self-convolution of the triangular sequence
a_j=(q-|j|)_+, normalized by the central convolution coefficient.
The degree is 2(q-1). Given requested K, choose q=floor(K/2)+1,
which yields degree at most K. Zero coefficients at the high end
(when K is odd) are intentional. This is not the standard Jackson
coefficient formula for real-line Chebyshev KPM. Applying the latter
directly to the circle would require an independent justification.

Plain Fejer, Jackson, Dirichlet, Lanczos and heat windows ignore the
positional radius r. Explicit Poisson-* compositions multiply by r^k.
A Jackson-smoothed reference and a Poisson-smoothed reference are
**different targets** and should not be conflated in accuracy plots.

## 3. Isometric Arnoldi and Gauss-Szego quadrature

An m-step isometric Arnoldi run on a unitary U generates a unitary
Hessenberg completion H_m. Let its eigenvalues be exp(i theta_l) and
normalized eigenvectors be columns of W. Define

    w_l=abs(W(1,l))^2,    sum_l w_l=1.

The weighted m-node circle quadrature satisfies, in exact arithmetic,

    v'*U^k*v = sum_l w_l*exp(i*k*theta_l),    |k| <= m-1,

where v is the normalized Arnoldi starting vector. Using equal node
weights is incorrect. Random isotropic starts and averaging yield a
trace-density estimator. For a trigonometric smoothing polynomial of
degree K<=m-1, the quadrature matches the corresponding Fourier-moment
estimate for the *same* starts, up to roundoff. An infinite Poisson
kernel incurs additional quadrature truncation error.

Current code uses twice-reorthogonalized Arnoldi and dense unitary
closure, not a specialized short-recurrence isometric Arnoldi method.
The Hermitian counterpart is Lanczos and Gaussian quadrature / SLQ.

## 4. Poisson resolvents and sketched Krylov projections

For z=r*exp(i*theta) in the open unit disk,

    rho_r(theta) = Re tr[(U+zI)(U-zI)^(-1)]/(2*pi*n)
                 = Re[1+(2*z/n)*tr((U-zI)^(-1))]/(2*pi).

With solver='auto' (the default), small operators of dimension at most
128 are assembled explicitly by applying the operator to the n canonical
basis vectors and solved using MATLAB backslash. The assembled matrix
costs n matvecs and O(n^2) storage. Each shift is refactorized, but
all stochastic probe right-hand sides are handled together; the resulting
solve error should be at the usual backward-stable direct-solver level.
Above this cutoff, the default is MATLAB GMRES, retaining matrix-free
storage. The cutoff is configurable via direct_threshold. Either method
can be selected explicitly via solver='backslash' or solver='gmres'.
Errors include probe variance and, for GMRES, iterative solve error.

The experimental sketch method constructs once per probe a Krylov basis
Q_m and AQ_m=U*Q_m, then for each z solves the projected random
least-squares problem

    c_z = argmin_c ||S*((AQ_m-z*Q_m)*c-v)||_2,
    y_z=Q_m*c_z,

where S is a Gaussian or sparse CountSketch embedding. The same Q_m,
AQ_m and S are reused at all phase angles, with sketch randomness
independent of trace-probe randomness. The code evaluates the **true**
residual ||(U-zI)*y_z-v||/||v|| from the cached AQ_m. An optional
fallback retries inaccurate projected solves with full GMRES.

This is a correctness-oriented research prototype, NOT an optimized
low-memory sGMRES implementation. Building the full Arnoldi basis still
costs O(n*m) storage and reorthogonalization; Gaussian sketches store
O(q*n) entries, CountSketch O(n) nonzeros. Small sketches may fail
to embed the Krylov space, and without fallback a projected density
may be inaccurate or negative. A low true residual controls solve
accuracy, but the resulting trace estimate retains probe error.

**Unitary-specific caveat:** for |z|<1, the convergent Neumann expansion
of (U-zI)^(-1) uses negative powers of U. A forward (positive-power)
Krylov space may approximate a resolvent poorly for spectrum wrapping
around the whole circle. Adjoint, rational or inverse Krylov spaces are
important future alternatives.

The Hermitian analogues are Gaussian/Lorentzian resolvent filters and
sketched GMRES/GMRES-like projection algorithms for linear systems.
General randomized-sketching theory does not automatically establish
uniform eigenphase-density error bounds.

## 5. Fair experimental comparison

Separate smoothing bias, Fourier truncation, randomized sampling,
quadrature approximation, linear-system residual, and finite-versus-limit
reference error. Report smoothing weights, number of matvecs, probes,
Arnoldi order, sketch rows and distribution, pointwise sample errors,
and *true* shifted-system residuals. A Poisson resolvent estimator
does not compute a Jackson-smoothed density.

For a matched Jackson experiment, select the same K in moments and
quadrature, with Arnoldi order m>=K+1. For a matched Poisson target,
use the same radius for moment, quadrature and resolvent routines,
accounting for moment and quadrature truncation.

## 6. Verified source literature

### Unit circle / unitary methods

- W. B. Gragg, *Positive Definite Toeplitz Matrices, the Arnoldi Process
  for Isometric Operators, and Gaussian Quadrature on the Unit Circle*,
  J. Comput. Appl. Math. 46 (1993), 183-198.
  https://doi.org/10.1016/0377-0427(93)90294-L
- A. Bultheel, M. J. Cantero and R. Cruz-Barroso, *Matrix Methods for
  Quadrature Formulas on the Unit Circle: A Survey*, J. Comput. Appl.
  Math. 284 (2015), 78-100.
  https://doi.org/10.1016/j.cam.2014.11.002
- B. Simon, *Orthogonal Polynomials on the Unit Circle*, AMS Colloquium
  Publications, vol. 54, Parts 1-2 (2005).
  https://doi.org/10.1090/coll054.1
- B. Simon, *OPUC on One Foot*, Bull. Amer. Math. Soc. 42 (2005),
  431-460. https://doi.org/10.1090/S0273-0979-05-01075-X
- Classical positive trigonometric approximation: Fejer kernels and
  squared-Fejer Jackson kernels, e.g. *Approximation of Univariate
  Functions*, section 1.2.4 (Jackson kernel).
- N. J. Higham and M. Mikaitis, *Anymatrix: An Extensible MATLAB
  Matrix Collection*, Numer. Algorithms 90 (2022), 1175-1196.
  https://doi.org/10.1007/s11075-021-01226-2

### Hermitian density, moments, stochastic traces and quadrature

- L. Lin, Y. Saad and C. Yang, *Approximating Spectral Densities of Large
  Matrices*, SIAM Rev. 58 (2016), 34-65.
  https://doi.org/10.1137/130934283
- A. Weisse, G. Wellein, A. Alvermann and H. Fehske, *The Kernel
  Polynomial Method*, Rev. Mod. Phys. 78 (2006), 275-306.
  https://doi.org/10.1103/RevModPhys.78.275
- G. H. Golub and G. Meurant, *Matrices, Moments and Quadrature with
  Applications*, Princeton University Press (2010).
  https://press.princeton.edu/books/ebook/9781400833887/matrices-moments-and-quadrature-with-applications-pdf
- S. Ubaru, J. Chen and Y. Saad, *Fast Estimation of tr(f(A)) via
  Stochastic Lanczos Quadrature*, SIAM J. Matrix Anal. Appl. 38
  (2017), 1075-1099.
- T. Chen, T. Trogdon and S. Ubaru, *Randomized Matrix-Free Quadrature:
  Unified and Uniform Bounds for Stochastic Lanczos Quadrature and
  the Kernel Polynomial Method* (2025).
  https://doi.org/10.1137/23M1600414
- M. F. Hutchinson, *A Stochastic Estimator of the Trace of the
  Influence Matrix for Laplacian Smoothing Splines* (1989).
  https://doi.org/10.1080/03610918908812806
- H. Avron and S. Toledo, *Randomized Algorithms for Estimating the
  Trace of an Implicit Symmetric Positive Semi-definite Matrix*,
  J. ACM 58 (2011), Article 8.
  https://doi.org/10.1145/1944345.1944349
  (PSD-specific error bounds do not directly apply to unitary U.)

### Projected systems and spectral filters

- Y. Nakatsukasa and J. A. Tropp, *Fast and Accurate Randomized
  Algorithms for Linear Systems and Eigenvalue Problems* (2024).
  https://doi.org/10.1137/23M1565413
- S. Guettel and M. Schweitzer, *Randomized Sketching for Krylov
  Approximations of Large-Scale Matrix Functions* (2023).
  https://doi.org/10.1137/22M1518062
- L. Lin, *Randomized Estimation of Spectral Densities of Large Matrices
  Made Accurate*, Numerische Mathematik (2017).
  https://doi.org/10.1007/s00211-016-0837-7

### Status and novelty

Established unitary mathematics includes OPUC/CMV, circle moments,
isometric Arnoldi, Szego quadrature and the Poisson/Caratheodory
transform. Stochastic trace identities, GMRES, and randomized
least-squares sketches are also established in their respective
general settings. Combining circle Jackson damping, randomized
Szego rules for empirical trace measures, or sketched Krylov solves
inside a unitary resolvent trace estimator is implemented here as an
experimental adaptation. We make **no claim that a combination is new**
without a dedicated literature search and proof, and we do not
transfer Hermitian-specific bounds without justification.
