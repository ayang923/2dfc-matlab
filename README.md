# 2DFC-MATLAB User Guide

> A comprehensive reference for understanding, setting up, and running the 2D Fourier Continuation code in this repository.

---

## Table of Contents

- [Overview](#overview)
- [Setup](#setup)
- [Repository Structure](#repository-structure)
- [How to Run the 2DFC Algorithm](#how-to-run-the-2dfc-algorithm)
  1. [Select Parameters and Load FC Matrices](#1-select-parameters-and-load-fc-matrices)
  2. [Define the Domain via `Curve_seq_obj`](#2-define-the-domain-via-curve_seq_obj)
  3. [Call `FC2D`](#3-call-fc2d)
- [Understanding the Output](#understanding-the-output)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)

---

## Overview

This repository implements the **2D Fourier Continuation (2DFC)** algorithm, which takes a function $f(x,y)$ defined on an arbitrary 2D domain bounded by a sequence of $C^2$ curves and produces a smooth, periodic Fourier series representation on a bounding Cartesian rectangle. This allows spectral-accuracy computations (derivatives, Poisson solvers, etc.) on non-rectangular domains.

The repository also includes a **Poisson/Laplace solver** built on top of 2DFC: the algorithm computes a spectral particular solution via the inverse Laplacian of the 2DFC extension, then solves the remaining homogeneous Laplace problem using a boundary integral equation (BIE) method.

The algorithm is described in detail in:

> Bruno, Oscar P., and Allen Yang. Two-Dimensional Fourier Continuation for Domains with Corners. Submitted manuscript, 2026.

---

## Setup

### Requirements

- MATLAB R2021a or later (tested on MATLAB R2024)
- Symbolic Math Toolbox (required only for precomputing FC matrices via `precomp_fc_data`)

### Installation

No installation is required. Add `src/` and `data/` to your MATLAB path and run any of the example scripts under `examples/`.

Precomputed FC matrices for are saved in `data/FC_data/`. If you need matrices for different parameters, use `generate_bdry_continuations` (see [Section 1](#1-select-parameters-and-load-fc-matrices)).

---

## Repository Structure

```
2dfc-matlab/
├── README.md                          This guide
├── data/
│   ├── generate_bdry_continuations.m  Generates and saves FC matrices
│   └── FC_data/                   Precomputed A and Q matrices (.mat and .txt)
│       ├── A_d{d}_C{C}_r{n_r}.mat
│       └── Q_d{d}_C{C}_r{n_r}.mat
├── src/                           Core library
│   ├── FC2D.m                     Main entry point
│   ├── FC2D_patches.m             Variant for pre-built patch arrays
│   ├── Curve_seq_obj.m            Circular linked list of boundary curves
│   ├── Curve_obj.m                Single C^2 curve with patch construction
│   ├── Q_patch_obj.m              Base class for all patch types
│   ├── S_patch_obj.m              Smooth (interior) boundary patch
│   ├── C2_patch_obj.m             Convex corner patch
│   ├── C1_patch_obj.m             Concave corner patch
│   ├── R_cartesian_mesh_obj.m     Cartesian mesh with spectral operators
│   ├── R_xi_eta_inversion.m       Batch Newton inversion of patch maps
│   ├── fcont_gram_blend_S.m       1D blending-to-zero FC (S-patches)
│   ├── fcont_gram_blend_C2.m      2D blending-to-zero FC (C2 corner)
│   ├── precomp_fc_data.m          Precomputes FC matrices A and Q
│   ├── barylag.m                  Barycentric Lagrange interpolation
│   ├── inpolygon_mesh.m           Fast in-polygon test for Cartesian grids
│   ├── newton_solve.m             General Newton root-finding
│   ├── mgs.m                      Modified Gram-Schmidt QR decomposition
│   └── shift_idx_mesh.m           Stencil shift utility for interpolation
└── examples/
    ├── 2DFC-examples/
    │   ├── boomerang_2D_FC.m      Single smooth closed curve
    │   ├── teardrop_2DFC.m        Teardrop with moderate tip
    │   ├── teardrop_sharp_2DFC.m  Very sharp teardrop (near-cusp)
    │   ├── heart_sharp_2DFC.m     Heart with near-cusp indentation
    │   └── guitarbase_2DFC.m      Four-curve guitar-body domain
    └── poisson-examples/
        ├── poisson_solver.m       Poisson solver (full resolution output)
        ├── poisson_solver_coarse.m Poisson solver (coarse evaluation grid)
        ├── laplace_solver.m       Laplace BVP via boundary integral equation
        └── ...
```

---

## How to Run the 2DFC Algorithm

Running the 2DFC algorithm requires three steps:

1. Select algorithm parameters and load (or precompute) the FC matrices
2. Define the domain boundary as a sequence of $C^2$ curves using `Curve_seq_obj`
3. Call `FC2D`

### 1. Select Parameters and Load FC Matrices

#### Parameter Reference

**Function and mesh parameters:**

| Parameter | Description |
|-----------|-------------|
| `f` | Function handle `f(x,y)` to be represented |
| `h` | Cartesian mesh step size |

**1D-FC parameters:**

| Parameter | Description |
|-----------|-------------|
| `d` | Number of Gram matching points for the 1D continuation; controls accuracy. Typical values: 6–10 |
| `C_S` | Number of continuation points for smooth (S-type) patches |
| `C_C` | Number of continuation points for corner (C1/C2-type) patches |
| `n_r` | Refinement factor: the continuation grid is `n_r` times finer than the patch mesh |

**Interpolation and Newton solver parameters:**

| Parameter | Description |
|-----------|-------------|
| `M` | Polynomial degree for interpolating patch values onto the Cartesian mesh; `d+3` is a good default |
| `eps_xi_eta` | Newton solver tolerance in parameter (ξ,η) space; `1e-13` to `1e-14` is typical |
| `eps_xy` | Newton solver tolerance in physical (x,y) space; set equal to `eps_xi_eta` in most cases |

**Optional arguments to `FC2D` (bounding box and grid size):**

After the required arguments through `M`, you may pass up to three optional trailing arguments. They control how the bounding rectangle for `R_cartesian_mesh_obj` is built from the patch extents.

| Parameter | Description |
|-----------|-------------|
| `n_x_padded` | With `n_y_padded`, overrides the default `R.n_x` / `x_end` from the bounding box (uniform spacing `h`). Must be **≥** the natural grid size for that box; see `R_cartesian_mesh_obj`. Padding is used only if **`nargin ≥ 16`** and **`n_x_padded` is non-empty**. |
| `n_y_padded` | Paired with `n_x_padded`: sets `R.n_y` and `y_end`. Supply **both** when you want a larger FFT grid than the minimal bounding rectangle. |
| `perturb` | Read only if **`nargin ≥ 17`**. **Logical:** how far to **inflate** the patch bounding box before placing the grid (reduces ambiguous inside/outside tests when a grid line coincides with the boundary polygon). **`true`**: expand each side by **`rand(1)*h`** (independent draws in x and y). **`false`** or **omit** (e.g. call with only 14 or 16 arguments): expand each side by the fixed amount **`h`** — **reproducible**. To set `perturb` without padding, call `FC2D(..., M, [], [], perturb)`. |

See [§3 Call `FC2D`](#3-call-fc2d) for example call patterns. The Poisson drivers `poisson_solver` / `poisson_solver_coarse` accept the same three optional arguments and forward them to `FC2D`.

#### Loading FC Matrices

The FC matrices `A` and `Q` encode the 1D continuation operation and must be precomputed for the chosen `(d, C, n_r)` triple. Precomputed matrices for common parameters are in `data/FC_data/`. Load them as follows:

```matlab
d = 8; C = 27; n_r = 6;

% Generate matrices if not already on disk
if ~exist(['FC_data/A_d',num2str(d),'_C',num2str(C),'_r',num2str(n_r),'.mat'])
    generate_bdry_continuations(d, C, C, 12, 20, 4, 256, n_r);
end

load(['FC_data/A_d',num2str(d),'_C',num2str(C),'_r',num2str(n_r),'.mat']);
load(['FC_data/Q_d',num2str(d),'_C',num2str(C),'_r',num2str(n_r),'.mat']);
A = double(A);
Q = double(Q);
```

The same matrices are used for both S-type patches (`A_S`, `Q_S`) and corner patches (`A_C`, `Q_C`). You can use different `C_S` and `C_C` values, but the precomputed matrices must match the respective `C` value.

---

### 2. Define the Domain via `Curve_seq_obj`

#### Background

The domain boundary must be expressible as an ordered, counter-clockwise sequence of $k$ $C^2$ curves $(c_1, c_2, \dots, c_k)$ such that the end of each curve meets the start of the next:

$$(\ell_1^i(1),\, \ell_2^i(1)) = (\ell_1^{i+1}(0),\, \ell_2^{i+1}(0)) \quad \text{for } i = 1, \dots, k$$

Each curve is parametrized as $(x, y) = (\ell_1(\theta),\, \ell_2(\theta))$ for $\theta \in [0,1]$.

#### Patches and Why They Are Needed

The 2DFC algorithm requires the function to be extended smoothly to zero near the boundary. It does this by tiling the boundary region with overlapping **patches** — small curvilinear coordinate systems in which a smooth blending-to-zero extension can be applied.

Two types of patches are constructed per curve:

- **S-patch (smooth patch):** covers most of the curve's interior region, away from the junctions at $\theta=0$ and $\theta=1$.
- **Corner patch:** centered at each junction point ($\theta=0$ of the current curve / $\theta=1$ of the previous curve, and vice versa). Corner patches handle the angular geometry at curve junctions, which the S-patch cannot.

The S-patch and corner patches overlap, and a **partition of unity** blends their contributions smoothly in the overlapping region.

> **Remark (convex C2 patches).** The **C2-type** (convex corner) patch construction and its partition-of-unity weights assume, for implementation simplicity, that **neighboring patches do not intersect** the region $\mathcal{M}_p^{\mathcal{C}_2}\bigl((d{-}1)h_\xi^{\mathcal{C}_2} \times (d{-}1)h_\eta^{\mathcal{C}_2}\bigr)$—the image in $(x,y)$ of the product $(d{-}1)h_\xi^{\mathcal{C}_2} \times (d{-}1)h_\eta^{\mathcal{C}_2}$ in $(\xi,\eta)$ under the parametrization $\mathcal{M}_p^{\mathcal{C}_2}$ (notations as in the paper). In practice, tune `frac_n_*` and `h_norm` so adjacent S/C1 patches avoid overlapping that set.

#### Patch Construction Parameters

To simplify patch construction, `add_curve` accepts fractional parameters that express patch sizes as fractions of the curve's total discretization count `n`, making them scale-invariant. Each curve is first **uniformly discretized** in parameter space: $\theta \in [0,1]$ is sampled at $n$ equally spaced points (or $n$ is chosen automatically from arc length and `h_norm`). All patch widths are then expressed as **fractions of $n$**, so the same relative layout works when you refine the mesh.

- A **fraction of the $n$ points** (given by `frac_n_C_0`) is allocated to the **corner patch centered at $\theta = 0$** — that patch uses the discretization near the start of the curve.
- Another **fraction** (`frac_n_C_1`) is allocated to the **corner patch centered at $\theta = 1$** — the discretization near the end of the curve.
- The remaining interior segment supports the **S-patch**; the `frac_n_S_0` and `frac_n_S_1` parameters then describe how the S-patch **overlaps** into each corner patch (for the partition of unity), again as fractions relative to the corner-patch segment sizes.

This “count points along $\theta$, then carve out corner regions by fractions” design keeps the interface simple and **scale-invariant**: you tune relative patch sizes without re-specifying physical lengths whenever $h$ changes.

Each call to `add_curve` accepts:

| Parameter | Description |
|-----------|-------------|
| `n` | **Uniform $\theta$-grid size.** Number of equally spaced samples of $\theta \in [0,1]$ (endpoints included), i.e. $n-1$ intervals in $\theta$. Pass **`0`** to auto-set `n = ceil(arc_length / h_norm) + 1`. |
| `frac_n_C_0` | **Corner patch at $\theta=0$.** Fraction of **`n`** used for the corner patch at the **start** of this segment (junction with the previous curve). Counts: `n_C_0 = ceil(frac_n_C_0 * n)` $\theta$-samples from the low-$\theta$ side. Pass **`0`** for default `ceil(n/10)`. |
| `frac_n_C_1` | **Corner patch at $\theta=1$.** Fraction of **`n`** for the corner patch at the **end** of this segment (junction with the next curve). Counts: `n_C_1 = ceil(frac_n_C_1 * n)` samples from the high-$\theta$ side. Pass **`0`** for default `ceil(n/10)`. |
| `frac_n_S_0` | **S-patch overlap at the $\theta=0$ corner.** Fraction of **`n_C_0`** (not of `n`) by which the smooth patch extends into the $\theta=0$ corner patch: `n_S_0 = ceil(frac_n_S_0 * n_C_0)` — this sets the partition-of-unity overlap width. Pass **`0`** for default `ceil(2/3 * n_C_0)`. |
| `frac_n_S_1` | **S-patch overlap at the $\theta=1$ corner.** Fraction of **`n_C_1`**: `n_S_1 = ceil(frac_n_S_1 * n_C_1)`. Pass **`0`** for default `ceil(2/3 * n_C_1)`. |
| `h_norm` | **Normal-direction resolution.** Physical step size in $(x,y)$ for the patch mesh **inward** from the boundary; with `n` it determines how the boundary strip is sampled perpendicular to the curve. |

Larger `frac_n_C_*` gives wider corner patches along $\theta$; larger `frac_n_S_*` widens the S-patch’s overlap into each corner patch (stronger S contribution in the blend).

> **More control:** The `add_curve` interface is intentionally simple. If you need finer control — e.g., setting exact arc-length widths, non-uniform discretizations, or directly constructing `S_patch_obj` and `C1/C2_patch_obj` objects — you can build patches manually as well. This is more involved but gives complete flexibility over the patch geometry.

#### Example: Single Smooth Curve

```matlab
curve_seq = Curve_seq_obj();
curve_seq.add_curve( ...
    l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, ...
    0,     ...  % n: auto-compute from arc length
    0.1,   ...  % frac_n_C_0: corner patch at theta=0 spans 10% of n
    0.1,   ...  % frac_n_C_1: corner patch at theta=1 spans 10% of n
    0.6,   ...  % frac_n_S_0: S-patch overlaps 60% into the theta=0 corner patch
    0.6,   ...  % frac_n_S_1: S-patch overlaps 60% into the theta=1 corner patch
    h_norm);
```

Call `curve_seq.plot_geometry(d)` after adding all curves to visualize the patch boundaries before running the full algorithm.

---

### 3. Call `FC2D`

**Minimal call** (fixed `h` expansion, auto grid size from bounding box):

```matlab
[R, interior_patches, FC_patches, fc_err] = FC2D( ...
    f, h, curve_seq, eps_xi_eta, eps_xy, ...
    d, C_S, n_r, A_S, Q_S, C_C, A_C, Q_C, M);
```

**Optional trailing arguments** (in order: `n_x_padded`, `n_y_padded`, `perturb`):

| Goal | Example |
|------|---------|
| Override grid dimensions only | `FC2D(..., M, nx, ny)` — larger `nx`,`ny` add **padding** around the domain (same `h`; see `R_cartesian_mesh_obj`) |
| Random bounding-box expansion | `FC2D(..., M, [], [], true)` — `rand(1)*h` per side |
| Padding + random expansion | `FC2D(..., M, nx, ny, true)` |
| Padding + fixed `h` expansion | `FC2D(..., M, nx, ny, false)` or `FC2D(..., M, nx, ny)` — last arg omitted defaults to fixed `h` |

Padding is enabled only when **`n_x_padded` is non-empty** (`nargin >= 16` and `~isempty(n_x_padded)`). To pass **`perturb`** without padding, use **placeholders**: `FC2D(..., M, [], [], perturb)`.

`FC2D` prints absolute max, relative max, and relative $L^2$ errors to the command window as a quick accuracy check.

---

## Understanding the Output

| Output | Type | Description |
|--------|------|-------------|
| `R` | `R_cartesian_mesh_obj` | The Cartesian mesh; contains `f_R` (function values), `fc_coeffs` (Fourier coefficients), and spectral operator methods |
| `interior_patches` | cell array | `{S_patch_1, C_patch_1, S_patch_2, C_patch_2, ...}` for each curve |
| `FC_patches` | cell array | Four `Q_patch_obj`s per curve: one S extension and three corner extensions |
| `fc_err` | double | Relative $L^2$ error of the Fourier approximation (evaluated on a $2\times$ finer grid) |

The most commonly used output is `R`. Key fields and methods:

- `R.fc_coeffs` — `(n_y × n_x)` array of Fourier coefficients in `fftshift` order (DC at center)
- `R.f_R` — `(n_y × n_x)` assembled continuation function values on the Cartesian grid
- `R.in_interior` — `(n_y × n_x)` logical mask of interior points
- `R.grad(f_mesh)` — spectral gradient
- `R.lap(f_mesh)` — spectral Laplacian
- `R.inv_lap()` — spectral inverse Laplacian (particular solution for Poisson)

---

## Examples

The code used to generate the 2D-FC examples considered in the paper are in `examples/2DFC-examples/` and the code used for the Poisson examples considered in the paper are in `examples/poisson-examples/`

---

## Troubleshooting

**Error is large or doesn't converge:**
- Call `curve_seq.plot_geometry(d)` to visually inspect patch placement; adjust `frac_n_C_0/1` and `frac_n_S_0/1` if patches overlap incorrectly
- Try increasing `h_norm` or lowering `n_curve` to make the blending-to-zero region larger
- Increasing the patch overlap regions

**Newton solver non-convergence warnings:**
- These usually arise at extreme patch boundaries; a small number are tolerable
- This can happen if patch decomposition is inconsistent or polygonal approximation of global boundary is not refined enough
- Loosen `eps_xi_eta` / `eps_xy` if the issue is widespread

**Need FC matrices for new parameters:**
- Run `generate_bdry_continuations(d, C, C, 12, 20, 4, 256, n_r)` (requires Symbolic Math Toolbox; takes several minutes for large `d`)
