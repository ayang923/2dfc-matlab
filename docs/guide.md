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

> Bruno, O. P., & Lyon, M. (2010). High-order unconditionally stable FC-AD solvers for general smooth domains I. Basic elements. *Journal of Computational Physics*, 229(6), 2009–2033.

---

## Setup

### Requirements

- MATLAB R2021a or later (tested on MATLAB R2024)
- Symbolic Math Toolbox (required only for precomputing FC matrices via `precomp_fc_data`)

### Installation

No installation is required. Add `src/` to your MATLAB path and run any of the example scripts under `examples/`.

Precomputed FC matrices for common parameter values are included in `data/FC_data/`. If you need matrices for different parameters, use `generate_bdry_continuations` (see [Section 1](#1-select-parameters-and-load-fc-matrices)).

---

## Repository Structure

```
2dfc-matlab/
├── docs/
│   └── guide.md                   This guide
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

The domain boundary must be expressible as an ordered, counter-clockwise sequence of $k$ $C^2$ curves $(c_1, c_2, \dots, c_k)$ such that:

$$(\ell_1^i(1),\, \ell_2^i(1)) = (\ell_1^{i+1}(0),\, \ell_2^{i+1}(0)) \quad \text{for } i = 1, \dots, k$$

$$(\ell_1^k(1),\, \ell_2^k(1)) = (\ell_1^1(0),\, \ell_2^1(0))$$

Each curve is parametrized as $(x, y) = (\ell_1(\theta),\, \ell_2(\theta))$ for $\theta \in [0,1]$. The algorithm constructs one **smooth (S-type) patch** covering the interior of each curve, and one **corner patch** at each junction $\theta = 1$ of the current curve / $\theta = 0$ of the next curve.

#### Patch Construction Parameters

Each call to `add_curve` requires the curve parametrization functions and the following parameters controlling how the patches are sized:

| Parameter | Description |
|-----------|-------------|
| `n` | Number of uniform discretization points along $\theta \in [0,1]$. If `0`, computed automatically as $\lceil L / h_\text{norm} \rceil + 1$ where $L$ is the arc length |
| `frac_n_C_0` | Fraction of `n` used to define the corner patch width at $\theta = 0$. Default (if `0`): $1/10$ |
| `frac_n_C_1` | Fraction of `n` used to define the corner patch width at $\theta = 1$. Default (if `0`): $1/10$ |
| `frac_n_S_0` | Fraction of the $\theta=0$ corner patch that overlaps the S-patch. Default (if `0`): $2/3$ |
| `frac_n_S_1` | Fraction of the $\theta=1$ corner patch that overlaps the S-patch. Default (if `0`): $2/3$ |
| `h_norm` | Normal step size (physical distance) for the S-patch mesh in the direction inward from the boundary |

Intuitively: `frac_n_C_0` and `frac_n_C_1` control the *size* of the two corner patches, while `frac_n_S_0` and `frac_n_S_1` control how far the S-patch *overlaps* into each corner patch (needed for the partition of unity).

#### Example: Single Smooth Curve

```matlab
curve_seq = Curve_seq_obj();
curve_seq.add_curve( ...
    l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, ...
    0,     ...  % n: auto-compute
    0.1,   ...  % frac_n_C_0
    0.1,   ...  % frac_n_C_1
    0.6,   ...  % frac_n_S_0
    0.6,   ...  % frac_n_S_1
    h_norm);
```

You can call `curve_seq.plot_geometry(d)` after adding all curves to visualize the patch boundaries before running the full algorithm.

---

### 3. Call `FC2D`

```matlab
[R, interior_patches, FC_patches, fc_err] = FC2D( ...
    f, h, curve_seq, eps_xi_eta, eps_xy, ...
    d, C_S, n_r, A_S, Q_S, C_C, A_C, Q_C, M);
```

`FC2D` will print absolute max, relative max, and relative $L^2$ errors to the command window as a quick accuracy check.

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

All examples are in `examples/2DFC-examples/`. Each follows the same three-step pattern: define `f`, build `curve_seq`, call `FC2D`.

| Script | Domain | Notes |
|--------|--------|-------|
| `boomerang_2D_FC.m` | Boomerang / figure-eight | Single smooth curve |
| `teardrop_2DFC.m` | Teardrop, moderate tip | Single smooth curve |
| `teardrop_sharp_2DFC.m` | Teardrop, near-cusp tip | Very fine mesh required |
| `heart_sharp_2DFC.m` | Heart with near-cusp | Near-cusp indentation |
| `guitarbase_2DFC.m` | Guitar-body shape | Four connected curves |

Poisson/Laplace solver examples are in `examples/poisson-examples/`.

---

## Troubleshooting

**Error is large or doesn't converge:**
- Refine the Cartesian mesh (reduce `h`)
- Call `curve_seq.plot_geometry(d)` to visually inspect patch placement; adjust `frac_n_C_0/1` and `frac_n_S_0/1` if patches overlap incorrectly
- Try increasing `d` or `M`

**Newton solver non-convergence warnings:**
- These usually arise at extreme patch boundaries; a small number are tolerable
- Tighten `eps_xi_eta` / `eps_xy` if the issue is widespread

**Out-of-memory errors:**
- Reduce `n_r`, or use a coarser mesh and larger `d` to maintain accuracy with fewer grid points

**Need FC matrices for new parameters:**
- Run `generate_bdry_continuations(d, C, C, 12, 20, 4, 256, n_r)` (requires Symbolic Math Toolbox; takes several minutes for large `d`)
