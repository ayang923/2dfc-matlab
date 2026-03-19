# How to Use This Code

> A comprehensive guide to understand, set up, and run the code in this repository.

---

## Table of Contents

- Overview
- Setup Instructions
  - Prerequisites
  - Installation
- Code Structure
- How to Run
  - Selecting parameters and loading pre-computed FC matrices
  - Defining domain and selecting discretization parameters using the Curve_seq_obj object
  - Calling the FC2D method
- Troubleshooting

---

## Overview

## Setup Instructions
All of this code was written and tested on MATLAB2021a but should work on other versions as well.

### Prerequisites
Code Tested on MATLAB 2024
Symbolic Math Toolbox

### Installation
N/A

### Code Structure
`FC_2D.m` contains "main" 2DFC function  
  
`Curve_seq_obj.m` contains class definition for "curve sequence" objects that stores the domain of the boundary. The algorithm assume that the domain it considers is the union of a finite number of $$C^2$$ curves.
`Curve_obj.m` contains class definition of "curve" objects that encodes each $$C^2$$ curve in the domain as well as parameters that deterimine how the domain is discretized.
  
`S_patch_obj.m`, `C1_patch_obj.m`, and `C2_patch_obj.m` contain class definitions for each patch type. S-type patches denote patches constructed along smooth portions of boundary (twice continuously differentiable), C1-type patches denote patches constructed along concave corners (interior angle > 180 degrees), and C2-type patches denote patches constructed along convex corners (interior angle < 180 degrees)  
`Q_patch_obj.m` contains class definition for square patches. The bulk of useable functions are defined in this class and each of the S-type, C1-type, and C2-type classes inherit from this class.  

`R_cartesian_mesh_obj.m` contains the class definition for the cartesian mesh that the FFT is computed on. The function values of each patch are interpolated on to this mesh to obtain the continuation function values on the Cartesian mesh.  

`boomerang_2D_FC.m`, `teardrop_2DFC.m`, and `polynomial_curve.m` contains examples of the 2D-FC algorithm for three different examples.

## How to Run
Running the 2DFC algorithm consists of three main steps:
1. Selecting parameters and loading pre-computed FC matrices
2. Defining domain and selecting discretization parameters using the Curve_seq_obj object
3. Calling the FC2D method

Below, I will elaborate in each step in more detal.

### 1. Selecting parameters and loading pre-computed FC matrices
#### Background
This step consists mostly of defining the `FC2D` method's function inputs. Here, I consider `f` a parameter even if it really is an input for the 2D-FC method in general.

#### Relevant Parameters
Domain and Function Parameters:
- `f` - Function handle for 2D function that is to be approximated with a Fourier expansion
- `curve_seq` - Curve_seq_obj describing domain of interest (more in second step)

1DFC Parameters:  
- `d` - 1D-FC degree
- `C_S` - Number of 1DFC continuation points used for smooth patches
- `C_C` - Number of 1DFC continuation points used for corner patches
- `n_r` - Refinement factor of 1DFC (for example, if n_r = 6, then the values of the continuation are given on a mesh that is 6 times finer than the mesh of the original function values)

Patch Parameters:
- `eps_xi_eta` - Precision each patch uses in parameter space (xi-eta space)
- `eps_xy` - Precision each patch uses in real space (x-y space); note this does not refer to the accuracy of the Fourier Series representation of the continuation function, it merely denotes what value is used for epsilon when using Newton's method to invert the parametrization of each patch.
- `M` - Degree of 1D polynomial interpolation used to compute function values in real space (on the basis of function values in parameter space of patch)

Mesh Parameters:
- `h` - Step size of Cartesian mesh used

#### Other Implementation Details
After selecting the 1DFC parameters described, either compute or load the continuation matrices `A_S`, `Q_S`, `A_C`, and `Q_C` as double precision matrices; the matrices `A_S` and `Q_S` are dependent on the parameter `C_S`, the matrices `A_C` and `Q_C` are dependent on the parameter `C_C`, and all four matrices are dependent on `d` and `n_r`. If the continuation matrices required aren't stored on disk already, the method `generate_bdry_continuations` can be used to compute these matrices and store them on disk. The parameter selection for this method will not be described here.

### 2. Defining domain and selecting discretization parameters using the Curve_seq_obj object
#### Background
This step consists mostly of defining the sequence of $$k$$ $$C^2$$ curves $$(c_1, c_2, \dots, c_k)$$ that define the domain of interest as well as selecting parameters that determine how each patch is constructed and discretized. Each curve $$c_i$$ is $$C^2$$ in the sense that it can be parameterized as the points $$(x, y) = (\ell_1^i(\theta), \ell_2^i(\theta))$$ for all $$\theta \in [0,1]$$ such that $$\ell_1$$ and $$\ell_2$$ are both twice differentiable functions. This algorithm assumes that the sequence $$c_1, c_2, \dots, c_k)$$ is in order, i.e. $$(\ell_1^1(1), \ell_2^1(1)) = (\ell_1^2(0), \ell_2^2(0)),\ (\ell_1^2(1), \ell_2^2(1)) = (\ell_1^3(0), \ell_2^3(0)),\ \dots,\ (\ell_1^k(1), \ell_2^k(1)) = (\ell_1^1(0), \ell_2^1(0))$$, and positively oriented with respect to $$\theta$$ (counter-clockwise). Then, for each curve, since the entire domain is the union of several of these curves, the corners of this domain would arise at $$\theta=0$$ and $$\theta=1$$. Thus, each curve has one smooth patch and two corner patches associated with it.
  
Further, this parametrization admits a natural discretization and method of patch construction via $$\theta$$. In detail, we can discretize the curve with respect to theta using a uniform mesh over the interval $$[0, 1]$$. Then, by choosing appropriate points from this discretization in addition to some other parameters, we can construct discretized meshes of each patch and their associated parametrizations.

#### Relevant Parameters
- `l_1` - Parametrization of x-coordinate of curve, i.e. $$\ell_1(\theta)$$
- `l_2` - Parametrization of y-coordinate of curve, i.e. $$\ell_2(\theta)$$
- `l_1_prime` - First derivative of `l_1`, i.e. $$\ell_1'(\theta)$$
- `l_2_prime` - First derivative of `l_2`, i.e. $$\ell_2'(\theta)$$
- `l_1_dprime` - Second derivative of `l_1`, i.e. $$\ell_1''(\theta)$$
- `l_2_dprime` - Second derivative of `l_2`, i.e. $$\ell_2''(\theta)$$
- `n` - Number of points discretization points for curve; if 0 is entered, then the length of the curve is calculated and `n` is set to length of the curve divided by `h_norm` (described below) + 1, the idea behind this "default value" is matching the step size in the tangential direction in real space to be around the same as the step size in the normal direction in real space
- `frac_n_C_0` - Fraction of total points that are used on this curve to construct the corner patch around $$(\ell_1(0), \ell_2(0))$$
- `frac_n_C_1` - Fraction of total points that are used on this curve to construct the corner patch around $$(\ell_1(1), \ell_2(1))$$
- `frac_n_S_0` - Fraction of points of the corner patch around $$(\ell_1(0), \ell_2(0))$$ that the smooth patch intersects
- `frac_n_S_1` - Fraction of points of the corner patch around $$(\ell_1(1), \ell_2(1))$$ that the smooth patch intersects
- `h_norm` - Step size in normal direction to boundary of discretization of smooth patch in real space

Essentially, `frac_n_C_0` and `frac_n_C_1` control the size of the two corner patches relative to the curve while `frac_n_S_0` and `frac_n_S_1` control how much the smooth patch intersects the two corner patches.

#### Other Implementation Details
The Curve_seq_obj object acts as a circular linked list, and each curve must be added individually (and in order) to this list using the `add_curve` method associated with this object. This method takes in the parameters described above. Making sure the domain is parametrized and inputted according to the assumptions described is essential for this method working. 

### 3. Calling the FC2D method
Once everything in step 1 and 2 are done, then running the 2DFC algorithm is as simple as running the `FC2D` function with the correct parameters.

### Expected Outputs
The `FC2D` function returns a `R_cartesian_mesh_obj` object and a cell array of the patch objects. In addition, the error of the Fourier expansion approximation should be displayed in the command window. The `R_cartesian_mesh_obj` object contains all relevant information, including the Cartesian mesh information, continuation function values on Cartesian mesh, and Fourier expansion coefficients. The cell array of patch objects is mainly used for debugging and confirming that the patches were constructed as expected. 

## Troubleshooting
If the error blows up or doesn't converge, consider:
- Refining mesh
- Plotting the patch meshes and playing around with the patch parameters associated with each curve
