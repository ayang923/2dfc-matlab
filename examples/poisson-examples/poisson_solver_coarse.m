function [u_num_mat, R, R_eval] = poisson_solver_coarse(curve_seq, f, u_boundary, h, G_cf, p, int_eps, eps_xi_eta, eps_xy, d, C, n_r, A, Q, M, h_eval, n_x_padded, n_y_padded, perturb)
% POISSON_SOLVER_COARSE  Solves Delta(u) = f using a coarser evaluation grid.
%
% Identical algorithm to poisson_solver, but the solution is evaluated on a
% separate coarser Cartesian mesh R_eval (step size h_eval > h) rather than
% on the fine FC grid R. This saves memory and time when visualizing or
% comparing against an exact solution at moderate resolution.
%
% Inputs:
%   curve_seq  - Curve_seq_obj describing the domain boundary
%   f          - Function handle for the Poisson right-hand side f(x,y)
%   u_boundary - Function handle for the Dirichlet boundary data u(x,y)
%   h          - Fine Cartesian mesh step size for the FC computation
%   G_cf       - Quadrature refinement factor for the IE boundary discretization
%   p          - Polynomial degree for the IE graded-mesh quadrature
%   int_eps    - Convergence tolerance for adaptive IE integration
%   eps_xi_eta - Newton inversion tolerance in parameter (xi,eta) space
%   eps_xy     - Newton inversion tolerance in physical (x,y) space
%   d          - Number of Gram matching points for 1D FC
%   C          - Number of continuation points for FC
%   n_r        - Refinement factor for the FC continuation grid
%   A          - Precomputed FC continuation matrix (n_r*C x d)
%   Q          - Precomputed FC Gram polynomial matrix (d x d)
%   M          - Polynomial interpolation degree
%   h_eval     - Step size of the coarse evaluation grid (h_eval >= h)
%   n_x_padded - (optional) Override for R grid size in x; pass [] to skip
%   n_y_padded - (optional) Override for R grid size in y; pass [] to skip
%   perturb    - (optional) If true, expand bounding box by rand(1)*h per side;
%                if false or omitted, expand by the fixed h
%
% Outputs:
%   u_num_mat - (n_y_eval x n_x_eval) solution values on R_eval; NaN outside domain
%   R         - Fine R_cartesian_mesh_obj (holds the particular solution in R.f_R)
%   R_eval    - Coarse R_cartesian_mesh_obj on which the solution is returned

    % Build the FC grid, forwarding any optional arguments to FC2D
    fc2d_args = {};
    if nargin >= 18
        fc2d_args{end+1} = n_x_padded;
        fc2d_args{end+1} = n_y_padded;
    end
    if nargin >= 19;  fc2d_args{end+1} = perturb;  end

    [R, ~, ~, ~] = FC2D(f, h, curve_seq, eps_xi_eta, eps_xy, d, C, n_r, A, Q, C, A, Q, M, fc2d_args{:});

    tic;

    R_eval = R_cartesian_mesh_obj(R.x_start, R.x_end - h_eval, R.y_start, R.y_end - h_eval, ...
        h_eval, R.boundary_X, R.boundary_Y);

    % Particular solution on the fine grid; interpolate onto the coarse grid
    R.f_R = R.inv_lap();
    R_eval.f_R = R_locally_compute_vec(R, R_eval.R_X, R_eval.R_Y, M);

    % Homogeneous BVP with modified boundary data u - u_p|boundary
    u_m_up_boundary = @(x, y) u_boundary(x, y) - R_locally_compute_vec(R, x, y, M);
    uh = laplace_solver_coarse(R, R_eval, curve_seq, u_m_up_boundary, G_cf, p, M, int_eps, eps_xi_eta, eps_xy, n_r);

    u_num_mat = uh + R_eval.f_R;
    toc;
end

% -------------------------------------------------------------------------

function u_num_mat = R_locally_compute_vec(R, X, Y, M)
% R_LOCALLY_COMPUTE_VEC  Evaluates R.locally_compute element-wise over arrays X, Y.
    u_num_mat = zeros(size(X));
    for i = 1:numel(X)
        u_num_mat(i) = R.locally_compute(X(i), Y(i), M);
    end
end
