function [u_num_mat, R] = poisson_solver(curve_seq, f, u_boundary, h, G_cf, p, int_eps, eps_xi_eta, eps_xy, d, C, n_r, A, Q, M)
%POISSON_SOLVER Summary of this function goes here
%   Detailed explanation goes here
[R, ~, ~, fc_err] = FC2D(f, h, curve_seq, eps_xi_eta, eps_xy, d, C, n_r, A, Q, C, A, Q, M);

% compute particular solution
R.f_R = R.inv_lap;

u_m_up_boundary = @(x, y) u_boundary(x, y) - R_locally_compute_vec(R, x, y, M);
uh = laplace_solver(R, curve_seq, u_m_up_boundary, G_cf, p, M, int_eps, eps_xi_eta, eps_xy, n_r);

u_num_mat = uh+R.f_R;
end

function u_num_mat = R_locally_compute_vec(R, X, Y, M)
    u_num_mat = zeros(size(X));
    for i = 1:numel(X)
        u_num_mat(i) = R.locally_compute(X(i), Y(i), M);
    end
end