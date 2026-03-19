function [u_num_mat] = laplace_solver(R, curve_seq, u_G, cfac, p, M, int_eps, eps_xi_eta, eps_xy, rho)
%LAPLACE_SOLVER_NEW Summary of this function goes here
%   Detailed explanation goes here

%% Computes the n discretization values
curve_n_rho1 = zeros(curve_seq.n_curves, 1);
curr = curve_seq.first_curve;
for i = 1:curve_seq.n_curves
    curve_n_rho1(i) = ceil((curr.n-1)*cfac);
    curr = curr.next_curve;
end
curve_param_rho1 = curve_param_obj(curve_n_rho1);

%% Constructs IE curve seq obj

IE_curve_seq = IE_curve_seq_obj(curve_seq, p);
[A_rho1, b_rho1] = IE_curve_seq.construct_A_b(curve_param_rho1, u_G);
gr_phi_rho1 = A_rho1 \ b_rho1;
gr_phi_fft_rho1 = fftshift(fft(gr_phi_rho1))/curve_param_rho1.n_total;

[s_patches, c_0_patches, c_1_patches] = IE_curve_seq.construct_interior_patches(curve_param_rho1, R.h, M, eps_xi_eta, eps_xy);
[well_interior_msk, s_patch_msks] = gen_R_msks(R, rho, s_patches, c_0_patches, c_1_patches);

%% Fills in points that are well within the domain
u_num_mat = zeros(size(R.f_R));
u_num_mat_fine = zeros(size(R.f_R));

disp('interior')
for idx = R.R_idxs(well_interior_msk)'
    u_num_mat(idx) = IE_curve_seq.u_num(R.R_X(idx), R.R_Y(idx), curve_param_rho1, gr_phi_rho1);
end

to_update = well_interior_msk;
rho_fine_IE = 2;

while sum(to_update, 'all') > 0
    [gr_phi_fine, curve_param_fine] = refine_gr_phi(curve_param_rho1, rho_fine_IE, gr_phi_fft_rho1);

    for idx = R.R_idxs(to_update)'
        u_num_mat_fine(idx) = IE_curve_seq.u_num(R.R_X(idx), R.R_Y(idx), curve_param_fine, gr_phi_fine);
    end

    to_update = to_update & abs(u_num_mat_fine - u_num_mat) > int_eps;
    u_num_mat(to_update) = u_num_mat_fine(to_update);
    rho_fine_IE = rho_fine_IE + 1;

    sum(to_update, 'all')
    rho_fine_IE
end

%% Evaluating points in s_patches
% refine one more time
rho_coarse_IE = rho_fine_IE;
gr_phi_coarse = gr_phi_fine;
curve_param_coarse = curve_param_fine;

rho_fine_IE = rho_coarse_IE + 1;
[gr_phi_fine, curve_param_fine] = refine_gr_phi(curve_param_rho1, rho_fine_IE, gr_phi_fft_rho1);

disp('s patches')
% excludes boundary
for eta_idx = M:-1:2
    for i = 1:curve_seq.n_curves
        s_patch = s_patches{i};
        [patch_X_s, patch_Y_s] = s_patch.xy_mesh;
        for xi_idx = 1:s_patch.n_xi
            x = patch_X_s(eta_idx, xi_idx); y = patch_Y_s(eta_idx, xi_idx);
            u_num_coarse_val = IE_curve_seq.u_num(x, y, curve_param_coarse, gr_phi_coarse);
            u_num_fine_val = IE_curve_seq.u_num(x, y, curve_param_fine, gr_phi_fine);
            
            while abs(u_num_coarse_val - u_num_fine_val) > int_eps
                rho_coarse_IE = rho_fine_IE;
                curve_param_coarse = curve_param_fine;
                gr_phi_coarse = gr_phi_fine;

                rho_fine_IE = rho_coarse_IE + 1
                [gr_phi_fine, curve_param_fine] = refine_gr_phi(curve_param_rho1, rho_fine_IE, gr_phi_fft_rho1);

                u_num_coarse_val = u_num_fine_val;
                u_num_fine_val = IE_curve_seq.u_num(x, y, curve_param_fine, gr_phi_fine);
            end
            s_patch.f_XY(eta_idx, xi_idx) = u_num_coarse_val;
        end
    end
end

% boundary value computation
for i = 1:curve_seq.n_curves
    s_patch = s_patches{i};
    [patch_X, patch_Y] = s_patch.xy_mesh;

    s_patch.f_XY(1, :) = u_G(patch_X(1, :), patch_Y(1, :));
end

for i = 1:curve_seq.n_curves
    s_patch = s_patches{i};
    in_patch = s_patch_msks{i};
    R_s_idxs = R.R_idxs(in_patch);

    [P_xi_s, P_eta_s] = R_xi_eta_inversion(R, s_patch, in_patch);

    for idx = 1:length(R_s_idxs)
        if s_patch.in_patch(P_xi_s(idx), P_eta_s(idx))
            u_num_mat(R_s_idxs(idx)) = s_patch.locally_compute(P_xi_s(idx), P_eta_s(idx), M);
        else
            s_patch_msks{i}(R_s_idxs(idx)) = false;
        end
    end
end

%% Evaluating points in corner patches using refinements
c_pts_msk = R.in_interior & ~well_interior_msk;

for i = 1:length(s_patch_msks)
    c_pts_msk = c_pts_msk & ~s_patch_msks{i};
end

R_idxs_c_left = R.R_idxs(c_pts_msk);

% estimate distance to boundary
dist_to_boundary = zeros(size(R_idxs_c_left));
for i = 1:length(R_idxs_c_left)
    dist_to_boundary(i) = sqrt(min((R.boundary_X-R.R_X(R_idxs_c_left(i))).^2+(R.boundary_Y-R.R_Y(R_idxs_c_left(i))).^2));
end

% compute points in decreasing distance order
[~, trav_order] = sort(dist_to_boundary, 'descend');
R_idxs_c_left = R_idxs_c_left(trav_order);

for i = 1:size(R_idxs_c_left, 1)
    x = R.R_X(R_idxs_c_left(i)); y = R.R_Y(R_idxs_c_left(i));

    u_num_coarse_val = IE_curve_seq.u_num(x, y, curve_param_coarse, gr_phi_coarse);
    u_num_fine_val = IE_curve_seq.u_num(x, y, curve_param_fine, gr_phi_fine);

    while abs(u_num_coarse_val - u_num_fine_val) > int_eps
        rho_coarse_IE = rho_fine_IE;
        curve_param_coarse = curve_param_fine;
        gr_phi_coarse = gr_phi_fine;

        rho_fine_IE = rho_coarse_IE + 1
        [gr_phi_fine, curve_param_fine] = refine_gr_phi(curve_param_rho1, rho_fine_IE, gr_phi_fft_rho1);

        u_num_coarse_val = u_num_fine_val;
        u_num_fine_val = IE_curve_seq.u_num(x, y, curve_param_fine, gr_phi_fine);
    end
    u_num_mat(R_idxs_c_left(i)) = u_num_coarse_val;
end

u_num_mat(~R.in_interior) = nan;
end

function [gr_phi_rho, curve_param_rho] = refine_gr_phi(curve_param_rho1, rho_IE, gr_phi_fft_rho1)
    curve_param_rho = curve_param_obj(curve_param_rho1.curve_n * rho_IE);
    
    padded_fft_coeffs = [zeros(ceil((curve_param_rho1.n_total*rho_IE-curve_param_rho1.n_total)/2), 1); gr_phi_fft_rho1; zeros(floor((curve_param_rho1.n_total*rho_IE-curve_param_rho1.n_total)/2), 1)];
    gr_phi_rho = rho_IE*curve_param_rho1.n_total*real(ifft(ifftshift(padded_fft_coeffs)));
end

function [well_interior_msk, s_patch_msks] = gen_R_msks(R, n_r, s_patches, c_0_patches, c_1_patches)
    well_interior_msk = R.in_interior;
    s_patch_msks = cell(size(s_patches));
    
    for i = 1:length(s_patch_msks)
        %S_patch
        s_patch = s_patches{i};
        [bound_X, bound_Y] = s_patch.boundary_mesh_xy(n_r, false);
        in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
        s_patch_msks{i} = in_patch;
        well_interior_msk = well_interior_msk & ~in_patch;

        c_0_patch = c_0_patches{i};
        c_1_patch = c_1_patches{i};
        
        if isobject(c_0_patch)
            [bound_X, bound_Y] = c_0_patch.boundary_mesh_xy(n_r, false);
            in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
            well_interior_msk = well_interior_msk & ~in_patch;
        end
        
        if  isobject(c_1_patch)
            [bound_X, bound_Y] = c_1_patch.boundary_mesh_xy(n_r, false);
            in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
            well_interior_msk = well_interior_msk & ~in_patch;
        end
    end

    for i = 1:length(s_patch_msks)
        M = s_patches{i}.n_eta;
        corner = s_patches{i}.M_p(1, 0)';
        well_interior_msk = well_interior_msk & ((corner(1)-R.R_X).^2 + (corner(2)-R.R_Y).^2 > (M*R.h).^2);
        s_patch_msks{i} = s_patch_msks{i} & ((corner(1)-R.R_X).^2 + (corner(2)-R.R_Y).^2 > (M*R.h).^2);
        s_patch_msks{mod(i, length(s_patch_msks))+1} = s_patch_msks{mod(i, length(s_patch_msks))+1} & ((corner(1)-R.R_X).^2 + (corner(2)-R.R_Y).^2 > (M*R.h).^2);
    end
end
