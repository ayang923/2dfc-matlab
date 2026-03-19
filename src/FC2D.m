function [R, interior_patches, FC_patches, fc_err] = FC2D(f, h, curve_seq, eps_xi_eta, eps_xy, d, C_S, n_r, A_S, Q_S, C_C, A_C, Q_C, M)
%FC2D Summary of this function goes here
%   Detailed explanation goes here
    interior_patches = curve_seq.construct_patches(f, d, eps_xi_eta, eps_xy);
    
    x_min = curve_seq.first_curve.l_1(0);
    x_max = curve_seq.first_curve.l_1(0);
    y_min = curve_seq.first_curve.l_2(0);
    y_max = curve_seq.first_curve.l_2(0);
    
    FC_patches = cell(4*curve_seq.n_curves, 1);
    for i = 1:curve_seq.n_curves
        FC_patches{4*i-3} = interior_patches{2*i-1}.FC(C_S, n_r, d, A_S, Q_S);
        [FC_patches{4*i-2}, FC_patches{4*i-1}, FC_patches{4*i}] = interior_patches{2*i}.FC(C_C, n_r, d, A_C, Q_C, M);
        for j = 0:3
            x_min = min(FC_patches{4*i-j}.x_min, x_min);
            x_max = max(FC_patches{4*i-j}.x_max, x_max);
            y_min = min(FC_patches{4*i-j}.y_min, y_min);
            y_max = max(FC_patches{4*i-j}.y_max, y_max);
        end
    end
    
    % for i = 1:length(FC_patches)
    %     [boundary_x, boundary_y] = FC_patches{i}.boundary_mesh_xy(1, false); plot(boundary_x, boundary_y)
    % end
    [boundary_X, boundary_Y] = curve_seq.construct_boundary_mesh(n_r*20);
    R = R_cartesian_mesh_obj(x_min-rand(1)*h, x_max+rand(1)*h, y_min-rand(1)*h, y_max+rand(1)*h, h, boundary_X, boundary_Y);

    for i=1:length(FC_patches)
        R.interpolate_patch(FC_patches{i}, n_r, M);
    end

    R.fill_interior(f);
    
    R.compute_fc_coeffs();
    [R_X_err, R_Y_err, f_interpolation, interior_idx] = R.ifft_interpolation(2);
    f_exact = f(R_X_err, R_Y_err);
    fc_err = max(abs(f_exact(interior_idx) - f_interpolation(interior_idx)));
    disp(['abs max error: ', num2str(fc_err)]);

    fc_err = max(abs(f_exact(interior_idx) - f_interpolation(interior_idx)))/max(abs(f_exact(interior_idx)));
    disp(['rel max error: ', num2str(fc_err)]);

    fc_err = sqrt(sum((f_exact(interior_idx) - f_interpolation(interior_idx)).^2)./sum((f_exact(interior_idx).^2)));
    disp(['rel 2 error: ', num2str(fc_err)])
end

