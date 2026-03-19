function FC2D_patches(interior_patches, R, f_interior, curve_seq, d, C_S, n_r, A_S, Q_S, C_C, A_C, Q_C, M)
%FC2D Summary of this function goes here
%   Detailed explanation goes here    
    FC_patches = cell(4*curve_seq.n_curves, 1);
    for i = 1:curve_seq.n_curves
        FC_patches{4*i-3} = interior_patches{2*i-1}.FC(C_S, n_r, d, A_S, Q_S);
        [FC_patches{4*i-2}, FC_patches{4*i-1}, FC_patches{4*i}] = interior_patches{2*i}.FC(C_C, n_r, d, A_C, Q_C, M);
    end

    for i=1:length(FC_patches)
        R.interpolate_patch(FC_patches{i}, n_r, M);
    end

    R.f_R(R.in_interior) = f_interior;
    
    R.compute_fc_coeffs();
end


