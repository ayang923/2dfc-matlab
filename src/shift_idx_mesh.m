function idx_mesh = shift_idx_mesh(idx_mesh, min_bound, max_bound)
    if idx_mesh(1) < min_bound
            idx_mesh = idx_mesh + min_bound - idx_mesh(1);
    end
    if idx_mesh(end) > max_bound
            idx_mesh = idx_mesh + max_bound - idx_mesh(end);
    end
end
