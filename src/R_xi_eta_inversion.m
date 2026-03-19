function [P_xi, P_eta] = R_xi_eta_inversion(R, patch, in_patch)
%XI_ETA_INVERSION Summary of this function goes here
%   Detailed explanation goes here
    R_patch_idxs = R.R_idxs(in_patch);
    in_patch = double(in_patch);

    % computing initial "proximity map" with floor and ceil
    % operator
    [XI, ETA] = patch.xi_eta_mesh();
    [patch_X, patch_Y] = patch.convert_to_XY(XI, ETA);

    floor_X_j = floor((patch_X-R.x_start)/R.h);
    ceil_X_j = ceil((patch_X-R.x_start)/R.h);
    floor_Y_j = floor((patch_Y-R.y_start)/R.h);
    ceil_Y_j = ceil((patch_Y-R.y_start)/R.h);

    n_in_patch = sum(in_patch, 'all');            
    P_xi = zeros(n_in_patch, 1);
    P_eta = zeros(n_in_patch, 1);

    for i = 1:length(R_patch_idxs)
        P_xi(i) = nan; P_eta(i) = nan;
        in_patch(R_patch_idxs(i)) = i;
    end


    for i = 1:size(patch_X, 1)
        for j = 1:size(patch_X, 2)
            neighbors = [floor_X_j(i, j), floor_X_j(i, j), ceil_X_j(i, j), ceil_X_j(i, j); floor_Y_j(i, j), ceil_Y_j(i, j), floor_Y_j(i, j), ceil_Y_j(i, j)];

            for neighbor_i = 1:size(neighbors, 2)
                neighbor = neighbors(:, neighbor_i) + 1;
                if any(neighbor > [R.n_x; R.n_y]) || any(neighbor < [1; 1])
                    continue;
                end
                patch_idx = sub2ind([R.n_y, R.n_x], neighbor(2), neighbor(1));

                if (in_patch(patch_idx) ~= 0) && isnan(P_xi(in_patch(patch_idx)))
                    [xi, eta, converged] = patch.inverse_M_p((neighbor(1)-1)*R.h + R.x_start, (neighbor(2)-1)*R.h + R.y_start, [XI(i, j); ETA(i, j)]);
                    if converged
                        P_xi(in_patch(patch_idx)) = xi;
                        P_eta(in_patch(patch_idx)) = eta;
                    else
                        warning("Nonconvergence in interpolation")
                    end
                end
            end
        end
    end

    % second pass for points that aren't touched, could
    % theoretically modify so that we continuously do this until
    % all points are touched

    nan_set = containers.Map('KeyType', 'int64', 'ValueType', 'logical');
    % first iterate through and construct nan set            
    for i=1:length(P_xi)
        if isnan(P_xi(i))
            nan_set(i) = true;
        end
    end

    % pass through nan set until empty
    while nan_set.Count > 0
        for key = keys(nan_set)
            [i, j] = ind2sub([R.n_y, R.n_x], R_patch_idxs(key{1}));

            neighbor_shifts = [1 -1 0 0 1 1 -1 -1; 0 0 -1 1 1 -1 1 -1];
            is_touched = false;
            for neighbor_shift_i = 1:size(neighbor_shifts, 2)
                neighbor_shift = neighbor_shifts(:, neighbor_shift_i);
                neighbor = sub2ind([R.n_y, R.n_x], i+ neighbor_shift(1), j + neighbor_shift(2));

                if (in_patch(neighbor) ~= 0) && ~isnan(P_xi(in_patch(neighbor)))
                    [xi, eta, converged] = patch.inverse_M_p((j-1)*R.h+R.x_start, (i-1)*R.h+R.y_start, [P_xi(in_patch(neighbor)); P_eta(in_patch(neighbor))]);
                    if converged
                        P_xi(key{1}) = xi; P_eta(key{1}) = eta;
                        is_touched = true;
                        break;
                    else
                        warning("Nonconvergence in interpolation")
                    end
                end
            end
            if is_touched
                remove(nan_set, key{1});
            end
        end
    end
end

