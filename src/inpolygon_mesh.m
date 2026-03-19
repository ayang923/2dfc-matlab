function [in] = inpolygon_mesh(R_X, R_Y, boundary_x, boundary_y)
    boundary_x_edge_1 = boundary_x(1:end-1);
    boundary_x_edge_2 = boundary_x(2:end);
    boundary_y_edge_1 = boundary_y(1:end-1);
    boundary_y_edge_2 = boundary_y(2:end);
    
    boundary_idxs = transpose(1:length(boundary_x_edge_1));
    
    x_start = R_X(1, 1);
    y_start = R_Y(1, 1);
    
    h_x = R_X(1, 2) - R_X(1, 1);
    h_y = R_Y(2, 1) - R_Y(1, 1);
    
    boundary_y_j = (boundary_y - y_start)/h_y;
    boundary_y_edge_1_j = boundary_y_j(1:end-1);
    boundary_y_edge_2_j = boundary_y_j(2:end);
    % rounding exact points if needed
    boundary_y_edge_1_j(abs(boundary_y_edge_1_j-round(boundary_y_edge_1_j)) < eps) = round(boundary_y_edge_1_j(abs(boundary_y_edge_1_j-round(boundary_y_edge_1_j)) < eps));
    boundary_y_edge_2_j(abs(boundary_y_edge_2_j-round(boundary_y_edge_2_j)) < eps) = round(boundary_y_edge_2_j(abs(boundary_y_edge_2_j-round(boundary_y_edge_2_j)) < eps));
    
    intersection_idxs = boundary_idxs(floor(boundary_y_edge_1_j) ~= floor(boundary_y_edge_2_j));

    in = false(size(R_X));
    for intersection_idx=intersection_idxs'
        x_edge_1 = boundary_x_edge_1(intersection_idx);
        x_edge_2 = boundary_x_edge_2(intersection_idx);
        y_edge_1 = boundary_y_edge_1(intersection_idx);
        y_edge_2 = boundary_y_edge_2(intersection_idx);
        y_edge_1_j = floor(boundary_y_edge_1_j(intersection_idx));
        y_edge_2_j = floor(boundary_y_edge_2_j(intersection_idx));
        
        intersection_mesh_y = ((min(y_edge_1_j, y_edge_2_j)+1):max(y_edge_1_j, y_edge_2_j))*h_y+y_start;
        intersection_x = x_edge_1 + (x_edge_2-x_edge_1).*(intersection_mesh_y-y_edge_1)./(y_edge_2-y_edge_1);

        mesh_intersection_idxs = sub2ind(size(in), round((intersection_mesh_y-y_start)/h_y)+1, floor((intersection_x-x_start)/h_x)+1);
        in(mesh_intersection_idxs) = ~in(mesh_intersection_idxs);
    end
    for vert_idx = 1:size(in, 1)
        in_interior = false;
        for horz_idx = 1:size(in, 2)
            if in(vert_idx, horz_idx) && ~in_interior
                in_interior = true;
                in(vert_idx, horz_idx) = false;
            elseif in(vert_idx, horz_idx) && in_interior
                in_interior = false;
            elseif in_interior
                in(vert_idx, horz_idx) = true;
            end
        end
    end
end

