classdef Curve_seq_obj < handle
    properties
        first_curve
        last_curve
        n_curves
    end
    
    methods
        function obj = Curve_seq_obj()
            obj.first_curve = [];
            obj.last_curve = [];
            obj.n_curves = 0;
        end
        
        function add_curve(obj, l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n, frac_n_C_0, frac_n_C_1, frac_n_S_0, frac_n_S_1, h_norm)
            new_curve = Curve_obj(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n, frac_n_C_0, frac_n_C_1, frac_n_S_0, frac_n_S_1, h_norm, obj.first_curve);
            if isempty(obj.first_curve)
                obj.first_curve = new_curve;
                obj.last_curve = new_curve;
            end
            
            obj.last_curve.next_curve = new_curve;
            obj.last_curve = new_curve;            
            obj.n_curves = obj.n_curves+1;
        end
        
        function patches = construct_patches(obj, f, d, eps_xi_eta, eps_xy)
            curr = obj.first_curve;
            prev_S_patch = nan;
            prev_C_patch = nan;
            patches = cell(obj.n_curves*2, 1);
            
            for i = 1:obj.n_curves
                S_patch = curr.construct_S_patch(f, d, eps_xi_eta, eps_xy);
                C_patch = curr.construct_C_patch(f, d, eps_xi_eta, eps_xy);

                if i ~= 1
                    prev_C_patch.apply_w_W(prev_S_patch)
                    prev_C_patch.apply_w_L(S_patch)
                end
                
                patches{2*i-1} = S_patch;
                patches{2*i} = C_patch;
                
                prev_S_patch = S_patch;
                prev_C_patch = C_patch;
                curr = curr.next_curve;
            end
            C_patch.apply_w_W(S_patch);
            C_patch.apply_w_L(patches{1});
        end
        
        function plot_geometry(obj, d)
            curr = obj.first_curve;
            
            figure;
            for i = 1:obj.n_curves
                S_patch = curr.construct_S_patch(@(x, y) zeros(size(x)), d, nan, nan);
                C_patch = curr.construct_C_patch(@(x, y) zeros(size(x)), d, nan, nan);
                
                [boundary_X, boundary_Y] = S_patch.Q.boundary_mesh_xy(1, false);
                plot(boundary_X(:), boundary_Y(:));%,patches{2*i-1}.Q.f_XY(:));
                
                hold on;
                [boundary_X, boundary_Y] = C_patch.L.boundary_mesh_xy(1, false);
                plot(boundary_X(:), boundary_Y(:));%,patches{2*i}.L.f_XY);
                
                [boundary_X, boundary_Y] = C_patch.W.boundary_mesh_xy(1, false);
                plot(boundary_X(:), boundary_Y(:));

                curr = curr.next_curve;
            end
            [boundary_X, boundary_Y] = obj.construct_boundary_mesh(1);
            plot(boundary_X, boundary_Y);
        end        
        
        function [boundary_X, boundary_Y] = construct_boundary_mesh(obj, n_r)
            curr = obj.first_curve;
            n_points = 0;
            for i = 1:obj.n_curves
                n_points = n_points + (curr.n-1)*n_r + 1;
                curr = curr.next_curve;
            end
                        
            boundary_X = zeros(n_points, 1);
            boundary_Y = zeros(n_points, 1);
            curr_idx = 1;
            curr = obj.first_curve;
            for i = 1:obj.n_curves
                boundary_X(curr_idx:curr_idx+(curr.n-1)*n_r) = curr.l_1(linspace(0, 1, (curr.n-1)*n_r + 1)');
                boundary_Y(curr_idx:curr_idx+(curr.n-1)*n_r) = curr.l_2(linspace(0, 1, (curr.n-1)*n_r + 1)');
                
                curr_idx = curr_idx + (curr.n-1)*n_r + 1;
                curr = curr.next_curve;                
            end
        end
    end
end

