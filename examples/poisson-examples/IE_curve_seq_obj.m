classdef IE_curve_seq_obj < handle
    %IE_CURVE_SEQ_OBJ Obj for integral equation solver
    
    properties
        n_curves
        first_curve
        last_curve
    end
    
    methods
        function obj = IE_curve_seq_obj(curve_seq, p)
            %IE_CURVE_SEQ_OBJ Construct an instance of this class
            
            obj.n_curves = curve_seq.n_curves;
            obj.first_curve = IE_curve_obj(curve_seq.first_curve, 1, p);
            
            curr_curve = curve_seq.first_curve;
            curr_ie_curve = obj.first_curve;
            for i = 2:obj.n_curves
                curr_curve= curr_curve.next_curve;
                curr_ie_curve.next_curve = IE_curve_obj(curr_curve, i, p);
                curr_ie_curve = curr_ie_curve.next_curve;
            end
            obj.last_curve = curr_ie_curve;
            curr_ie_curve.next_curve = obj.first_curve;
        end
        
        function [A, b] = construct_A_b(obj, curve_param, f)
            curve_n = curve_param.curve_n;
            n_total = curve_param.n_total;
            start_idx = curve_param. start_idx;
            end_idx = curve_param.end_idx;
            
            A = zeros(n_total, n_total);
            b = zeros(n_total, 1);


            curr_y = obj.first_curve;
            for i = 1:obj.n_curves
                [s_mesh_y, ds_y] = curr_y.s_mesh(curve_n(i));
                theta_mesh_y = curr_y.w(s_mesh_y);
                b(start_idx(i):end_idx(i)) = curr_y.w_prime(s_mesh_y).*f(curr_y.l_1(theta_mesh_y), curr_y.l_2(theta_mesh_y));

                curr_x = obj.first_curve;
                for j = 1:obj.n_curves
                    [s_mesh_x, ds_x] = curr_x.s_mesh(curve_n(j));
                    theta_mesh_x = curr_x.w(s_mesh_x);
                    [theta_2, theta_1] = meshgrid(theta_mesh_y, theta_mesh_x);

                    A_local = curr_x.w_prime(s_mesh_x).*curr_x.K_boundary(theta_1, theta_2, curr_y).*sqrt(curr_y.l_1_prime(theta_mesh_y').^2+curr_y.l_2_prime(theta_mesh_y').^2)*ds_y;

                    if i == j
                        msk_diagonals = logical(diag(ones(curve_n(j), 1)));
                        A_local(msk_diagonals) = curr_x.w_prime(s_mesh_x).*curr_x.K_boundary_same_point(theta_mesh_y).*sqrt(curr_x.l_1_prime(theta_mesh_x).^2+curr_x.l_2_prime(theta_mesh_x).^2)*ds_x + 1/2;
                    end

                    A(start_idx(j):end_idx(j), start_idx(i):end_idx(i)) = A_local;

                    curr_x = curr_x.next_curve;
                end

                curr_y = curr_y.next_curve;
            end
            A(isnan(A)) = 0;
        end

        function [A, b] = construct_A_b_unstable(obj, curve_param, f)
            curve_n = curve_param.curve_n;
            n_total = curve_param.n_total;
            start_idx = curve_param. start_idx;
            end_idx = curve_param.end_idx;
            
            A = zeros(n_total, n_total);
            b = zeros(n_total, 1);


            curr_y = obj.first_curve;
            for i = 1:obj.n_curves
                [s_mesh_y, ds_y] = curr_y.s_mesh(curve_n(i));
                theta_mesh_y = curr_y.w(s_mesh_y);
                b(start_idx(i):end_idx(i)) = f(curr_y.l_1(theta_mesh_y), curr_y.l_2(theta_mesh_y));
        
                curr_x = obj.first_curve;
                for j = 1:obj.n_curves
                    [s_mesh_x, ds_x] = curr_x.s_mesh(curve_n(j));
                    theta_mesh_x = curr_x.w(s_mesh_x);
                    [theta_2, theta_1] = meshgrid(theta_mesh_y, theta_mesh_x);
        
                    A_local = curr_y.w_prime(s_mesh_y').*curr_x.K_boundary(theta_1, theta_2, curr_y).*sqrt(curr_y.l_1_prime(theta_mesh_y').^2+curr_y.l_2_prime(theta_mesh_y').^2)*ds_y;
        
                    if i == j
                        msk_diagonals = logical(diag(ones(curve_n(j), 1)));
                        A_local(msk_diagonals) = curr_y.w_prime(s_mesh_y).*curr_x.K_boundary_same_point(theta_mesh_y).*sqrt(curr_x.l_1_prime(theta_mesh_x).^2+curr_x.l_2_prime(theta_mesh_x).^2)*ds_x + 1/2;
                    end
        
                    A(start_idx(j):end_idx(j), start_idx(i):end_idx(i)) = A_local;
        
                    curr_x = curr_x.next_curve;
                end
        
                curr_y = curr_y.next_curve;
            end
            A(isnan(A)) = 0;
        end

        function u_num = u_num(obj, x, y, curve_param, gr_phi)
           curr = obj.first_curve;
           u_num = 0;
           for i = 1:obj.n_curves
                u_num =  u_num + curr.u_num_curve(x, y, curve_param, gr_phi);
                curr = curr.next_curve;
            end
        end
        
        function [s_patches, c_0_patches, c_1_patches] = construct_interior_patches(obj, curve_param, h_norm, M, eps_xi_eta, eps_xy)
            % Computing theta thresholds
            curr = obj.first_curve;
            for i = 1:obj.n_curves
                curr_n = curve_param.curve_n(i);
                next_n = curve_param.curve_n(mod(i, obj.n_curves)+1);
                curr_v = [curr.l_1(1); curr.l_2(1)] - [curr.l_1(1-1/curr_n); curr.l_2(1-1/curr_n)];
                next_v = [curr.next_curve.l_1(1/next_n); curr.next_curve.l_2(1/next_n)] - [curr.l_1(1); curr.l_2(1)];

                % C2-type patch means cross product is positive
                if curr_v(1)*next_v(2) - curr_v(2)*next_v(1) >= 0
                    curr.C2_corner = true;
                    [theta_1, theta_2] = compute_normal_intersection(curr, curr.next_curve, M, h_norm, eps_xy, [1; 0]);
                    curr.c_1_theta_thresh = floor(theta_1 * curr_n) / curr_n;
                    curr.next_curve.c_0_theta_thresh =  ceil(theta_2 * next_n) / next_n;
                 % C1-type patch otherwise
                else
                    curr.C2_corner = false;
                    curr.c_1_theta_thresh = 1;
                    curr.next_curve.c_0_theta_thresh = 0;
                end
                curr = curr.next_curve;
            end
            
            s_patches = cell(obj.n_curves, 1);
            c_0_patches = cell(obj.n_curves, 1);
            c_1_patches = cell(obj.n_curves, 1);
            curr = obj.first_curve;
            for i = 1:obj.n_curves
                [s_patches{i}, c_0_patches{i}, c_1_patches{i}] = curr.construct_interior_patch(curve_param, h_norm, M, eps_xi_eta, eps_xy);
                curr = curr.next_curve;
            end
        end
        
        function int_num = int_segment_general(obj, curve_param, x, y, gr_phi, curve_idx_interval, idx_interval)
            curr = obj.first_curve;
            
            % go to first curve
            for i = 2:curve_idx_interval(1)
                curr = curr.next_curve;
            end
            
            local_start_idx = idx_interval(1);
            
            % wraps arround
            if curve_idx_interval(1) > curve_idx_interval(2) || (curve_idx_interval(1) == curve_idx_interval(2) && idx_interval(2) < idx_interval(1))
                curve_end_idx = curve_idx_interval(2) + obj.n_curves;
            else
                curve_end_idx = curve_idx_interval(2);
            end
            
            i = curve_idx_interval(1);
            int_num = 0;
            
            local_end_idx = curve_param.curve_n(curr.curve_idx);
            while i < curve_end_idx
                int_num = int_num + curr.int_segment_general(curve_param, x, y, gr_phi, [local_start_idx; local_end_idx]);
                
                curr = curr.next_curve;
                i = i+1;
                
                local_start_idx = 1;
                local_end_idx = curve_param.curve_n(curr.curve_idx);
            end
            
            local_end_idx = idx_interval(2);
            int_num = int_num + curr.int_segment_general(curve_param, x, y, gr_phi, [local_start_idx; local_end_idx]);
        end
        
        function int_num = int_segment_boundary(obj, curve_param, theta_target, curve_target, gr_phi, curve_idx_interval, idx_interval)
            curr = obj.first_curve;
            
            % go to first curve
            for i = 2:curve_idx_interval(1)
                curr = curr.next_curve;
            end
            
            local_start_idx = idx_interval(1);
            
            % wraps arround
            if curve_idx_interval(1) > curve_idx_interval(2) || (curve_idx_interval(1) == curve_idx_interval(2) && idx_interval(2) < idx_interval(1))
                curve_end_idx = curve_idx_interval(2) + obj.n_curves;
            else
                curve_end_idx = curve_idx_interval(2);
            end
            
            i = curve_idx_interval(1);
            int_num = 0;
            
            local_end_idx = curve_param.curve_n(curr.curve_idx);
            while i < curve_end_idx
                int_num = int_num + curr.int_segment_boundary(curve_param, theta_target, curve_target, gr_phi, [local_start_idx; local_end_idx]);
                
                curr = curr.next_curve;
                i = i+1;
                
                local_start_idx = 1;
                local_end_idx = curve_param.curve_n(curr.curve_idx);
            end
            
            local_end_idx = idx_interval(2);
            int_num = int_num + curr.int_segment_boundary(curve_param, theta_target, curve_target, gr_phi, [local_start_idx; local_end_idx]);
        end
        
        function [interp_int_seg] = interp_int_seg(obj, xi, eta, patch, target_curve, curve_param, gr_phi, curve_idx_interval, idx_interval, u_G, gr_phi_C, curve_idx_C, idx_interval_C)            
            eta_mesh = patch.eta_mesh;
            
            interpol_nodes = patch.M_p(xi, eta_mesh);
            interpol_vals = zeros(size(eta_mesh));
            for i = 2:length(interpol_vals)
                interpol_vals(i) = obj.int_segment_general(curve_param, interpol_nodes(i, 1), interpol_nodes(i, 2), gr_phi, curve_idx_interval, idx_interval);
            end
            interpol_vals(1) = u_G(interpol_nodes(1, 1), interpol_nodes(1, 2)) - obj.int_segment_boundary(curve_param, xi, target_curve, gr_phi_C, curve_idx_C, idx_interval_C);
            
            interp_int_seg = barylag([eta_mesh, interpol_vals], eta);
        end
    end
end


function  [theta_1, theta_2] = compute_normal_intersection(curve_1, curve_2, M, h_norm, eps_xy, initial_guess)
    nu_norm = @(theta, curve) sqrt(curve.l_1_prime(theta).^2 + curve.l_2_prime(theta).^2);
    err_guess_x = @(theta_1, theta_2) curve_1.l_1(theta_1) - (M-1)*h_norm*(curve_1.l_2_prime(theta_1)./nu_norm(theta_1, curve_1)) - (curve_2.l_1(theta_2) - (M-1)*h_norm*(curve_2.l_2_prime(theta_2)./nu_norm(theta_2, curve_2)));
    err_guess_y = @(theta_1, theta_2) curve_1.l_2(theta_1) + (M-1)*h_norm*(curve_1.l_1_prime(theta_1)./nu_norm(theta_1, curve_1)) - (curve_2.l_2(theta_2) + (M-1)*h_norm*(curve_2.l_1_prime(theta_2)./nu_norm(theta_2, curve_2)));

    derr_x_d1 = @(theta_1, theta_2) curve_1.l_1_prime(theta_1)-(M-1)*h_norm*(curve_1.l_2_dprime(theta_1).*nu_norm(theta_1, curve_1).^2-curve_1.l_2_prime(theta_1).*(curve_1.l_2_dprime(theta_1).*curve_1.l_2_prime(theta_1)+curve_1.l_1_dprime(theta_1).*curve_1.l_1_prime(theta_1)))./nu_norm(theta_1, curve_1).^3;
    derr_y_d1 = @(theta_1, theta_2) curve_1.l_2_prime(theta_1)+(M-1)*h_norm*(curve_1.l_1_dprime(theta_1).*nu_norm(theta_1, curve_1).^2-curve_1.l_1_prime(theta_1).*(curve_1.l_2_dprime(theta_1).*curve_1.l_2_prime(theta_1)+curve_1.l_1_dprime(theta_1).*curve_1.l_1_prime(theta_1)))./nu_norm(theta_1, curve_1).^3;
    derr_x_d2 = @(theta_1, theta_2) -(curve_2.l_1_prime(theta_2)-(M-1)*h_norm*(curve_2.l_2_dprime(theta_2).*nu_norm(theta_2, curve_2).^2-curve_2.l_2_prime(theta_2).*(curve_2.l_2_dprime(theta_2).*curve_2.l_2_prime(theta_2)+curve_2.l_1_dprime(theta_2).*curve_2.l_1_prime(theta_2)))./nu_norm(theta_2, curve_2).^3);
    derr_y_d2 = @(theta_1, theta_2) -( curve_2.l_2_prime(theta_2)+(M-1)*h_norm*(curve_2.l_1_dprime(theta_2).*nu_norm(theta_2, curve_2).^2-curve_2.l_1_prime(theta_2).*(curve_2.l_2_dprime(theta_2).*curve_2.l_2_prime(theta_2)+curve_2.l_1_dprime(theta_2).*curve_2.l_1_prime(theta_2)))./nu_norm(theta_2, curve_2).^3);

    err_guess = @(v) [err_guess_x(v(1), v(2)); err_guess_y(v(1), v(2))];
    J_err = @(v) [derr_x_d1(v(1), v(2)) derr_x_d2(v(1), v(2)); derr_y_d1(v(1), v(2)) derr_y_d2(v(1), v(2))];

    [v_guess, converged] = newton_solve(err_guess, J_err, initial_guess, eps_xy, 100);
    theta_1 = v_guess(1); theta_2 = v_guess(2);
    if ~converged
        warning('Nonconvergence in normal-boundary intersection')
    end
end

