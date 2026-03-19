classdef IE_curve_obj < handle
    %IE_CURVE_OBJ Curve obj for integral equation solver
    
    properties
        l_1
        l_2
        l_1_prime
        l_2_prime
        l_1_dprime
        l_2_dprime
        
        w
        w_prime
                
        c_0_theta_thresh
        c_1_theta_thresh
        C2_corner
        
        curve_idx
        next_curve
    end
    
    methods
        function obj = IE_curve_obj(curve, curve_idx, p)
            %IE_CURVE_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj.l_1 = curve.l_1;
            obj.l_2 = curve.l_2;
            obj.l_1_prime = curve.l_1_prime;
            obj.l_2_prime = curve.l_2_prime;
            obj.l_1_dprime = curve.l_1_dprime;
            obj.l_2_dprime = curve.l_2_dprime;
                        
            obj.curve_idx = curve_idx;
            
            v = @(s) (1/p-1/2)*(1-2*s).^3+1/p*(2*s-1)+1/2;
            v_prime = @(s) 2/p-6*(1/p-1/2)*(1-2*s).^2;

            obj.w = @(s) v(s).^p./(v(s).^p+v(1-s).^p);
            obj.w_prime = @(s) p*((v_prime(s).*v(s).^(p-1))./(v(s).^p+v(1-s).^p)-(v(s).^(p-1).*v_prime(s)-v(1-s).^(p-1).*v_prime(1-s)).*v(s).^p./(v(s).^p+v(1-s).^p).^2);
        end
        
        function [s_mesh, ds] = s_mesh(obj, n)
            ds = 1/n;
            s_mesh = linspace(0, 1-ds, n)';
        end
        
        function theta_mesh = theta_mesh(obj, n)
            [s_mesh, ~] = obj.s_mesh(n);
            theta_mesh = obj.w(s_mesh);
        end
        
        function u_num_curve = u_num_curve(obj, x, y, curve_param, gr_phi)
            curve_n = curve_param.curve_n(obj.curve_idx);
            start_idx = curve_param.start_idx(obj.curve_idx);
            end_idx = curve_param.end_idx(obj.curve_idx);
            [s_mesh, ds] =  obj.s_mesh(curve_n);
            theta_mesh = obj.w(s_mesh);

            u_num_curve = transpose(obj.K_general([x; y], theta_mesh).*gr_phi(start_idx:end_idx).*sqrt(obj.l_1_prime(theta_mesh).^2+obj.l_2_prime(theta_mesh).^2)) * ones(curve_n, 1) * ds;
        end
        
        function [s_patch, c_0_patch, c_1_patch] = construct_interior_patch(obj, curve_param, h_norm, M, eps_xi_eta, eps_xy)
            curve_n = curve_param.curve_n(obj.curve_idx);
            xi_diff = 1;
            xi_tilde = @(xi) xi;

            nu_norm = @(theta) sqrt(obj.l_1_prime(theta).^2 + obj.l_2_prime(theta).^2);

            M_p_1 = @(xi, eta) obj.l_1(xi_tilde(xi)) - eta.*obj.l_2_prime(xi_tilde(xi))./nu_norm(xi_tilde(xi));
            M_p_2 = @(xi, eta) obj.l_2(xi_tilde(xi)) + eta.*obj.l_1_prime(xi_tilde(xi))./nu_norm(xi_tilde(xi));
            M_p = @(xi, eta) [M_p_1(xi, eta), M_p_2(xi, eta)];

            dM_p_1_dxi = @(xi, eta) xi_diff * (obj.l_1_prime(xi_tilde(xi))-eta*(obj.l_2_dprime(xi_tilde(xi)).*nu_norm(xi_tilde(xi)).^2-obj.l_2_prime(xi_tilde(xi)).*(obj.l_2_dprime(xi_tilde(xi)).*obj.l_2_prime(xi_tilde(xi))+obj.l_1_dprime(xi_tilde(xi)).*obj.l_1_prime(xi_tilde(xi))))./nu_norm(xi_tilde(xi)).^3);
            dM_p_2_dxi = @(xi, eta) xi_diff * (obj.l_2_prime(xi_tilde(xi))+eta*(obj.l_1_dprime(xi_tilde(xi)).*nu_norm(xi_tilde(xi)).^2-obj.l_1_prime(xi_tilde(xi)).*(obj.l_2_dprime(xi_tilde(xi)).*obj.l_2_prime(xi_tilde(xi))+obj.l_1_dprime(xi_tilde(xi)).*obj.l_1_prime(xi_tilde(xi))))./nu_norm(xi_tilde(xi)).^3);
            dM_p_1_deta = @(xi, eta) -obj.l_2_prime(xi_tilde(xi)) ./ nu_norm(xi_tilde(xi));
            dM_p_2_deta = @(xi, eta) obj.l_1_prime(xi_tilde(xi)) ./ nu_norm(xi_tilde(xi));

            J = @(v) [dM_p_1_dxi(v(1), v(2)), dM_p_1_deta(v(1), v(2)); dM_p_2_dxi(v(1), v(2)), dM_p_2_deta(v(1), v(2))];

            n_xi_interior = max(round((obj.c_1_theta_thresh - obj.c_0_theta_thresh)*curve_n) + 1, M);
            s_patch = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi_interior, M, obj.c_0_theta_thresh, obj.c_1_theta_thresh, 0, (M-1)*h_norm, zeros(M, n_xi_interior));

            if obj.c_0_theta_thresh ~= 0
                n_xi_corner =  max(round(obj.c_0_theta_thresh*curve_n) + 1, M);
                c_0_patch = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi_corner, M, 0, obj.c_0_theta_thresh, 0, (M-1)*h_norm, zeros(M, n_xi_corner));
            else
                c_0_patch = nan;
            end
            if obj.c_1_theta_thresh ~= 1
                n_xi_corner = max(round((1-obj.c_1_theta_thresh)*curve_n) + 1, M);
                c_1_patch = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi_corner, M, obj.c_1_theta_thresh, 1, 0, (M-1)*h_norm, zeros(M, n_xi_corner));
            else
                c_1_patch = nan;
            end
        end
        
        function int_num = int_segment_general(obj, curve_param, x, y, gr_phi, idx_interval)
            [s_mesh, ds] = obj.s_mesh(curve_param.curve_n(obj.curve_idx)); s_mesh = s_mesh(idx_interval(1):idx_interval(2));
            theta_mesh = obj.w(s_mesh);

            gr_phi_idxs = curve_param.start_idx(obj.curve_idx) - 1 + (idx_interval(1):idx_interval(2))';
            int_num = transpose(obj.K_general([x; y], theta_mesh).*gr_phi(gr_phi_idxs).*sqrt(obj.l_1_prime(theta_mesh).^2+obj.l_2_prime(theta_mesh).^2)) * ones(length(s_mesh), 1) * ds;
        end
        
        % Assuming s_target is not within this interval
        function int_num = int_segment_boundary(obj, curve_param, theta_target, curve_target, gr_phi, idx_interval)
            [s_mesh, ds] = obj.s_mesh(curve_param.curve_n(obj.curve_idx)); s_mesh = s_mesh(idx_interval(1):idx_interval(2));
            theta_mesh = obj.w(s_mesh);
            
            gr_phi_idxs = curve_param.start_idx(obj.curve_idx) - 1 + (idx_interval(1):idx_interval(2))';
            int_num = transpose(curve_target.K_boundary(theta_target, theta_mesh, obj).*gr_phi(gr_phi_idxs).*sqrt(obj.l_1_prime(theta_mesh).^2+obj.l_2_prime(theta_mesh).^2)) * ones(length(s_mesh), 1) * ds;
        end
        
        function [gr_phi_U] = c_0_apply_POU(obj, curve_param, U, gr_phi)
            theta_mesh = obj.theta_mesh(curve_param.curve_n(obj.curve_idx));
            U_start_idx = sum(theta_mesh < 2*obj.c_0_theta_thresh) + 1;
            U_end_idx = U_start_idx + length(U) - 1;
            curve_param.U_c_0_intervals(obj.curve_idx, :) = [U_start_idx, U_end_idx];
            
            gr_phi_start_idx = curve_param.start_idx(obj.curve_idx) + U_start_idx - 1;
            gr_phi_end_idx = curve_param.start_idx(obj.curve_idx) + U_end_idx - 1;     
            gr_phi_U = gr_phi; gr_phi_U(gr_phi_start_idx:gr_phi_end_idx) = gr_phi_U(gr_phi_start_idx:gr_phi_end_idx) .* (1-U);
        end
        
        function [gr_phi_U] = c_1_apply_POU(obj, curve_param, U, gr_phi)
            theta_mesh = obj.theta_mesh(curve_param.curve_n(obj.curve_idx));
            U_end_idx = sum(theta_mesh < 2*obj.c_1_theta_thresh-1);
            U_start_idx = U_end_idx - length(U) + 1;
            curve_param.U_c_1_intervals(obj.curve_idx, :) = [U_start_idx, U_end_idx];
            
            gr_phi_start_idx = curve_param.start_idx(obj.curve_idx) + U_start_idx - 1;
            gr_phi_end_idx = curve_param.start_idx(obj.curve_idx) + U_end_idx - 1;     
            gr_phi_U = gr_phi; gr_phi_U(gr_phi_start_idx:gr_phi_end_idx) = gr_phi_U(gr_phi_start_idx:gr_phi_end_idx) .* U;
        end

        function K_general_eval = K_general(obj, x, theta)
            K_general_eval = -1 / (2 * pi) * ( ...
                (x(1) - obj.l_1(theta)) .* obj.l_2_prime(theta) - ...
                (x(2) - obj.l_2(theta)) .* obj.l_1_prime(theta) ...
            ) ./ ( ...
                sqrt(obj.l_1_prime(theta).^2 + obj.l_2_prime(theta).^2) .* ...
                ( ...
                    (x(1) - obj.l_1(theta)).^2 + ...
                    (x(2) - obj.l_2(theta)).^2 ...
                ) ...
            );
        end
        
        function  K_boundary_same_point_eval = K_boundary_same_point(obj, theta)
            K_boundary_same_point_eval = -1 / (4 * pi) * ( ...
                obj.l_2_prime(theta) .* obj.l_1_dprime(theta) - ...
                obj.l_1_prime(theta) .* obj.l_2_dprime(theta) ...
            ) ./ ( ...
                (obj.l_1_prime(theta).^2 + obj.l_2_prime(theta).^2).^(3/2) ...
            );
        end
        
        function K_boundary_eval = K_boundary(obj, theta_1, theta_2, curve_2)
            K_boundary_eval = -1 / (2 * pi) * ( ...
                (obj.l_1(theta_1) - curve_2.l_1(theta_2)) .* curve_2.l_2_prime(theta_2) - ...
                (obj.l_2(theta_1) - curve_2.l_2(theta_2)) .* curve_2.l_1_prime(theta_2) ...
            ) ./ ( ...
                sqrt(curve_2.l_1_prime(theta_2).^2 + curve_2.l_2_prime(theta_2).^2) .* ...
                ( ...
                    (obj.l_1(theta_1) - curve_2.l_1(theta_2)).^2 + ...
                    (obj.l_2(theta_1) - curve_2.l_2(theta_2)).^2 ...
                ) ...
            );
        end
    end

end


