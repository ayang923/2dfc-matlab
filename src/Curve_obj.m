classdef Curve_obj < handle
    properties
        l_1
        l_2
        l_1_prime
        l_2_prime
        l_1_dprime
        l_2_dprime
        
        n
        n_C_0
        n_C_1
        n_S_0
        n_S_1
        
        h_tan
        h_norm
        
        next_curve
    end
    
    methods
        function obj = Curve_obj(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n, frac_n_C_0, frac_n_C_1, frac_n_S_0, frac_n_S_1, h_norm, next_curve)
            obj.l_1 = l_1;
            obj.l_2 = l_2;
            obj.l_1_prime = l_1_prime;
            obj.l_2_prime = l_2_prime;
            obj.l_1_dprime = l_1_dprime;
            obj.l_2_dprime = l_2_dprime;
            
            if n == 0
                % use default value (h is around h_norm)
                L = obj.compute_length();
                obj.n = ceil(L ./ h_norm) + 1;
            else
                obj.n = n;
            end
            if frac_n_C_0 == 0
                obj.n_C_0 = ceil(1/10*obj.n);
            else
                obj.n_C_0 = ceil(frac_n_C_0*obj.n);
            end
            if frac_n_C_1 == 0
                obj.n_C_1 = ceil(1/10*obj.n);
            else
                obj.n_C_1 = ceil(frac_n_C_1*obj.n);
            end
            if frac_n_S_0 == 0
                obj.n_S_0 = ceil(2/3*obj.n_C_0);
            else
                obj.n_S_0 = ceil(frac_n_S_0*obj.n_C_0);
            end
            if frac_n_S_1 == 0
                obj.n_S_1 = ceil(2/3*obj.n_C_1);
            else
                obj.n_S_1 = ceil(frac_n_S_1*obj.n_C_1);
            end

            obj.h_tan = 1/(obj.n-1);
            obj.h_norm = h_norm;
            
            if isempty(next_curve)
                obj.next_curve = obj;
            else
                obj.next_curve = next_curve;
            end
        end
        
        function S_patch = construct_S_patch(obj, f, d, eps_xi_eta, eps_xy)
            xi_diff = 1-(obj.n_C_1-obj.n_S_1)*obj.h_tan - (obj.n_C_0-obj.n_S_0)*obj.h_tan;
            xi_0 = (obj.n_C_0 - obj.n_S_0) * obj.h_tan;
            xi_tilde = @(xi) xi_diff*xi +xi_0;
            
            nu_norm = @(theta) sqrt(obj.l_1_prime(theta).^2 + obj.l_2_prime(theta).^2);
            
            M_p_1 = @(xi, eta) obj.l_1(xi_tilde(xi)) - eta.*obj.l_2_prime(xi_tilde(xi))./nu_norm(xi_tilde(xi));
            M_p_2 = @(xi, eta) obj.l_2(xi_tilde(xi)) + eta.*obj.l_1_prime(xi_tilde(xi))./nu_norm(xi_tilde(xi));
            M_p = @(xi, eta) [M_p_1(xi, eta), M_p_2(xi, eta)];
            
            dM_p_1_dxi = @(xi, eta) xi_diff * (obj.l_1_prime(xi_tilde(xi))-eta*(obj.l_2_dprime(xi_tilde(xi)).*nu_norm(xi_tilde(xi)).^2-obj.l_2_prime(xi_tilde(xi)).*(obj.l_2_dprime(xi_tilde(xi)).*obj.l_2_prime(xi_tilde(xi))+obj.l_1_dprime(xi_tilde(xi)).*obj.l_1_prime(xi_tilde(xi))))./nu_norm(xi_tilde(xi)).^3);
            dM_p_2_dxi = @(xi, eta) xi_diff * (obj.l_2_prime(xi_tilde(xi))+eta*(obj.l_1_dprime(xi_tilde(xi)).*nu_norm(xi_tilde(xi)).^2-obj.l_1_prime(xi_tilde(xi)).*(obj.l_2_dprime(xi_tilde(xi)).*obj.l_2_prime(xi_tilde(xi))+obj.l_1_dprime(xi_tilde(xi)).*obj.l_1_prime(xi_tilde(xi))))./nu_norm(xi_tilde(xi)).^3);
            dM_p_1_deta = @(xi, eta) -obj.l_2_prime(xi_tilde(xi)) ./ nu_norm(xi_tilde(xi));
            dM_p_2_deta = @(xi, eta) obj.l_1_prime(xi_tilde(xi)) ./ nu_norm(xi_tilde(xi));
            
            J = @(v) [dM_p_1_dxi(v(1), v(2)), dM_p_1_deta(v(1), v(2)); dM_p_2_dxi(v(1), v(2)), dM_p_2_deta(v(1), v(2))];
            
            S_patch = S_patch_obj(M_p, J, obj.h_norm, eps_xi_eta, eps_xy, obj.n-(obj.n_C_1-obj.n_S_1)-(obj.n_C_0-obj.n_S_0), d, nan);
            [X, Y] = S_patch.Q.xy_mesh;
            S_patch.Q.f_XY = f(X, Y);
        end
        
        function C_patch = construct_C_patch(obj, f, d, eps_xi_eta, eps_xy)
            % compute angle between two curves to determine what type of
            % corner patch
            curr_v = [obj.l_1(1); obj.l_2(1)] - [obj.l_1(1-1/(obj.n-1)); obj.l_2(1-1/(obj.n-1))];
            next_v = [obj.next_curve.l_1(1/(obj.next_curve.n-1)); obj.next_curve.l_2(1/(obj.next_curve.n-1))] - [obj.l_1(1); obj.l_2(1)];
            
            % C2-type patch means cross product is positive
            if curr_v(1)*next_v(2) - curr_v(2)*next_v(1) >= 0
                xi_diff = -(obj.n_C_1-1)*obj.h_tan;
                eta_diff = (obj.next_curve.n_C_0-1)*obj.next_curve.h_tan;
                xi_0 = 1;
                eta_0 = 0;
                
                xi_tilde = @(xi) xi_diff*xi + xi_0;
                eta_tilde = @(eta) eta_diff*eta +eta_0;
                
                M_p_1 = @(xi, eta) obj.l_1(xi_tilde(xi)) + obj.next_curve.l_1(eta_tilde(eta)) - obj.l_1(1);
                M_p_2 = @(xi, eta) obj.l_2(xi_tilde(xi)) + obj.next_curve.l_2(eta_tilde(eta)) - obj.l_2(1);
                M_p = @(xi, eta) [M_p_1(xi, eta), M_p_2(xi, eta)];
                
                dM_p_1_dxi = @(xi, eta) xi_diff*obj.l_1_prime(xi_tilde(xi));
                dM_p_2_dxi = @(xi, eta) xi_diff*obj.l_2_prime(xi_tilde(xi));
                dM_p_1_deta = @(xi, eta) eta_diff*obj.next_curve.l_1_prime(eta_tilde(eta));
                dM_p_2_deta = @(xi, eta) eta_diff*obj.next_curve.l_2_prime(eta_tilde(eta));
                
                J = @(v) [dM_p_1_dxi(v(1), v(2)), dM_p_1_deta(v(1), v(2)); dM_p_2_dxi(v(1), v(2)), dM_p_2_deta(v(1), v(2))];
                
                C_patch = C2_patch_obj(M_p, J, eps_xi_eta, eps_xy, obj.n_C_1, obj.next_curve.n_C_0, d, nan, nan);
                [X_L, Y_L] = C_patch.L.xy_mesh;
                C_patch.L.f_XY = f(X_L, Y_L);
                [X_W, Y_W] = C_patch.W.xy_mesh;
                C_patch.W.f_XY = f(X_W, Y_W);
            else
                xi_diff = 2*(obj.n_C_1-1)*obj.h_tan;
                eta_diff = -2*(obj.next_curve.n_C_0-1)*obj.next_curve.h_tan;
                xi_0 = 1-(obj.n_C_1-1)*obj.h_tan;
                eta_0 = (obj.next_curve.n_C_0-1)*obj.next_curve.h_tan;
                
                xi_tilde = @(xi) xi_diff*xi + xi_0;
                eta_tilde = @(eta) eta_diff*eta +eta_0;
                
                M_p_1 = @(xi, eta) obj.l_1(xi_tilde(xi)) + obj.next_curve.l_1(eta_tilde(eta)) - obj.l_1(1);
                M_p_2 = @(xi, eta) obj.l_2(xi_tilde(xi)) + obj.next_curve.l_2(eta_tilde(eta)) - obj.l_2(1);
                M_p = @(xi, eta) [M_p_1(xi, eta), M_p_2(xi, eta)];
                
                dM_p_1_dxi = @(xi, eta) xi_diff*obj.l_1_prime(xi_tilde(xi));
                dM_p_2_dxi = @(xi, eta) xi_diff*obj.l_2_prime(xi_tilde(xi));
                dM_p_1_deta = @(xi, eta) eta_diff*obj.next_curve.l_1_prime(eta_tilde(eta));
                dM_p_2_deta = @(xi, eta) eta_diff*obj.next_curve.l_2_prime(eta_tilde(eta));
                
                J = @(v) [dM_p_1_dxi(v(1), v(2)), dM_p_1_deta(v(1), v(2)); dM_p_2_dxi(v(1), v(2)), dM_p_2_deta(v(1), v(2))];

                C_patch = C1_patch_obj(M_p, J, eps_xi_eta, eps_xy, obj.n_C_1*2-1, obj.next_curve.n_C_0*2-1, d, nan, nan);
                [X_L, Y_L] = C_patch.L.xy_mesh;
                C_patch.L.f_XY = f(X_L, Y_L);
                [X_W, Y_W] = C_patch.W.xy_mesh;
                C_patch.W.f_XY = f(X_W, Y_W);
            end
        end
        
        function curve_length = compute_length(obj)
            curve_length = integral(@(theta) sqrt(obj.l_1_prime(theta).^2+obj.l_2_prime(theta).^2), 0, 1);
        end
    end
end

