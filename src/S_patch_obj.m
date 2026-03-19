classdef S_patch_obj < handle
    %S_PATCH_OBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Q
        h
    end
    
    methods
        function obj = S_patch_obj(M_p, J, h, eps_xi_eta, eps_xy, n_xi, d, f_XY)
            %S_PATCH_OBJ Construct an instance of this class
            %   Detailed explanation goes here            
            obj.Q = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, d, 0, 1, 0, h*(d-1), f_XY);
            
            obj.h = h;
        end
        
        function S_patch_copy = copy_patch(obj)
            S_patch_copy = S_patch_obj(obj.Q.M_p, obj.Q.J, obj.h, obj.Q.eps_xi_eta, obj.Q.eps_xy, obj.Q.n_xi, obj.Q.n_eta, obj.Q.f_XY);
        end
        
        function S_fcont_patch = FC(obj, C, n_r, d, A, Q)
            h_eta = (obj.Q.eta_end-obj.Q.eta_start)/(obj.Q.n_eta-1);
            fcont = fcont_gram_blend_S(obj.Q.f_XY, d, A, Q);
            
            S_fcont_patch = Q_patch_obj(obj.Q.M_p, obj.Q.J, obj.Q.eps_xi_eta, obj.Q.eps_xy, obj.Q.n_xi, C*n_r+1, obj.Q.xi_start, obj.Q.xi_end, obj.Q.eta_start-C*h_eta, obj.Q.eta_start, fcont);
        end
    end
end

