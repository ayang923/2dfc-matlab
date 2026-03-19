% C2_PATCH_OBJ - Wrapper function for C2_patch type objects
%
% Properties:
%   L - "Long" Q_patch_obj corresponding to the region [0, (d-1)h] x [(d-1)h, 1].
%   W - "Wide" Q_patch_obj corresponding to the region [0, 1] x [0, (d-1h)].
%
% Methods:
%   C2_patch_obj - Class constructor.
%   FC - Computes extension values for the C2 patch given.
%   compute_w_normalization - Computes normalization factor for the
%   partition of unity associated with C2 patch and two window patches
%
% Author: Allen Yang
% Email: aryang@caltech.edu

classdef C2_patch_obj < handle
    properties
        L Q_patch_obj
        W Q_patch_obj
    end
    
    methods
        function obj = C2_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, d, f_L, f_W)
            % C2_PATCH_OBJ Constructor for the class.
            %    obj = C2_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi,
            %    n_eta, d, f_L, f_W) intializes the object with the given
            %    properties
            
            h_xi = 1./(n_xi-1);
            h_eta = 1./(n_eta-1);
            obj.W = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, d, 0, 1, 0, (d-1)*h_eta, f_W);            
            obj.L = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, d, n_eta-d+1, 0, (d-1)*h_xi, (d-1)*h_eta, 1, f_L);
        end
        
        function C2_patch_copy = copy_patch(obj)
            d = obj.L.n_xi;
            n_xi = obj.W.n_xi;
            n_eta = obj.L.n_eta+d-1;
            C2_patch_copy = C2_patch_obj(obj.L.M_p, obj.L.J, obj.L.eps_xi_eta, obj.L.eps_xy, n_xi, n_eta, d, obj.L.f_XY, obj.W.f_XY);
        end
                
        function [C2_fcont_patch_L, C2_fcont_patch_W, C2_fcont_patch_corner] = FC(obj, C, n_r, d, A, Q, M)
            % FC Computes the blending-to-zero extension values for this
            % patch and returns them as three Q-patch type objects
            %
            % Input Parameters:
            %    C (int): number of unrefined extension values generated
            %       in one dimensional blending-to-zero procedure
            %	 n_r (int): refinement factor of one dimensional
            %       blending-to-zero procedure
            %    d (int): number of "matching points" used for one
            %       dimensional blending-to-zero procedure
            %    A (double): precomputed A (n_r*C x d) matrix in one
            %       dimensional blending-to-zero procedure
            %    Q (double): precomputed Q (d x d) matrix in one dimensional
            %       blending-to-zero procedure
            %    L_norm (double): partition of unity
            %       normalization values for L
            %    W_norm (double): partition of unity
            %       normalization values for W
            % Output:
            %	 C1_fcont_patch_L (Q_patch_obj): blending-to-zero values
            %       for L
            %    C1_fcont_patch_W (Q_patch_obj): blending-to-zero
            %       values for W
            %    C1_fcont_patch_corner (Q_patch_obj): blending-to-zero
            %       values in corner region
            
            [h_xi, h_eta] = obj.W.h_mesh; %h mesh is the same for both W and L
            f_XY_W = obj.W.f_XY;
            f_XY_L = obj.L.f_XY;
            
            fcont_W = fcont_gram_blend_S(f_XY_W, d, A, Q);
            fcont_corner = transpose(fcont_gram_blend_S(fcont_W', d, A, Q));
            fcont_L = transpose(fcont_gram_blend_S([f_XY_W(:, 1:d); f_XY_L(2:end, :)]', d, A, Q));
            
            
            C2_fcont_patch_W =  Q_patch_obj(obj.W.M_p, obj.W.J, obj.W.eps_xi_eta, obj.W.eps_xy, obj.W.n_xi, C*n_r+1, 0, 1, -(C)*h_eta, 0, fcont_W);
            C2_fcont_patch_L = Q_patch_obj(obj.L.M_p, obj.L.J, obj.L.eps_xi_eta, obj.L.eps_xy, C*n_r+1, obj.L.n_eta + d-1, -(C)*h_xi, 0, 0, 1, fcont_L);
            C2_fcont_patch_corner = Q_patch_obj(obj.W.M_p, obj.W.J, obj.W.eps_xi_eta, obj.W.eps_xy, C*n_r+1, C*n_r+1, -(C)*h_xi, 0, -C*h_eta, 0, fcont_corner);
        end
        
        function apply_w_W(obj, window_patch_W)
            % compute_w_normalization computes normalization constants for
            % partition of unity associated with L, W, and their respective
            % window patches
            %
            % Input Parameters:
            %    window_patch_L (S_patch_obj): S-type patch that acts as a
            %       "window patch" for L
            %    window_patch_W (S_patch_obj): S-type patch that acts as
            %       "window patch" for W
            % Output:
            %	 C1_L_norm (Q_patch_obj): partition of unity normalization
            %       values for L
            %    C1_W_norm (Q_patch_obj): partition of unity normalization
            %       values for W
            %    window_L_norm (Q_patch_obj): partition of unity normalization
            %       values for window_patch_L
            %    window_W_norm (Q_patch_obj): partition of unity normalization
            %       values for window_patch_eta
            obj.W.apply_w_normalization_xi_right(window_patch_W.Q);
        end
        
        function apply_w_L(obj, window_patch_L)
            obj.L.apply_w_normalization_eta_up(window_patch_L.Q);
        end
    end
end





