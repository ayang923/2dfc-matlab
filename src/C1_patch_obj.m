% C1_PATCH_OBJ - Wrapper function for C1_patch type objects
%
% Properties:
%   L - "Long" Q_patch_obj corresponding to the region [1/2, 1/2+(d-1)h] x [0, 1/2+(d-1)h].
%   W - "Wide" Q_patch_obj corresponding to the region [0, 1/2] x [1/2, 1/2+(d-1)h].
%
% Methods:
%   C1_patch_obj - Class constructor.
%   FC - Computes extension values for the C1 patch given.
%   refine_W - Refines function values in W for FC.
%   compute_w_normalization - Computes normalization factor for the
%   partition of unity associated with C1 patch and two window patches
%
% Author: Allen Yang
% Email: aryang@caltech.edu

classdef C1_patch_obj < handle    
    properties
        L Q_patch_obj
        W Q_patch_obj
    end
    
    methods
        function obj = C1_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, d, f_L, f_W)
            % C1_PATCH_OBJ Constructor for the class.
            %    obj = C1_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi,
            %    n_eta, d, f_L, f_W) intializes the object with the given
            %    properties
            assert(mod(n_xi, 2) == 1 && mod(n_eta, 2) == 1, "n_xi and n_eta must be odd"); % n_xi and n_eta must be odd for parametrization and discretization to work
            
            h_xi = 1./(n_xi-1);
            h_eta = 1./(n_eta-1);
            
            obj.L = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, d, (n_eta+1)/2+(d-1), 1/2, 1/2+(d-1)*h_xi, 0, 1/2+(d-1)*h_eta, f_L);
            obj.W = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, (n_xi+1)/2, d, 0, 1/2, 1/2, 1/2+(d-1)*h_eta, f_W);
        end
        
        function C2_patch_copy = copy_patch(obj)
            d = obj.L.n_xi;
            n_xi = obj.W.n_xi*2-1;
            n_eta = (obj.L.n_eta-(d-1))*2-1;
            C2_patch_copy = C2_patch_obj(obj.L.M_p, obj.L.J, obj.L.eps_xi_eta, obj.L.eps_xy, n_xi, n_eta, d, obj.L.f_XY, obj.W.f_XY);
        end
        
        function  [C1_fcont_patch_L, C1_fcont_patch_W_refined, C1_fcont_patch_W_unrefined] = FC(obj, C, n_r, d, A, Q, M)
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
            %    M (int): number of points used in polynomial interpolation to refine
            %       the function values of W
            %    L_norm (double): partition of unity
            %       normalization values for L
            %    W_norm (double): partition of unity
            %       normalization values for W
            % Output:
            %	 C1_fcont_patch_L (Q_patch_obj): blending-to-zero values
            %       for L
            %    C1_fcont_patch_W_refined (Q_patch_obj): blending-to-zero
            %       values for W in region that overlaps with the
            %       blending-to-zero values for L
            %    C1_fcont_patch_W_unrefined(Q_patch_obj): blending-to-zero
            %       values for W in region that doesn't overlap with
            %       blending-to-zero values for L
            L_f_XY = obj.L.f_XY;            
            W_f_XY = obj.W.f_XY;
            
            [h_xi, h_eta] = obj.L.h_mesh;

            L_fcont = transpose(fcont_gram_blend_S(L_f_XY', d, A, Q));
            C1_fcont_patch_L = Q_patch_obj(obj.L.M_p, obj.L.J, obj.L.eps_xi_eta, obj.L.eps_xy, C*n_r+1,  obj.L.n_eta - (d-1), 1/2-C*h_xi, 1/2, 0, 1/2, L_fcont(1:(obj.L.n_eta - (d-1)), :));
            
            if 1/2-C*h_xi < 0
                error("Please refine mesh of C1 patch");
            else
                [W_unrefined_f_XY, W_refined_f_XY] = obj.refine_W(W_f_XY, C, n_r, M);
                W_minus_fcont = W_refined_f_XY - L_fcont(end-(d-1):end, :);
                
                W_fcont_refined = fcont_gram_blend_S(W_minus_fcont, d, A, Q);
                W_fcont_unrefined = fcont_gram_blend_S(W_unrefined_f_XY, d, A, Q);
                
                C1_fcont_patch_W_refined = Q_patch_obj(obj.W.M_p, obj.W.J, obj.W.eps_xi_eta, obj.W.eps_xy, C*n_r+1, C*n_r+1, 1/2-C*h_xi, 1/2, 1/2-C*h_eta, 1/2, W_fcont_refined);
                C1_fcont_patch_W_unrefined = Q_patch_obj(obj.W.M_p, obj.W.J, obj.W.eps_xi_eta, obj.W.eps_xy, obj.W.n_xi-C, C*n_r+1, 0, 1/2-C*h_xi, 1/2-C*h_eta, 1/2, W_fcont_unrefined);
            end
        end
        
        
        
        function [W_unrefined_f_XY, W_refined_f_XY] = refine_W(obj, W_f_XY, C, n_r, M)
            % refine_W refines W for use in FC using polynomial
            % interpolation
            %
            % Input Parameters:
            %    W_f_XY (double): function values associated with W
            %    C (int): number of unrefined extension values generated
            %       in one dimensional blending-to-zero procedure
            %	 n_r (int): refinement factor of one dimensional
            %       blending-to-zero procedure
            %    M (int): number of "matching points" used for one
            %       dimensional blending-to-zero procedure
            %    A (double): precomputed A (n_r*C x d) matrix in one
            %       dimensional blending-to-zero procedure
            %    Q (double): precomputed Q (d x d) matrix in one dimensional
            %       blending-to-zero procedure
            %    M (int): number of points used in polynomial interpolation to refine
            %       the function values of W
            % Output:
            %	 W_unrefined_f_XY (Q_patch_obj): function values of W that
            %       don't overlap with blending-to-zero extension values of L
            %    W_refined_f_XY (Q_patch_obj): function values of W that
            %       are defined on a refined mesh that do overlap with
            %       blending-to-zero extension values of L
            [h_xi, ~] = obj.W.h_mesh;
            
            W_unrefined_f_XY = W_f_XY(:, 1:(obj.W.n_xi-C));
            
            W_refined_xi_mesh = transpose((1/2 - C*h_xi):(h_xi/n_r):1/2);
                                    
            W_refined_f_XY = zeros(obj.W.n_eta, length(W_refined_xi_mesh));
            half_M = floor(M/2);
            for eta_j = 0:obj.W.n_eta-1
                for xi_j = obj.W.n_xi-C-1:(obj.W.n_xi-2)
                    if mod(M, 2) ~= 0
                        interpol_xi_j_mesh = transpose(xi_j-half_M:xi_j+half_M);
                    else
                        interpol_xi_j_mesh = transpose(xi_j-half_M+1:xi_j+half_M);
                    end

                    interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 0, obj.W.n_xi-1);
                    interpol_xi_mesh = h_xi*interpol_xi_j_mesh;
                    interpol_val = W_f_XY(eta_j+1, interpol_xi_j_mesh+1)';
                    
                    xi_eval_j = ((xi_j-(obj.W.n_xi-C-1))*n_r):((xi_j-(obj.W.n_xi-C-1)+1)*n_r);
                    
                    W_refined_f_XY(eta_j+1, xi_eval_j+1) = barylag([interpol_xi_mesh, interpol_val], W_refined_xi_mesh(xi_eval_j+1));
                end
            end
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
            obj.W.apply_w_normalization_xi_left(window_patch_W.Q);
        end
        
        function apply_w_L(obj, window_patch_L)
            obj.L.apply_w_normalization_eta_down(window_patch_L.Q);
        end
    end
end


