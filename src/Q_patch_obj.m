% Q_PATCH_OBJ - Generic patch objects used by S-type, C1-type, and C2-type
% patches. Assumes a rectangular domain in parameter space.
%
% Properties:
%   M_p - parametrization that maps Q from parameter space to real space,
%   assumed to be vectorized
%   J - Jacobian of M_p - 
%   n_xi - number of points xi axis is discretized into
%   n_eta - number of points eta axis is discretized into
%   xi_start - minimum value of xi in Q
%   xi_end - maximum value of xi in Q
%   eta_start - minimum value of eta in Q
%   eta_end - maximum value of eta in Q
%   f_XY - function values associated with patch on discretized mesh
%   x_min - minimum x value of M_p(Q)
%   x_max - maximum x value of M_p(Q)
%   y_min - minimum y alue of M_p(Q)
%   y_max - maximum y value of M_p(Q)
%   w_1D - one dimensional window function used to construct partition of
%       unity
%   w - unnormalized partition of unity function for this patch
%   eps_xi_eta - error tolerance in xi-eta space
%   eps_xy - error tolerance in x-y space, should be dependent on
%       eps_xi_eta and the maximum value of the J in Q
%
% Methods:
%   Q_patch_obj - Class constructor.
%   h_mesh - returns meshsize for xi and eta for object
%   xi_mesh - returns discretized xi mesh for object
%   eta_mesh - returns discretized eta mesh for object
%   xi_eta_mesh - returns entire discretzied (xi, eta) mesh associated with
%       patch -- i.e. the domain of the patch
%   xy_mesh - returns M_p(xi_eta_mesh)
%   boundary_mesh - returns discretized mesh of boundary in parameter space
%   boundary_mesh_xy - returns M_p(boundary_mesh_xy)
%   convert_to_XY - converts M_p(xi, eta) for given vector/matrices M_p,
%   in_patch - returns whether a given mesh in xi-eta space is in the domain
%   round_boundary_points - rounds points within the prescribed error tolerance near the boundary to exact
%       boundary points
%   inverse_M_p - computes the inverse of M_p numerically using Newton's
%       method
%   locally_compute - computes the function value of some xi-eta point not
%       in the domain of the patch using polynomial interpolation
%   compute_w_normalization_xi_right - computes partition of unity
%       normalization values for a "main" Q patch and "window" Q patch where
%       the window Q patch is to the right of the main patch
%   compute_w_normalization_xi_left - computes partition of unity
%       normalization values where window patch is to the left of the
%       main patch
%   compute_w_normalization_eta_up - computes partition of unity
%       normalization values where window patch is above the main patch
%   compute_w_normalization_eta_down - computes partition of unity
%       normalization values where window patch is below the main patch
%
%
%   Note the compute_w_normalization functions operate on this
%   assumptions about the window patches:
%       - the window patch is parametrized such that the "bounding" edge of
%       the window patch (edge of window patch that is"contained" within
%       main patch) corresponds to eta-axis in the window patch's parameter
%       space
%
% Author: Allen Yang
% Email: aryang@caltech.edu

classdef Q_patch_obj < handle
    properties
        M_p
        J
        n_xi
        n_eta
        xi_start
        xi_end
        eta_start
        eta_end
        f_XY
        x_min
        x_max
        y_min
        y_max
        w_1D
        
        eps_xi_eta
        eps_xy
    end
    
    methods
        function obj = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY)
            % Q_PATCH_OBJ Constructor for the class.
            %    obj = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi,
            %    n_eta, xi_start, xi_end, eta_start, eta_end, f_XY) intializes the object with the given
            %    properties
            %    
            % Input parameters:
            %    M_p (function handle) - parametrization from xi-eta space
            %       to x-y space. Takes in a column vector of xi values and
            %       column vector of eta values and returns a 2xlength(xi) =
            %       2xlength(eta) matrix of x-y values
            %    J (function handle) - Jacobian of M_p. Takes in a 2x1
            %       vector [xi; eta] and returns the 2x2 Jacobian for that
            %       given xi-eta
            %    eps_xi_eta (double) - error tolerance in xi-eta space
            %    eps_xy (double) - error tolerance in x-y space, is
            %       dependent on eps_xi_eta and J
            %    n_xi (int) - number of values xi axis is discretized into
            %    n_eta (int) - number of values eta axis is discretized into
            %    xi_start (double) - minimum xi value of rectangle in parameter space 
            %    xi_end (double) - maximum xi value of rectangle in
            %       parameter space
            %    eta_start (double) - minimum eta value of rectangle in
            %       paramter space
            %    eta_end (double) -maximum eta value of rectangle in
            %       parameter space
            %    f_XY (double) - n_eta x n_xi sized matrix containing
            %       function values associated with patch
            
            obj.M_p = M_p;
            obj.J = J;
            
            obj.eps_xi_eta = eps_xi_eta;
            obj.eps_xy = eps_xy;
            
            obj.n_eta = n_eta;
            obj.n_xi = n_xi;
            obj.xi_start = xi_start;
            obj.xi_end = xi_end;
            obj.eta_start = eta_start;
            obj.eta_end = eta_end;
            obj.f_XY = f_XY;
            
            % sets XY bounds
            [XI, ETA] = obj.xi_eta_mesh();
            XY = M_p(XI(:), ETA(:));
            obj.x_min = min(XY(:, 1));
            obj.x_max = max(XY(:, 1));
            obj.y_min = min(XY(:, 2));
            obj.y_max = max(XY(:, 2));
            
            obj.w_1D = @(x) erfc(6*(-2*x+1))/2;
        end

        function [h_xi, h_eta] = h_mesh(obj)
            % h_mesh returns the step sizes for xi and eta
            h_xi =(obj.xi_end-obj.xi_start)./ (obj.n_xi-1);
            h_eta = (obj.eta_end-obj.eta_start) ./ (obj.n_eta-1);
        end
        
        function mesh = xi_mesh(obj)
            % xi_mesh returns the discretized xi mesh of patch
            mesh = transpose(linspace(obj.xi_start, obj.xi_end, obj.n_xi));
        end
        
        function mesh = eta_mesh(obj)
            % eta_mesh returns the discretized eta mesh of patch
            mesh = transpose(linspace(obj.eta_start, obj.eta_end, obj.n_eta));
        end
        
        function [XI, ETA] = xi_eta_mesh(obj)
            % xi_eta_mesh returns entire xi-eta cartesian mesh of patch
            [XI, ETA] = meshgrid(obj.xi_mesh(), obj.eta_mesh());
        end
        
        function [X, Y] = xy_mesh(obj)
            % xy_mesh returns entire x-y cartesian mesh of patch
            %    note this is equivalent to M_p(XI, ETA)
            [XI, ETA] = obj.xi_eta_mesh();
            [X, Y] = obj.convert_to_XY(XI, ETA);
        end
        
        function [boundary_mesh_xi, boundary_mesh_eta] = boundary_mesh(obj, n_r, pad_boundary)
            % boundary_mesh returns a discrete mesh of the boundary of patch in xi-eta
            %   space
            % 
            % Input parameters:
            %    pad_boundary (boolean) - pads boundary mesh with a little
            %       space if true, used to compensate for points lost
            %       possibly lost by discretized boundary
            % TODO: let number of points being used in discretization be
            %   input parameter
            if pad_boundary
                [h_xi, h_eta] = obj.h_mesh;
            else
                h_xi = 0;
                h_eta = 0;
            end
            
            n_r_xi_mesh = transpose(linspace(obj.xi_start, obj.xi_end, obj.n_xi*n_r));
            n_r_eta_mesh = transpose(linspace(obj.eta_start, obj.eta_end, obj.n_eta*n_r));
            boundary_mesh_xi = [ones(obj.n_eta*n_r, 1)*(obj.xi_start-h_xi); n_r_xi_mesh; ones(obj.n_eta*n_r, 1)*(obj.xi_end+h_xi); flip(n_r_xi_mesh); obj.xi_start-h_xi];
            boundary_mesh_eta = [n_r_eta_mesh; ones(obj.n_xi*n_r, 1)*(obj.eta_end+h_eta); flip(n_r_eta_mesh); ones(obj.n_xi*n_r, 1)*(obj.eta_start-h_eta); obj.eta_start-h_eta];
        end
        
        function [boundary_mesh_x, boundary_mesh_y] = boundary_mesh_xy(obj, n_r, pad_boundary)
            % boundary_mesh_xy returns the x-y coordinates of the
            %   discretized boundary mesh from boundary_mesh
            [boundary_mesh_xi, boundary_mesh_eta] = obj.boundary_mesh(n_r, pad_boundary);
            [boundary_mesh_x, boundary_mesh_y] = obj.convert_to_XY(boundary_mesh_xi, boundary_mesh_eta);
        end
        
        function [X, Y] = convert_to_XY(obj, XI, ETA)
            % convert_to_XY converts a certain XI, ETA mesh to X, Y
            %   coordinates, XI and ETA are assumed to be the same size
            % 
            % Input parameters:
            %   XI (double) - matrix of xi values
            %   ETA (double) - matrix of eta values
            %
            % Returns:
            %   X (double) - matrix of converted X values returned in the
            %       same size as XI and ETA
            %   Y (double) - matrix of converted Y values returned in the
            %       same size as XI and ETA
            XY = obj.M_p(XI(:), ETA(:));
            X = reshape(XY(:, 1), size(XI));
            Y = reshape(XY(:, 2), size(ETA));
        end
        
        function patch_msk = in_patch(obj, xi, eta)
            % in_patch returns a boolean mask for given xi, eta matrices
            %   determining wheter the coordinates are within the bounds of
            %   the patch
            patch_msk = xi >= obj.xi_start & xi <= obj.xi_end & eta >= obj.eta_start & eta <= obj.eta_end;
        end
        
        function [xi, eta] = round_boundary_points(obj, xi, eta)
            % round_boundary_points takes given xi, eta and rounds points
            %   that are within the xi, eta error tolerance to exact boundary values 
            xi(abs(xi - obj.xi_start) < obj.eps_xi_eta) = obj.xi_start;
            xi(abs(xi-obj.xi_end) < obj.eps_xi_eta) = obj.xi_end;
            eta(abs(eta - obj.eta_start) < obj.eps_xi_eta) = obj.eta_start;
            eta(abs(eta-obj.eta_end) < obj.eps_xi_eta) = obj.eta_end;
        end
        
        function [xi, eta, converged] = inverse_M_p(obj, x, y, initial_guesses)
            % inverse_M_p uses Newton's method to compute the inverse of
            %   M_p for scalar x and y
            %
            % Input parameters:
            %    x (double) - scalar x value
            %    y (double) - scalar y value
            %    initial_guesses (double) - 2 x n matrix of n initial
            %       guesses to use for Newton's method. If nan, the method
            %       will use a boundary mesh as the initial guesses
            
            if isnan(initial_guesses)
                initial_guesses = obj.default_initial_guesses(20);      
            end

            err_guess = @(x, y, v) transpose(obj.M_p(v(1), v(2))) - [x; y];

            for k = 1:size(initial_guesses, 2)
                initial_guess = initial_guesses(:, k);
                [v_guess, converged] = newton_solve(@(v) err_guess(x, y, v), obj.J, initial_guess, obj.eps_xy, 1000);
                
                [xi, eta] = obj.round_boundary_points(v_guess(1), v_guess(2));                                                
                if converged && obj.in_patch(xi, eta)
                    return
                end
            end
        end
        
        function initial_guesses = default_initial_guesses(obj, N)
            N_segment = ceil(N/4);
                
            xi_mesh = transpose(linspace(obj.xi_start, obj.xi_end, N_segment+1));
            eta_mesh = transpose(linspace(obj.eta_start, obj.eta_end, N_segment+1));

            xi_initial = [xi_mesh(1:end-1); xi_mesh(1:end-1); ones(N_segment, 1)*obj.xi_start; ones(N_segment, 1)*obj.xi_end];
            eta_initial = [ones(N_segment, 1)*obj.eta_start; ones(N_segment, 1)*obj.eta_end; eta_mesh(1:end-1); eta_mesh(1:end-1)];

            initial_guesses = [xi_initial'; eta_initial'];
        end

        function [f_xy, in_range] = locally_compute(obj, xi, eta, M)
            % locally_compute uses two step one dimensional polynomial interpolation to estimate the
            %   value of f for any xi, eta within the bounds of the patch (not necessarily on discrete mesh)
            %
            % Input parameters:
            %    xi (double) - scalar xi value within bounds of patch in
            %       parameter space
            %    eta (double) - scalar eta value within bounds of patch in
            %       parameter space
            %   M (int) - number of interpolation points used for each one
            %       dimensional procedure
            
            % check if xi, eta are within bounds of Q
            if ~obj.in_patch(xi, eta)
                f_xy = nan;
                in_range = false;
                warning('locally computing for point not in patch');
                return
            end
            
            in_range = true;
            [h_xi, h_eta] = obj.h_mesh();

            % j to the immediate left of point
            xi_j = floor((xi-obj.xi_start)/h_xi);
            eta_j = floor((eta-obj.eta_start)/h_eta);
            
            half_M =  floor(M/2);
            if mod(M, 2) ~= 0
                interpol_xi_j_mesh = transpose(xi_j-half_M:xi_j+half_M);
                interpol_eta_j_mesh = transpose(eta_j-half_M:eta_j+half_M);
            else
                interpol_xi_j_mesh = transpose(xi_j-half_M+1:xi_j+half_M);
                interpol_eta_j_mesh = transpose(eta_j-half_M+1:eta_j+half_M);
            end

            interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 0, obj.n_xi-1);
            interpol_eta_j_mesh = shift_idx_mesh(interpol_eta_j_mesh, 0, obj.n_eta-1);
            
            interpol_xi_mesh = h_xi*interpol_xi_j_mesh + obj.xi_start;
            interpol_eta_mesh = h_eta*interpol_eta_j_mesh + obj.eta_start;
            
            % first 1D interpolation
            interpol_xi_exact = zeros(M, 1);
            for horz_idx = 1:M
                interpol_val = obj.f_XY(interpol_eta_j_mesh(horz_idx)+1, interpol_xi_j_mesh+1)';
                interpol_xi_exact(horz_idx) = barylag([interpol_xi_mesh, interpol_val], xi);
            end
             % second 1D interpolation
            f_xy = barylag([interpol_eta_mesh, interpol_xi_exact], eta);
        end
        
        function apply_w_normalization_xi_right(obj, window_patch)
            % compute_w_normalization_xi_right computes the partition of
            %   unity normalization values for when the window patch is to
            %   the right of the patch
            %
            % Input parameter:
            %   window_patch (Q_patch_obj) - window patch
            % 
            % Returns:
            %   main_patch_w (double) - normalization constants for the
            %       domain of main patch, has the same shape as obj.f_XY
            %   window_patch_w (double) - normalization constants for the
            %       domain of window patch, has the same shape as
            %       window_patch.f_XY
            
            main_xi_corner = compute_xi_corner(obj, window_patch, true, window_patch.xi_end, true);
            window_xi_corner = compute_xi_corner(window_patch, obj, true, obj.xi_end, false);

            main_R_xi = main_xi_corner - obj.xi_end;
            main_xi_0 = obj.xi_end;
            main_w = @(xi, eta) obj.w_1D((xi-main_xi_0)/main_R_xi);
            
            window_w = compute_w_normalization_window(obj, main_w, window_patch, window_xi_corner, false);

            [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_xi_overlap_mesh(obj, main_xi_corner, true);
            [X_overlap, Y_overlap] = obj.convert_to_XY(XI_overlap, ETA_overlap);    
            
            w_unnormalized =  main_w(XI_overlap, ETA_overlap);
            initial_guesses = [linspace(window_xi_corner, window_patch.xi_end, 20); window_patch.eta_start * ones(1, 20)];
            obj.apply_w(w_unnormalized, window_patch, true, window_w, X_overlap, Y_overlap, XI_j, ETA_j, initial_guesses);
        end
        
        function apply_w_normalization_xi_left(obj, window_patch)
            % compute_w_normalization_xi_left computes the partition of
            %   unity normalization values for when the window patch is to
            %   the left of the patch
            %
            % Input parameter:
            %   window_patch (Q_patch_obj) - window patch
            % 
            % Returns:
            %   main_patch_w (double) - normalization constants for the
            %       domain of main patch, has the same shape as obj.f_XY
            %   window_patch_w (double) - normalization constants for the
            %       domain of window patch, has the same shape as
            %       window_patch.f_XY
            main_xi_corner = compute_xi_corner(obj, window_patch, true, window_patch.xi_end, false);
            window_xi_corner = compute_xi_corner(window_patch, obj, true, obj.xi_start, false);

            main_R_xi = main_xi_corner - obj.xi_start;
            main_xi_0 = obj.xi_start;
            main_w = @(xi, eta) obj.w_1D((xi-main_xi_0)/main_R_xi);
            
            window_w = compute_w_normalization_window(obj, main_w, window_patch, window_xi_corner, false);

            [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_xi_overlap_mesh(obj, main_xi_corner, false);
            [X_overlap, Y_overlap] = obj.convert_to_XY(XI_overlap, ETA_overlap);                
            
            w_unnormalized = main_w(XI_overlap, ETA_overlap);
            initial_guesses = [linspace(window_xi_corner, window_patch.xi_end, 20); window_patch.eta_start * ones(1, 20)];
            obj.apply_w(w_unnormalized, window_patch, true, window_w, X_overlap, Y_overlap, XI_j, ETA_j, initial_guesses);        
        end
    
        function apply_w_normalization_eta_up(obj, window_patch)
            % compute_w_normalization_eta_up computes the partition of
            %   unity normalization values for when the window patch above
            %   the patch
            %
            % Input parameter:
            %   window_patch (Q_patch_obj) - window patch
            % 
            % Returns:
            %   main_patch_w (double) - normalization constants for the
            %       domain of main patch, has the same shape as obj.f_XY
            %   window_patch_w (double) - normalization constants for the
            %       domain of window patch, has the same shape as
            %       window_patch.f_XY
            main_eta_corner = compute_eta_corner(obj, window_patch, true, window_patch.xi_start, true);
            window_xi_corner = compute_xi_corner(window_patch, obj, false, obj.eta_end, false);

            main_R_eta = main_eta_corner - obj.eta_end;
            main_eta_0 = obj.eta_end;
            main_w = @(xi, eta) obj.w_1D((eta-main_eta_0)/main_R_eta);
            
            window_w = compute_w_normalization_window(obj, main_w, window_patch, window_xi_corner, true);

            [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_eta_overlap_mesh(obj, main_eta_corner, true);
            [X_overlap, Y_overlap] = obj.convert_to_XY(XI_overlap, ETA_overlap);         
            
            w_unnormalized = main_w(XI_overlap, ETA_overlap);
            initial_guesses = [linspace(window_patch.xi_start, window_xi_corner, 20);window_patch.eta_start * ones(1, 20)];
            obj.apply_w(w_unnormalized, window_patch, true, window_w, X_overlap, Y_overlap, XI_j, ETA_j, initial_guesses);
        end
        
        function apply_w_normalization_eta_down(obj, window_patch)
            % compute_w_normalization_eta_down computes the partition of
            %   unity normalization values for when the window patch is
            %   below the patch
            %
            % Input parameter:
            %   window_patch (Q_patch_obj) - window patch
            % 
            % Returns:
            %   main_patch_w (double) - normalization constants for the
            %       domain of main patch, has the same shape as obj.f_XY
            %   window_patch_w (double) - normalization constants for the
            %       domain of window patch, has the same shape as
            %       window_patch.f_XY
            main_eta_corner = compute_eta_corner(obj, window_patch, true, window_patch.xi_start, false);
            window_xi_corner = compute_xi_corner(window_patch, obj, false, obj.eta_start, false);

            main_R_eta = main_eta_corner - obj.eta_start;
            main_eta_0 = obj.eta_start;
            main_w = @(xi, eta) obj.w_1D((eta-main_eta_0)/main_R_eta);
            
            [window_w] = compute_w_normalization_window(obj, main_w, window_patch, window_xi_corner, true);

            [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_eta_overlap_mesh(obj, main_eta_corner, false);
            [X_overlap, Y_overlap] = obj.convert_to_XY(XI_overlap, ETA_overlap);        
            
            w_unnormalized = main_w(XI_overlap, ETA_overlap);
            initial_guesses = [linspace(window_patch.xi_start, window_xi_corner, 20);window_patch.eta_start * ones(1, 20)];
            obj.apply_w(w_unnormalized, window_patch, true, window_w, X_overlap, Y_overlap, XI_j, ETA_j, initial_guesses);
        end
        
        function apply_w(obj, w_unnormalized, window_patch, window_patch_bound_xi, window_w, overlap_X, overlap_Y, overlap_XI_j, overlap_ETA_j, initial_guesses)            
            for i = 1:size(overlap_X, 1)
                if mod(i, 2) == 1
                    j_lst = 1:size(overlap_X, 2);
                else
                    j_lst = size(overlap_X, 2):-1:1;
                end

                for j = j_lst
                    if isnan(initial_guesses)
                        initial_guesses = obj.default_initial_guesses(20);
                    end
                    
                    for initial_guess = initial_guesses
                        [window_patch_xi, window_patch_eta, converged] = window_patch.inverse_M_p(overlap_X(i, j), overlap_Y(i, j), initial_guesses);
                        if i == 1 && j ==1
                            in_VpR = window_patch.in_patch(window_patch_xi, window_patch_eta);
                            if ~in_VpR
                                warning("First Point when computing overlaps for POU should be in patch, something may be wrong with geometry");
                            end
                        else
                            if window_patch_bound_xi
                                in_VpR = window_patch_xi >= window_patch.xi_start && window_patch_xi <= window_patch.xi_end;
                            else
                                in_VpR = window_patch_eta >= window_patch.eta_start && window_patch_eta <= window_patch.eta_end;
                            end
                        end

                        if converged && in_VpR
                            break;
                        end
                    end

                    if converged && in_VpR
                        xi_j = overlap_XI_j(i, j) + 1;
                        eta_j = overlap_ETA_j(i, j) + 1;
                        obj.f_XY(eta_j, xi_j) = obj.f_XY(eta_j, xi_j) .* w_unnormalized(i, j) ./ (window_w(window_patch_xi, window_patch_eta)+w_unnormalized(i, j));
                    elseif ~converged
                        warning("Nonconvergence in computing C2_norm")
                    end

                    if converged && in_VpR
                        initial_guesses = [window_patch_xi; window_patch_eta]; % using previous point as next initial guess
                    end
                end
            end
        end
    end
end

function [main_xi_corner] = compute_xi_corner(main_patch, window_patch, window_fix_xi, window_fixed_edge, window_patch_right)
%COMPUTE_XI_CORNER Summary of this function goes here
%   Detailed explanation goes here
        [h_xi_window, h_eta_window] = window_patch.h_mesh;
        
        if window_fix_xi
            window_xi_edge = window_fixed_edge;
            window_eta_edge = window_patch.eta_start;
        else
            window_xi_edge = window_patch.xi_start;
            window_eta_edge = window_fixed_edge;
        end
                
        if window_patch_right
            main_xi_corner = main_patch.xi_end;
        else
            main_xi_corner = main_patch.xi_start;
        end
        
        first_iter = true;
        while true
            window_xy_edge = window_patch.M_p(window_xi_edge, window_eta_edge);
            if first_iter
                [main_xi_edge, main_eta_edge, converged] = main_patch.inverse_M_p(window_xy_edge(1, 1), window_xy_edge(1, 2), nan);
            else
                [main_xi_edge, main_eta_edge, converged] = main_patch.inverse_M_p(window_xy_edge(1, 1), window_xy_edge(1, 2), [main_xi_edge; main_eta_edge]);
            end
            first_iter = false;
            
            if ~converged
                 warning("Nonconvergence in computing boundary mesh values");
                 break;
            end
            if (main_xi_edge < main_xi_corner && window_patch_right) || (main_xi_edge > main_xi_corner && ~window_patch_right)
                main_xi_corner = main_xi_edge;
            end
            if main_eta_edge > main_patch.eta_end || main_eta_edge < main_patch.eta_start
                break;
            end
            if window_fix_xi
                window_eta_edge = window_eta_edge + h_eta_window;
            else
                window_xi_edge = window_xi_edge + h_xi_window;
            end
        end
end

function [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_xi_overlap_mesh(main_patch, xi_corner, window_patch_right)
%COMPUTE_XI_OVERLAP_MESH Summary of this function goes here
%   Detailed explanation goes here
    [h_xi, h_eta] = main_patch.h_mesh;
    
    if window_patch_right
        xi_corner_j = ceil((xi_corner - main_patch.xi_start) / h_xi);
        [XI_j, ETA_j] = meshgrid((main_patch.n_xi-1):-1:xi_corner_j, 0:(main_patch.n_eta-1));        
    else
        xi_corner_j = floor((xi_corner - main_patch.xi_start) / h_xi);
        [XI_j, ETA_j] = meshgrid(0:xi_corner_j, 0:(main_patch.n_eta-1));
    end
    
    XI_overlap = XI_j * h_xi + main_patch.xi_start;
    ETA_overlap = ETA_j * h_eta + main_patch.eta_start;
end

function [main_eta_corner] = compute_eta_corner(main_patch, window_patch, window_fix_xi, window_fixed_edge, window_patch_up)
%COMPUTE_XI_CORNER Summary of this function goes here
%   Detailed explanation goes here
        [h_xi_window, h_eta_window] = window_patch.h_mesh;
        
        if window_fix_xi
            window_xi_edge = window_fixed_edge;
            window_eta_edge = window_patch.eta_start;
        else
            window_xi_edge = window_patch.xi_start;
            window_eta_edge = window_fixed_edge;
        end
                
        if window_patch_up
            main_eta_corner = main_patch.eta_end;
        else
            main_eta_corner = main_patch.eta_start;
        end
        
        first_iter = true;
        while true
            window_xy_edge = window_patch.M_p(window_xi_edge, window_eta_edge);
            if first_iter
                [main_xi_edge, main_eta_edge, converged] = main_patch.inverse_M_p(window_xy_edge(1, 1), window_xy_edge(1, 2), nan);
            else
                [main_xi_edge, main_eta_edge, converged] = main_patch.inverse_M_p(window_xy_edge(1, 1), window_xy_edge(1, 2), [main_xi_edge; main_eta_edge]);
            end
            
            first_iter = false;
            if ~converged
                 warning("Nonconvergence in computing boundary mesh values");
                 break;
            end
            if (main_eta_edge < main_eta_corner && window_patch_up) || (main_eta_edge > main_eta_corner && ~window_patch_up)
                main_eta_corner = main_eta_edge;
            end
            if main_xi_edge > main_patch.xi_end || main_xi_edge < main_patch.xi_start
                break;
            end
            if window_fix_xi
                window_eta_edge = window_eta_edge + h_eta_window;
            else
                window_xi_edge = window_xi_edge + h_xi_window;
            end
        end
end

function [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_eta_overlap_mesh(main_patch, eta_corner, window_patch_up)
%COMPUTE_XI_OVERLAP_MESH Summary of this function goes here
%   Detailed explanation goes here
    [h_xi, h_eta] = main_patch.h_mesh;
    
    if window_patch_up
        eta_corner_j = ceil((eta_corner - main_patch.eta_start) / h_eta);
        [XI_j, ETA_j] = meshgrid(0:(main_patch.n_xi-1), (main_patch.n_eta-1):-1:eta_corner_j);        
    else
        eta_corner_j = floor((eta_corner - main_patch.eta_start) / h_eta);
        [XI_j, ETA_j] = meshgrid(0:(main_patch.n_xi-1), 0:eta_corner_j);
    end
    
    XI_overlap = XI_j * h_xi + main_patch.xi_start;
    ETA_overlap = ETA_j * h_eta + main_patch.eta_start;
end

function [window_w] = compute_w_normalization_window(main_patch, main_w, window_patch, window_xi_corner, up_down)
    if up_down
        window_R_xi = window_xi_corner - window_patch.xi_start;
        window_xi_0 = window_patch.xi_start;
        window_w = @(xi, eta) window_patch.w_1D((xi-window_xi_0)/window_R_xi);
    else
        window_R_xi = window_xi_corner - window_patch.xi_end;
        window_xi_0 = window_patch.xi_end;
        window_w = @(xi, eta) window_patch.w_1D((xi-window_xi_0)/window_R_xi);
    end

    [XI_overlap_window, ETA_overlap_window, XI_j_window, ETA_j_window] = compute_xi_overlap_mesh(window_patch, window_xi_corner, ~up_down);
    [X_overlap_window, Y_overlap_window] = window_patch.convert_to_XY(XI_overlap_window, ETA_overlap_window);
    
    w_unnormalized = window_w(XI_overlap_window, ETA_overlap_window);
    window_patch.apply_w(w_unnormalized, main_patch, ~up_down, main_w, X_overlap_window, Y_overlap_window, XI_j_window, ETA_j_window, nan);
end
