% R_CARTESIAN_MESH_OBJ - Generic patch objects used by S-type, C1-type, and C2-type
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
%
% Author: Allen Yang
% Email: aryang@caltech.edu

classdef R_cartesian_mesh_obj < handle
    %Cartesian mesh being interpolated onto
    
    properties
        x_start
        x_end
        y_start
        y_end
        h
        n_x
        n_y
        x_mesh
        y_mesh
        R_X
        R_Y
        R_idxs
        
        boundary_X
        boundary_Y
        in_interior
        interior_idxs
        f_R
        
        fc_coeffs
    end
    
    methods
        function obj = R_cartesian_mesh_obj(x_start,x_end, y_start, y_end, h, boundary_X, boundary_Y)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.x_start = x_start;
            obj.y_start = y_start;
            obj.x_end = ceil((x_end-x_start)/h)*h+x_start;
            obj.y_end = ceil((y_end-y_start)/h)*h+y_start;
            
            obj.h = h;
            
            obj.n_x = round((obj.x_end-obj.x_start)/h)+1;
            obj.n_y = round((obj.y_end-obj.y_start)/h)+1;  
            
            % reduces round off error
            obj.x_mesh = transpose(round((linspace(obj.x_start, obj.x_end, obj.n_x)-x_start)/h)*h+x_start);
            obj.y_mesh = transpose(round((linspace(obj.y_start, obj.y_end, obj.n_y)-y_start)/h)*h+y_start);
            
            [obj.R_X, obj.R_Y] = meshgrid(obj.x_mesh, obj.y_mesh);
            obj.R_idxs = reshape(1:numel(obj.R_X), size(obj.R_X));

            obj.boundary_X = boundary_X;
            obj.boundary_Y = boundary_Y;

            obj.in_interior = inpolygon_mesh(obj.R_X, obj.R_Y, boundary_X,  boundary_Y);
            
            obj.interior_idxs = obj.R_idxs(obj.in_interior);
            obj.f_R = zeros(obj.n_y, obj.n_x);
        end
        
        function [R_patch_idxs] = interpolate_patch(obj, patch, n_r, M)
            % constructs vector of idxs of points that are in both patch
            % and cartesian mesh
            [bound_X, bound_Y] = patch.boundary_mesh_xy(n_r, false);
            in_patch = inpolygon_mesh(obj.R_X, obj.R_Y, bound_X, bound_Y) & ~obj.in_interior;
            R_patch_idxs = obj.R_idxs(in_patch);
            
            [P_xi, P_eta] = R_xi_eta_inversion(obj, patch, in_patch);
                        
            f_R_patch = zeros(size(R_patch_idxs));
            for i = 1:length(R_patch_idxs)
                [interior_val, in_range] = patch.locally_compute(P_xi(i), P_eta(i), M);
                if in_range
                     f_R_patch(i) = interior_val;                    
                end
            end
            obj.f_R(R_patch_idxs) = obj.f_R(R_patch_idxs) + f_R_patch;
        end
        
        function fill_interior(obj, f)
        % fills interior of function
            obj.f_R(obj.interior_idxs) = f(obj.R_X(obj.interior_idxs), obj.R_Y(obj.interior_idxs));        
        end
        
        function [f_xy, in_range] = locally_compute(obj, x, y, M)
            
            if x > obj.x_end || x < obj.x_start || y > obj.y_end || y < obj.y_start
                f_xy = nan;
                in_range = false;
                warning('(x, y) not in range');
                return
            end
            
            in_range = true;
           
            % j to the immediate left of point
            x_j = floor((x-obj.x_start)/obj.h);
            y_j = floor((y-obj.y_start)/obj.h);
            
            half_M =  floor(M/2);
            if mod(M, 2) ~= 0
                interpol_x_j_mesh = transpose(x_j-half_M:x_j+half_M);
                interpol_y_j_mesh = transpose(y_j-half_M:y_j+half_M);
            else
                interpol_x_j_mesh = transpose(x_j-half_M+1:x_j+half_M);
                interpol_y_j_mesh = transpose(y_j-half_M+1:y_j+half_M);
            end

            interpol_x_j_mesh = shift_idx_mesh(interpol_x_j_mesh, 0, obj.n_x-1);
            interpol_y_j_mesh = shift_idx_mesh(interpol_y_j_mesh, 0, obj.n_y-1);
            
            interpol_x_mesh = obj.h*interpol_x_j_mesh + obj.x_start;
            interpol_y_mesh = obj.h*interpol_y_j_mesh + obj.y_start;
            
            % first 1D interpolation
            interpol_x_exact = zeros(M, 1);
            for horz_idx = 1:M
                interpol_val = obj.f_R(interpol_y_j_mesh(horz_idx)+1, interpol_x_j_mesh+1)';
                interpol_x_exact(horz_idx) = barylag([interpol_x_mesh, interpol_val], x);
            end
             % second 1D interpolation
            f_xy = barylag([interpol_y_mesh, interpol_x_exact], y);
        end
        
        function compute_fc_coeffs(obj)
            obj.fc_coeffs = fftshift(fft2(obj.f_R) / numel(obj.f_R));
        end
        
        function [R_X_err, R_Y_err, f_interpolation, interior_idx] = ifft_interpolation(obj, rho_err)
            h_err = obj.h/rho_err;
            n_x_err = round((obj.x_end+obj.h-h_err-obj.x_start)/h_err)+1;
            n_y_err = round((obj.y_end+obj.h-h_err-obj.y_start)/h_err)+1;
            
            x_err_mesh = transpose(round((linspace(obj.x_start, obj.x_end+obj.h-h_err, n_x_err)-obj.x_start)/h_err)*h_err+obj.x_start);
            y_err_mesh = transpose(round((linspace(obj.y_start, obj.y_end+obj.h-h_err, n_y_err)-obj.y_start)/h_err)*h_err+obj.y_start);
            
            [R_X_err, R_Y_err] = meshgrid(x_err_mesh, y_err_mesh);

            n_x_diff = n_x_err - obj.n_x;
            n_y_diff = n_y_err - obj.n_y;
            
            padded_fc_coeffs = [zeros(ceil(n_y_diff/2), n_x_err); zeros(size(obj.fc_coeffs, 1), ceil(n_x_diff/2)) obj.fc_coeffs zeros(size(obj.fc_coeffs, 1), floor(n_x_diff/2)); zeros(floor(n_y_diff/2), n_x_err)];
            f_interpolation = real((n_x_err*n_y_err)*ifft2(ifftshift(padded_fc_coeffs)));
            
            idxs = reshape(1:numel(R_X_err), size(R_X_err));
            interior_idx = idxs(inpolygon_mesh(R_X_err, R_Y_err, obj.boundary_X, obj.boundary_Y));
        end
        
        function [grad_X, grad_Y, f_hat_x, f_hat_y] = grad(obj, f_mesh)
            if mod(obj.n_x, 2) == 0
                kx = (2*pi./(obj.x_end-obj.x_start+obj.h)).*(-obj.n_x/2:(obj. n_x/2-1));
            else
                kx =  (2*pi/(obj.x_end-obj.x_start+obj.h)).*((-(obj.n_x-1)/2):((obj.n_x-1)/2));
            end
            
            if mod(obj.n_y, 2) == 0
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*(-obj.n_y/2:(obj. n_y/2-1));
            else
                
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*((-(obj.n_y-1)/2):((obj.n_y-1)/2));
            end
            
            % Create 2D mesh for wavenumbers
            [KX, KY] = meshgrid(kx, ky);
            
            fft_coeffs = fftshift(fft2(f_mesh) / numel(obj.f_R));

            % Compute spectral derivatives
            f_hat_x = 1i * KX .* fft_coeffs;  % Fourier transform of ∂f/∂x
            f_hat_y = 1i * KY .* fft_coeffs;  % Fourier transform of ∂f/∂y

            % Transform back to real space
            grad_X = numel(obj.f_R)*real(ifft2(ifftshift(f_hat_x)));
            grad_Y = numel(obj.f_R)*real(ifft2(ifftshift(f_hat_y)));
        end
        
       function [f_div] = div(obj, f_1_mesh, f_2_mesh)
            if mod(obj.n_x, 2) == 0
                kx = (2*pi./(obj.x_end-obj.x_start+obj.h)).*(-obj.n_x/2:(obj. n_x/2-1));
            else
                kx =  (2*pi/(obj.x_end-obj.x_start+obj.h)).*((-(obj.n_x-1)/2):((obj.n_x-1)/2));
            end
            
            if mod(obj.n_y, 2) == 0
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*(-obj.n_y/2:(obj. n_y/2-1));
            else
                
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*((-(obj.n_y-1)/2):((obj.n_y-1)/2));
            end
            
            % Create 2D mesh for wavenumbers
            [KX, KY] = meshgrid(kx, ky);
            
            fft_coeffs_1 = fftshift(fft2(f_1_mesh) / numel(obj.f_R));
            fft_coeffs_2 = fftshift(fft2(f_2_mesh) / numel(obj.f_R));

            % Compute spectral derivatives
            f_hat = 1i* KX .* fft_coeffs_1 + 1i * KY .* fft_coeffs_2;  % Fourier transform of ∂f/∂x

            % Transform back to real space
            f_div = numel(obj.f_R)*real(ifft2(ifftshift(f_hat)));
       end
       
        
      function [f_lap, f_hat] = lap(obj, f_mesh)
            if mod(obj.n_x, 2) == 0
                kx = (2*pi./(obj.x_end-obj.x_start+obj.h)).*(-obj.n_x/2:(obj. n_x/2-1));
            else
                kx =  (2*pi/(obj.x_end-obj.x_start+obj.h)).*((-(obj.n_x-1)/2):((obj.n_x-1)/2));
            end
            
            if mod(obj.n_y, 2) == 0
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*(-obj.n_y/2:(obj. n_y/2-1));
            else
                
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*((-(obj.n_y-1)/2):((obj.n_y-1)/2));
            end
            
            % Create 2D mesh for wavenumbers
            [KX, KY] = meshgrid(kx, ky);
            
            fft_coeffs = fftshift(fft2(f_mesh) / numel(obj.f_R));

            % Compute spectral derivatives
            f_hat = -1* (KX.^2 + KY.^2).* fft_coeffs;

            % Transform back to real space
            f_lap = numel(obj.f_R)*real(ifft2(ifftshift(f_hat)));
      end
       function [f_lap, f_hat] = lap_coeff(obj, fft_coeffs)
            if mod(obj.n_x, 2) == 0
                kx = (2*pi./(obj.x_end-obj.x_start+obj.h)).*(-obj.n_x/2:(obj. n_x/2-1));
            else
                kx =  (2*pi/(obj.x_end-obj.x_start+obj.h)).*((-(obj.n_x-1)/2):((obj.n_x-1)/2));
            end
            
            if mod(obj.n_y, 2) == 0
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*(-obj.n_y/2:(obj. n_y/2-1));
            else
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*((-(obj.n_y-1)/2):((obj.n_y-1)/2));
            end
            
            % Create 2D mesh for wavenumbers
            [KX, KY] = meshgrid(kx, ky);
            
            % Compute spectral derivatives
            f_hat = -1* (KX.^2 + KY.^2).* fft_coeffs;

            % Transform back to real space
            f_lap = numel(obj.f_R)*real(ifft2(ifftshift(f_hat)));
       end
       
       function [f_inv_lap] = inv_lap(obj)
           if mod(obj.n_x, 2) == 0
                kx = (2*pi./(obj.x_end-obj.x_start+obj.h)).*(-obj.n_x/2:(obj. n_x/2-1));
            else
                kx =  (2*pi/(obj.x_end-obj.x_start+obj.h)).*((-(obj.n_x-1)/2):((obj.n_x-1)/2));
            end
            
            if mod(obj.n_y, 2) == 0
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*(-obj.n_y/2:(obj. n_y/2-1));
            else
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*((-(obj.n_y-1)/2):((obj.n_y-1)/2));
            end
            
            % Create 2D mesh for wavenumbers
            [KX, KY] = meshgrid(kx, ky);
            fc_coeffs_modified = ifftshift(obj.fc_coeffs); fc_coeff_0_0 = fc_coeffs_modified(1, 1); 
            fc_coeffs_modified = -fc_coeffs_modified./(ifftshift(KX.^2) + ifftshift(KY.^2));
            fc_coeffs_modified(1, 1) = 0; 
            f_inv_lap = numel(obj.f_R)*ifft2(fc_coeffs_modified) + fc_coeff_0_0*(obj.R_X.^2+obj.R_Y.^2)/4;
       end
       
       function filter = fft_filter(obj)
           if mod(obj.n_x, 2) == 0
                kx = (2*pi./(obj.x_end-obj.x_start+obj.h)).*(-obj.n_x/2:(obj. n_x/2-1));
            else
                kx =  (2*pi/(obj.x_end-obj.x_start+obj.h)).*((-(obj.n_x-1)/2):((obj.n_x-1)/2));
            end
            
            if mod(obj.n_y, 2) == 0
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*(-obj.n_y/2:(obj. n_y/2-1));
            else
                
                ky = (2*pi/(obj.y_end-obj.y_start+obj.h)).*((-(obj.n_y-1)/2):((obj.n_y-1)/2));
            end
            
            % Create 2D mesh for wavenumbers
            [KX, KY] = meshgrid(kx, ky);
            alpha = -log(1e-16);
            p = 18;
            filter = exp(-alpha*(2*KX/obj.n_x).^(2*p)).*exp(-alpha*(2*KY/obj.n_y).^(2*p));            
       end
    end
end

