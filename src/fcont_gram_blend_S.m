% fcont_gram_blend_S : This routine computes the Fourier continuation in
% parameter space for a S type patch, for a function over the unit square
% domain. The function blends towards 0 upwards because the top row usually
% corresponds with the x axis.
%
% Inputs: 
%   fx : (real) values of the functions f over the unit square domain given
%   a 2D matrix
%   d : number of Gram interpolation points 
%   A : matrix containing the values of the extended Gram polynomials at 
%   the continuation points
%   Q : matrix containing the values of the Gram polynomials at the first 
%        d points of the interval. 

%
% Outputs:
%   fcont : (real) vector of size lenfgth(fx) + C containing the values of
%   the fx and its continuation on the extended grid
%   fc: real vector of size C containing only the values of the continuation

function [fcont] = fcont_gram_blend_S(fx, d, A, Q)
fl = fx(1:d, :);

fc = flipud(A * (Q .' * flipud(fl)));
fc = double(fc);

fcont = [fc; fx(1, :)];
end

