function [xn, converged] = newton_solve(f,J,x0,tau,Nmax)
% Solves the nonlinear equation given by f(x) = 0 via Newton iteration
% Inputs: 
%    f         function to find root of
%    J         Jacobian of f
%    x0        initial guess
%    tau       convergence tolerance
%    Nmax      maximum number of iterations
% Outputs:
%    xn        matrix of newton iterates (xn(:,end) is the final answer)

xn = x0;

for i = 2:Nmax
    xprev = xn;
    xn = xn - J(xn)\f(xn);

    if max(abs(f(xn))) < tau && max(abs(xprev-xn)) < tau
        break
    end
end

if i == Nmax
    converged = false;
else
    converged =  true;
end
end
