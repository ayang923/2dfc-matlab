clc; clear; close all;

f = @(x, y) (3*pi)^2*sin(3*pi*x-1).*sin(3*pi*y-1);
u_boundary = @(x, y) -1/2.*sin(3*pi*x-1).*sin(3*pi*y-1);

alph = 1/2;
bet = tan(alph*pi/2);

h = 0.001;
curve_seq = Curve_seq_obj();

l_1 = @(theta) 2*sin(0.25*theta*pi);
l_2 = @(theta) -bet*sin(0.25*theta*2*pi);
l_1_prime = @(theta) 2*0.25*pi*cos(0.25*theta*pi);
l_2_prime = @(theta) -2*0.25*bet*pi*cos(0.25*theta*2*pi);
l_1_dprime = @(theta) -2*0.25^2*pi^2*sin(0.25*theta*pi);
l_2_dprime = @(theta) 4*bet*0.25^2*pi^2*sin(0.25*theta*2*pi);

h_tan = 2*h;
h_norm = h_tan;
n_curve = 0;

n_frac_C_0 = 0.1;
n_frac_C_1 = 0.2;
n_frac_S_0 = 0.6;
n_frac_S_1 = 0.7;

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n_curve, n_frac_C_0, n_frac_C_1, n_frac_S_0, n_frac_S_1, h_norm);

x_1 = [l_1(1); l_2(1)];
x_2 = [1; 0];
n = [x_2(2) - x_1(2); x_1(1) - x_2(1)]; n = n./norm(n);

c = 1/2*(x_1+x_2) + 0.1*n;

l_1 = @(theta) (1-theta).^2.*x_1(1)+2*(1-theta).*theta.*c(1)+theta.^2.*x_2(1);
l_2 = @(theta) (1-theta).^2.*x_1(2)+2*(1-theta).*theta.*c(2)+theta.^2.*x_2(2);
l_1_prime = @(theta) 2*((theta-1)*x_1(1)+(1-2*theta)*c(1)+theta*x_2(1));
l_2_prime = @(theta) 2*((theta-1)*x_1(2)+(1-2*theta)*c(2)+theta*x_2(2));
l_1_dprime = @(theta) 2*x_1(1)-4*c(1)+2*x_2(1);
l_2_dprime = @(theta) 2*x_1(2)-4*c(2)+2*x_2(2);

h_tan = 2*h;
h_norm = h_tan;
n_curve = 0;

n_frac_C_0 = 0.3;
n_frac_C_1 = 0.3;
n_frac_S_0 = 0.7;
n_frac_S_1 = 0.7;

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n_curve, n_frac_C_0, n_frac_C_1, n_frac_S_0, n_frac_S_1, h_norm);

l_1 = @(theta) 2*sin(theta*pi);
l_2 = @(theta) -bet*sin(theta*2*pi);

x_1 = [1; 0];
x_2 = [l_1(0.75); l_2(0.75)];
n = [x_2(2) - x_1(2); x_1(1) - x_2(1)]; n = n./norm(n);

c = 1/2*(x_1+x_2) - 0.1*n;

l_1 = @(theta) (1-theta).^2.*x_1(1)+2*(1-theta).*theta.*c(1)+theta.^2.*x_2(1);
l_2 = @(theta) (1-theta).^2.*x_1(2)+2*(1-theta).*theta.*c(2)+theta.^2.*x_2(2);
l_1_prime = @(theta) 2*((theta-1)*x_1(1)+(1-2*theta)*c(1)+theta*x_2(1));
l_2_prime = @(theta) 2*((theta-1)*x_1(2)+(1-2*theta)*c(2)+theta*x_2(2));
l_1_dprime = @(theta) 2*x_1(1)-4*c(1)+2*x_2(1);
l_2_dprime = @(theta) 2*x_1(2)-4*c(2)+2*x_2(2);

h_tan = 2*h;
h_norm = h_tan;
n_curve = 0;

n_frac_C_0 = 0.3;
n_frac_C_1 = 0.3;
n_frac_S_0 = 0.7;
n_frac_S_1 = 0.7;

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n_curve, n_frac_C_0, n_frac_C_1, n_frac_S_0, n_frac_S_1, h_norm);

l_1 = @(theta) 2*sin((0.25*theta+0.75)*pi);
l_2 = @(theta) -bet*sin((0.25*theta+0.75)*2*pi);
l_1_prime = @(theta) 2*0.25*pi*cos((0.25*theta+0.75)*pi);
l_2_prime = @(theta) -2*0.25*bet*pi*cos((0.25*theta+0.75)*2*pi);
l_1_dprime = @(theta) -2*0.25^2*pi^2*sin((0.25*theta+0.75)*pi);
l_2_dprime = @(theta) 4*bet*0.25^2*pi^2*sin((0.25*theta+0.75)*2*pi);

h_tan = 2*h;
h_norm = h_tan;
n_curve = 0;

n_frac_C_0 = 0.3;
n_frac_C_1 = 0.1;
n_frac_S_0 = 0.7;
n_frac_S_1 = 0.7;

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n_curve, n_frac_C_0, n_frac_C_1, n_frac_S_0, n_frac_S_1, h_norm);

d = 7;
C = 27;
n_r = 6;
M = d+3;
p = M;


if(exist(['FC_data/A_d',num2str(d),'_C', num2str(C), '_r', num2str(n_r), '.mat']) == 0 | ...
   exist(['FC_data/Q_d',num2str(d),'_C', num2str(C),  '_r', num2str(n_r), '.mat']) == 0)
    disp('FC data not found. Generating FC operators... \n');
    generate_bdry_continuations(d, C, C, 12, 20, 4, ...
        256, n_r);
end

load(['FC_data/A_d',num2str(d),'_C', num2str(C), '_r', num2str(n_r), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C', num2str(C),  '_r', num2str(n_r), '.mat']);

A = double(A);
Q = double(Q);

h_eval = 0.02;
[u_num_mat, R, R_eval] = poisson_solver_coarse(curve_seq, f, u_boundary, h, 1, p, 1e-14, 1e-13, 1e-13, d, C, n_r, A, Q, M, h_eval);

% Plotting using triangulation to reduce boundary stepping
% Interior vertices (from your eval grid)
Xint = R_eval.R_X(:);
Yint = R_eval.R_Y(:);
Zint = u_num_mat(:);

% Boundary vertices
[Xb, Yb] = curve_seq.construct_boundary_mesh(1);
Zb = u_boundary(Xb, Yb);

% Make sure boundary is a closed loop (append first point if needed)
if ~(Xb(1)==Xb(end) && Yb(1)==Yb(end))
    Xb = [Xb; Xb(1)];
    Yb = [Yb; Yb(1)];
    Zb = [Zb; Zb(1)];
end

% Build polygon for inside/outside test (nonconvex OK)
P = polyshape(Xb, Yb);

% Combine points
X = [Xint; Xb];
Y = [Yint; Yb];
Z = [Zint; Zb];

% Constrained edges for boundary (connect consecutive boundary points)
nint = numel(Xint);
nb   = numel(Xb);
bidx = (nint+1):(nint+nb);
E    = [bidx(1:end-1)' bidx(2:end)'];   % boundary edges (already closed)

% Constrained Delaunay triangulation in (x,y)
DT = delaunayTriangulation(X, Y, E);
T  = DT.ConnectivityList;

% Keep only triangles whose centroid lies inside the polygon
xc = mean(X(T), 2);
yc = mean(Y(T), 2);
inside = isinterior(P, xc, yc);
T_in = T(inside, :);

% Plot triangulated surface embedded in 3D (boundary is part of the mesh)
figure;
trisurf(T_in, X, Y, Z, 'EdgeColor','k');  % change 'k' -> 'none' for no mesh lines
shading interp;

hold on;

plot(R_eval.boundary_X, R_eval.boundary_Y, 'LineWidth', 2)