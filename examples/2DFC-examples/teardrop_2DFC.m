clc; clear; close all;

f = @(x, y) -exp(0.5*(x.^2+y.^2)).*(sin(3*pi*x)+cos(2.5*pi*y));
% f = @(x, y) 100*sin(10*(x-1)).*sin(10*(y-1));

alph = 0.5;
bet = tan(alph*pi/2);

l_1 = @(theta) 2*sin(theta*pi);
l_2 = @(theta) -bet*sin(theta*2*pi);
l_1_prime = @(theta) 2*pi*cos(theta*pi);
l_2_prime = @(theta) -2*bet*pi*cos(theta*2*pi);
l_1_dprime = @(theta) -2*pi^2*sin(theta*pi);
l_2_dprime = @(theta) 4*bet*pi^2*sin(theta*2*pi);

d = 8;
C = 27;
n_r = 6;
M = d+3;

h = 0.01;

h_norm = 1.5*h;
n_frac_C = 0.1;
n_frac_S = 0.6;
n_curve = 0; % default

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

curve_seq = Curve_seq_obj();
curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n_curve, n_frac_C, n_frac_C, n_frac_S, n_frac_S, h_norm);

curve_seq.plot_geometry(d)
R = FC2D(f, h, curve_seq, 5e-15, 5e-15, d, C, n_r, A, Q, C, A, Q, M);

