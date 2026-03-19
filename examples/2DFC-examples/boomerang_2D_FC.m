clc; clear; close all;

f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));

d = 7;
C = 27;
n_r = 6;

M = d+3;
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

alpha = 3/2;
beta = tan(alpha.*pi./2);

l_1 = @(theta) -2/3*sin(theta*3*pi);
l_2 = @(theta) beta*sin(theta*2*pi);
l_1_prime = @(theta) -2*pi*cos(theta*3*pi);
l_2_prime = @(theta) 2*pi*beta*cos(theta*2*pi);
l_1_dprime = @(theta) 6*pi^2*sin(theta*3*pi);
l_2_dprime = @(theta) -4*pi^2*beta*sin(theta*2*pi);

h = 0.01;

curve_seq = Curve_seq_obj();

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, 0, 0.1, 0.1, 0.67, 0.67, h);

[R, ~, FC_patches, ~] = FC2D(f, h, curve_seq, 1e-13, 1e-13, d, C, n_r, A, Q, C, A, Q, M)
