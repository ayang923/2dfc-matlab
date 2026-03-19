clc; clear; close all;

f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));
% f = @(x, y) 100*sin(10*(x-1)).*sin(10*(y-1));

alph = 1.99;
bet = tan(alph*pi/2);

l_1 = @(theta) bet*cos((1+alph)*pi*theta) - sin((1+alph)*pi*theta) - bet;
l_2 = @(theta) bet*sin((1+alph)*pi*theta) + cos((1+alph)*pi*theta) - cos(pi*theta);

l_1_prime  = @(theta) -bet*(1+alph)*pi*sin((1+alph)*pi*theta) ...
                      - (1+alph)*pi*cos((1+alph)*pi*theta);
l_2_prime  = @(theta)  bet*(1+alph)*pi*cos((1+alph)*pi*theta) ...
                      - (1+alph)*pi*sin((1+alph)*pi*theta) ...
                      + pi*sin(pi*theta);        

l_1_dprime = @(theta) -bet*(1+alph)^2*pi^2*cos((1+alph)*pi*theta) ...
                      + (1+alph)^2*pi^2*sin((1+alph)*pi*theta);
l_2_dprime = @(theta) -bet*(1+alph)^2*pi^2*sin((1+alph)*pi*theta) ...
                      - (1+alph)^2*pi^2*cos((1+alph)*pi*theta) ...
                      + pi^2*cos(pi*theta);        

d = 7;
C = 27;
n_r = 6;
M = d+3;


h = 0.0002;
h_tan =5.5*h;
h_norm = h_tan;

n_frac_C = 0.1;
n_frac_S = 0.6;

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
n_curve = ceil(3.5/h_tan);
curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n_curve, n_frac_C, n_frac_C, n_frac_S, n_frac_S, h_tan);

curve_seq.plot_geometry(d)
R = FC2D(f, h, curve_seq, 5e-15, 5e-15, d, C, n_r, A, Q, C, A, Q, M);

