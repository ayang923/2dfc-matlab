
% A routine used to build save the matrices necessary for the
% generation of the Fourier Continuation.
%
%
% Author:
%
%   Daniel Leibovici
%   email: dleibovi@caltech.edu
%

function generate_bdry_continuations(d, C, E, Z, n_ovr, modes_to_reduce, ...
    num_digits, n_r)

tic;
fprintf('Performing precomputations...\n');
[Q, A] =  ... 
    precomp_fc_data(d, C, Z, E, n_ovr, modes_to_reduce, num_digits, n_r);

this_file = mfilename('fullpath');
this_folder = fileparts(this_file);

save(fullfile(this_folder, ['FC_data/A_d',num2str(d),'_C', num2str(C), '_r', num2str(n_r), '.mat']), 'A');
save(fullfile(this_folder, ['FC_data/Q_d',num2str(d),'_C', num2str(C),  '_r', num2str(n_r), '.mat']), 'Q');

% Saving to Ascii
A = double(A);
Q = double(Q);
dlmwrite(fullfile(this_folder, ['FC_data/A_d',num2str(d),'_C', num2str(C), '_r', num2str(n_r),'.txt']), A, 'precision', 20);
dlmwrite(fullfile(this_folder, ['FC_data/Q_d',num2str(d),'_C', num2str(C), '_r', num2str(n_r), '.txt']), Q, 'precision', 20)

% save(['FC_data/ArQr_d',num2str(d),'_C', num2str(C), '.mat'], 'ArQr');
% save(['FC_data/AlQl_d',num2str(d),'_C', num2str(C), '.mat'], 'AlQl');
% save(['FC_data/ArQ_tilder_d',num2str(d),'_C', num2str(C), '.mat'], 'ArQ_tilder');
% save(['FC_data/AlQ_tildel_d',num2str(d),'_C', num2str(C), '.mat'], 'AlQ_tildel');
% 
% 
% dlmwrite(['FC_data/A', num2str(d), 'C', num2str(C), '.dat'], A, ...
%          'delimiter', ',', 'precision', 30);
% dlmwrite(['FC_data/Q', num2str(d), 'C', num2str(C), '.dat'], Q, ...
%          'delimiter', ',', 'precision', 30);
% dlmwrite(['FC_data/Q_tilde', num2str(d), 'C', num2str(C), '.dat'], Q_tilde, ...
%  'delimiter', ',', 'precision', 30);
% dlmwrite(['FC_data/ArQr', num2str(d), 'C', num2str(C), '.dat'], ArQr, ...
%      'delimiter', ',', 'precision', 30);
% dlmwrite(['FC_data/AlQl', num2str(d), 'C', num2str(C), '.dat'], AlQl, ...
%      'delimiter', ',', 'precision', 30);   
% dlmwrite(['FC_data/ArQ_tilder', num2str(d), 'C', num2str(C), '.dat'], ArQ_tilder, ...
%      'delimiter', ',', 'precision', 30); 
% dlmwrite(['FC_data/AlQ_tildel', num2str(d), 'C', num2str(C), '.dat'], AlQ_tildel, ...
%      'delimiter', ',', 'precision', 30);       
    
toc;

end
