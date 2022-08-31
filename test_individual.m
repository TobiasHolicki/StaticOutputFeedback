%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : test_individual.m                                             %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 05                                                            %
% Date    : 08.08.2022                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is for the purpose of testing some LMI based iteration for 
% static output-feedback Hinfty design.

% Clean up.
clc
clear


%% System data
% Load system data
% Performance for fo_iteration is poor, e.g., for AC3, AC4
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('AC6');

% System
sys = ss(A, [B1, B], [C1; C], [D11, D12; D21, zeros(ny, nu)]);
%sys = c2d(sys, 0.02); % For discrete-time tests

[~, ~, ga] = hinfsyn(sys, ny, nu);


[con, ga] = hinftest(sys, ny, nu);
ga
hinfnorm(lft(sys, con))

return

%% Construct some controllers


% *hinfstruct*
disp('====================== hinfstruct ======================')
% We search static controllers
Khis = tunableSS('con', 0, nu, ny, sys.Ts); 
opt  = hinfstructOptions;
opt.Randomstart = 0; % Turning the algorithm deterministic for comparison
hinfstruct(sys, Khis, opt);

% For testing other algorithms with the same initially stabilizing
% controller. (Adjust opt.MaxIter if hinfstruct does not manage to find 
% one) 
opt.MaxIter = 1;
con = hinfstruct(sys, Khis, opt);
con = ss(con);


% *svar iteration*
disp('====================== svar iteration ======================')
%svar_iteration(sys, ny, nu, disp=true, con=con);
[~, ga] = svar_iteration(sys, ny, nu, con=con, max_ite_ph2=50);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 50);
[~, ga] = svar_iteration(sys, ny, nu, con=con, max_ite_ph2=100);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 100);


% *pk iteration*
disp('====================== pk iteration ======================')
%pk_iteration(sys, ny, nu, disp=true, con=con);
[~, ga] = pk_iteration(sys, ny, nu, con=con, max_ite_ph2=50);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 50);
[~, ga] = pk_iteration(sys, ny, nu, con=con, max_ite_ph2=100);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 100);


% *fo iteration*
disp('====================== fo iteration ======================')
%fo_iteration(sys, ny, nu, disp=true, con=con);
[~, ga] = fo_iteration(sys, ny, nu, con=con, max_ite_ph2=50);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 50);
[~, ga] = fo_iteration(sys, ny, nu, con=con, max_ite_ph2=100);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 100);
[~, ga] = fo_iteration(sys, ny, nu, con=con, max_ite_ph2=200);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 200);


% *dual iteration*
disp('====================== dual iteration ======================')
%dual_iteration(sys, ny, nu, disp=true, con=con);
[~, ga] = dual_iteration(sys, ny, nu, con=con, max_ite_ph2=5);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 5);
[~, ga] = dual_iteration(sys, ny, nu, con=con, max_ite_ph2=10);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 10);


% *primal iteration*
disp('====================== primal iteration ======================')
%dual_iteration(sys, ny, nu, disp=true, con=con);
[~, ga] = primal_iteration(sys, ny, nu, con=con, max_ite_ph2=50);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 50);
[~, ga] = primal_iteration(sys, ny, nu, con=con, max_ite_ph2=100);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 100);


