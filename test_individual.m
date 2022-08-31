%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : test_individual.m                                             %
%                                                                         %
% Author  : Tobias Holicki                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is for the purpose of testing some LMI based iteration for 
% static output-feedback Hinfty design.
%  ! Note that the default solvers is SDPT3 called from Yalmip, but     !
%  ! using Mosek leads to much better results                           !


% Clean up.
clc
clear


%% System data
% Load system data
% Performance for fo_iteration is poor, e.g., for AC3, AC4
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('DIS1');

% System
sys = ss(A, [B1, B], [C1; C], [D11, D12; D21, zeros(ny, nu)]);
%sys = c2d(sys, 0.02); % For discrete-time tests




%% Construct some static controllers

% *hinfsyn* 
% The produced upper bound by this function is a lower bound on the
% closed-loop energy gain that can be achieved by the static designs.
[~, ~, ga] = hinfsyn(sys, ny, nu);
disp(['Lower bound by dynamic output-feedback: ', num2str(ga)])


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


% *primal dual iteration*
disp('==================== primal dual iteration ======================')
%primal_dual_iteration(sys, ny, nu, disp=true, con=con);
[~, ga] = primal_dual_iteration(sys, ny, nu, con=con, max_ite_ph2=5);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 5);
[~, ga] = primal_dual_iteration(sys, ny, nu, con=con, max_ite_ph2=10);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 10);


% *primal iteration*
disp('====================== primal iteration ======================')
[~, ga] = primal_iteration(sys, ny, nu, con=con, max_ite_ph2=50);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 50);
[~, ga] = primal_iteration(sys, ny, nu, con=con, max_ite_ph2=100);
fprintf('Final: Peak gain = %0.3f, Iterations = %i\n', ga, 100);


