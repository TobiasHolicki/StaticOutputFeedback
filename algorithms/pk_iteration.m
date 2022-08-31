%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : pk_iteration.m                                                %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 02                                                            %
% Date    : 02.08.2022                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a given continuous- or discrete-time LTI system
% dot/sigma x = A x +   Bd +  B2u
%           e = C x +   Dd + D12u
%           y = C2x + D21d  
% this function tries to computes a static output-feedback controller
%     u = Ky
% such that the energy gain of the closed-loop system of the channel 
% d -> e is as small as possible. 
% 
% Instead of searching for static controllers one can also search other 
% fixed order controllers whose describing matrices K := [Ak, Bk; Ck, Dk]
% satisfy arbitrary affine inequality or equality constraints:
%            A * K(:) <= b     and/or     Aeq * K(:) = beq.
%
% The design procedure is essentially the classical P-K (or D-K) iteration
% which is based on iteratively solving the closed-loop LMIs corresponding
% to the bounded real lemma. The idea is to alternatingly fix the involved 
% Lyapunov matrix and the describing matrices of the controller. 
%
% See also svar_iteration, fo_iteration, primal_dual_iteration, 
% primal_iteration
%
% ----- Input ---------------------------------------------------------- 
%   
% sys           - The given LTI system
% mea           - Number of measurements, i.e., number of rows in C2.
% act           - Number of actuators, i.e., number of columns in B2.
% lxc           - Desired state dimension of the controller.
% op            - Struct
%   max_ite_ph1 - Maximum number of iterations during the feasibility 
%                 phase of the algorithm.
%   max_ite_ph2 - Maximum number of iterations during the second phase
%                 of the algorithm aiming at minimizing the energy gain.
%   stp_slw_prg - The second phase terminates if the computed upper bounds
%                 do not decrease by more than 1% in the last stp_slw_prg
%                 iterations.
%   opt         - Yalmip solver options.
%   disp        - Displays progress and some additional information.
%   eps         - All inequalities are rendered strict with this number.
%   Aep, beq    - These specify the equality constraints on the controller.
%   A, b        - These specify the inequality constraints on the
%                 controller matrices.
%   con         - Some a priori known stabilizing controller.
% ----- Output ---------------------------------------------------------
%   con         - Computed controller
%   ga          - Computed upper bound on the closed-loop energy gain.
%  
function [con, ga] = pk_iteration(sys, mea, act, lxc, op)
    % Some sanity checks
    arguments
        sys {mustBeA(sys, "ss")}
        mea (1, 1) {mustBeInteger, mustBeNonnegative}
        act (1, 1) {mustBeInteger, mustBeNonnegative}
        lxc (1, 1) {mustBeInteger, mustBeNonnegative} = 0
        op.max_ite_ph1 (1, 1) {mustBeInteger, mustBeNonnegative} = 200
        op.max_ite_ph2 (1, 1) {mustBeInteger, mustBeNonnegative} = 300
        op.stp_slw_prg (1, 1) {mustBeInteger, mustBeNonnegative} = 100
        op.opt = sdpsettings('solver', 'sdpt3', 'verbose', 0);
        op.disp {mustBeA(op.disp, "logical")} = false
        op.eps (1, 1) {mustBeNumeric} = 1e-6
        op.Aeq = []
        op.beq = []
        op.A   = []
        op.b   = []
        op.con {mustBeA(op.con, "ss")} = ss([]);
    end
    
    % Override display function of no display is desired
    if op.disp
        displ = @(x) disp(x);
    else
        displ = @(x) 0; 
    end
    
    % Skip first phase if a stabilizing controller is provided
    if ~isempty(op.con)
        op.max_ite_ph1 = 0; 
        Kv = [op.con.a, op.con.b; op.con.c, op.con.d];
    end
    
    %% Abbreviations

    % *Some lengths*
    lx  = size(sys.a, 1);     % Number of system states
    lcl = lx + lxc;           % Number of closed-loop system states
    err = size(sys, 1) - mea; % Size of error signal
    dis = size(sys, 2) - act; % Size of generalized disturbance

    % *Closed-loop matrices*
    [A, B, C, D] = sssdata(sys, [dis, act], [err, mea]);
    Ac = @(K) blkdiag(A, zeros(lxc)) + blkodiag(B{2}, eye(lxc)) * K * ...
                                                  blkodiag(eye(lxc), C{2});
    Bc = @(K) [B{1}; zeros(lxc, dis)] + blkodiag(B{2}, eye(lxc)) * K * ...
                                                [zeros(lxc, dis); D{2, 1}];
    Cc = @(K) [C{1}, zeros(err, lxc)] + [zeros(err, lxc), D{1, 2}] * K *...
                                                  blkodiag(eye(lxc), C{2});
    Dc = @(K) D{1, 1} + [zeros(err, lxc), D{1, 2}] * K * ...
                                                [zeros(lxc, dis); D{2, 1}];
    if norm(D{2,2}) > 1e-10
       error(['The closed-loop must depend on the controller matrices', ...
              ' in an affine fashion, i.e., D22 = 0']); 
    end
    
    %% Preparations for both phases

    % *Define variables*
    r  = sdpvar(1, 1); % Bound to feasibility
    ga = sdpvar(1, 1); % Upper bound on the energy gain
    X  = sdpvar(lcl);  % Lyapunov certificate
    K  = sdpvar(lxc+act, lxc+mea, 'full'); % Controller parameters

    % *Constraints for the problems with fixed controller*
    if sys.Ts == 0 % Continuous time
        Te = [Ac(K)' * X + X * Ac(K) - r * eye(lcl), X * Bc(K), Cc(K)'; ...
              (X * Bc(K))', -ga * eye(dis), Dc(K)'; ...
              Cc(K), Dc(K), -ga * eye(err)];
        Con2 = [Te <= -op.eps * eye(lcl+err+dis)];
    else
        Te = [-X, X * Ac(K), X * Bc(K), zeros(lcl, err); ...
              (X * Ac(K))', -X - r * eye(lcl), zeros(lcl, dis), Cc(K)'; ...
              (X * Bc(K))', zeros(dis, lcl), -ga * eye(dis), Dc(K)'; ...
              zeros(err, lcl), Cc(K), Dc(K), -ga * eye(err)];
        Con2 = [Te <= -op.eps * eye(2*lcl+err+dis)];
    end
    Con2 = [Con2; X  >=  op.eps * eye(lcl)];

    % *Constraints for the problems with fixed Lyapunov certificate*
    Con1 = Con2;
    % General affine constraints on the controller matrices
    if ~isempty(op.A) && ~isempty(op.b)
        Con1 = [Con1; op.A * K(:) <= op.b];
    end
    if ~isempty(op.Aeq) && ~isempty(op.beq)
        Con1 = [Con1; op.Aeq * K(:) == op.beq];
    end

    % *Optimization problems*
    % For the first phase we use a very small penalty on the upper 
    % bound ga in order to avoid starting the second phase with a
    % gigantic energy gain.
    YOP1 = optimizer(Con1, r+0.001*ga, op.opt, X, {r, ga, K});
    YOP2 = optimizer(Con2, r+0.001*ga, op.opt, K, {r, ga, X});
    YOP3 = optimizer(Con1, ga, op.opt, {r, X}, {ga, K});
    YOP4 = optimizer(Con2, ga, op.opt, {r, K}, {ga, X});
    
    %% Phase 1 (Trying to find a stabilizing controller)
    
    displ('**Starting phase 1**');
    
    % *Iteratively solve LMIs*
    Xv  = eye(lcl) + 0.001*ones(lcl); % Initial certificate
    ind = 1;
    while ind <= op.max_ite_ph1
        % Solve an LMI problem with fixed controller
        sol = YOP1(Xv);
        rv  = sol{1};
        Kv  = sol{3};

        % Display progress
        displ(sprintf('Step: %i, Bound to feasibility: %0.5f', ind, rv));
        ind = ind + 1;

        % Stop if we have found a stabilizing controller
        if rv < 0; break; end

        % Solve an LMI problem with fixed controller matrices
        sol = YOP2(Kv);
        rv  = sol{1};
        Xv  = sol{3};

        % Display progress
        displ(sprintf('Step: %i, Bound to feasibility: %0.5f', ind, rv));
        ind = ind + 1;

        % Stop if we have found a stabilizing controller
        if rv < 0; break; end

        if ind >= op.max_ite_ph1 % Then we were unsuccessfull
            error(['No stabilizing controller satisfying the ', ...
                   'given constraints was found. Try to weaken ', ...
                   'those or to increase max_ite_ph1']);
        end
    end
    % *Compute energy gain and get a reasonable Lyapunov certificate*    
    sol = YOP4({0, Kv});
    gav = sol{1};
    Xv  = sol{2}; 
    
    displ(sprintf(['Found a stabilizing controller satisfying the ', ...
           'constraints \nwith the closed-loop energy-gain: ', ...
           num2str(gav)]));
       
    %% Phase 2 (Minimize the closed-loop energy gain)
    
    displ('**Starting phase 2**');
  
    % *Iteratively solve LMIs*
    ind = 1;
    gav = NaN * ones(op.max_ite_ph2, 1);
    while ind <= op.max_ite_ph2
        % Solve an LMI problem with fixed slack variable
        sol      = YOP3({0, Xv});
        gav(ind) = sol{1};
        Kv       = sol{2};

        % Display progress
        displ(sprintf('gamma(%i) = %0.5f', ind, gav(ind)));
        ind = ind + 1;

        % Solve an LMI problem with fixed controller matrices       
        sol      = YOP4({0, Kv});
        gav(ind) = sol{1};
        Xv       = sol{2};

        displ(sprintf('gamma(%i) = %0.5f', ind, gav(ind)));
        ind = ind + 1;

        % Stop earlier when the progress is too slow
        if ind > op.stp_slw_prg && ...
              abs(gav(ind-op.stp_slw_prg) - gav(ind-1)) / ...
                 gav(ind-op.stp_slw_prg) < 0.01             
            displ(['Stopping since the progress is below 1% in ', ...
                   sprintf('the last %i iterations', op.stp_slw_prg)]);
            break;
        end
    end
 
    % Final upper bound
    ga = gav(ind-1);
    
    % Build final controller
    con = mat2ss(Kv, lxc, sys.Ts);
    
    % Perform a posteriori analysis
    displ(['The closed-loop energy gain is ', ...
           num2str(hinfnorm(lft(sys, con)))]);   
end
