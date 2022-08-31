%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : fo_iteration.m                                                %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 05                                                            %
% Date    : 03.08.2022                                                    %
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
% The design procedure is a variant of the one proposed in [1]. It relies
% on applying elimination twice in the direction which introduces
% additional (slack) variables [2]. 
%
% [1] Felipe, A. and Oliveira, R. C. L. F., An LMI-Based Algorithm to 
%     Compute Robust Stabilizing Feedback Gains Directly as Optimization 
%     Variables, 2021.
% [2] Ebihara, Y. et al., S-Variable Approach to LMI-Based Robust Control, 
%     2015.
%
% See also pk_iteration, svar_iteration, dual_iteration, primal_iteration
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
function [con, ga] = fo_iteration(sys, mea, act, lxc, op)
    % Some sanity checks
    arguments
        sys {mustBeA(sys, "ss")}
        mea (1, 1) {mustBeInteger, mustBeNonnegative}
        act (1, 1) {mustBeInteger, mustBeNonnegative}
        lxc (1, 1) {mustBeInteger, mustBeNonnegative} = 0
        op.max_ite_ph1 (1, 1) {mustBeInteger, mustBeNonnegative} = 200
        op.max_ite_ph2 (1, 1) {mustBeInteger, mustBeNonnegative} = 100
        op.stp_slw_prg (1, 1) {mustBeInteger, mustBeNonnegative} = 50
        op.opt = sdpsettings('solver', 'mosek', 'verbose', 0);
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
    r  = sdpvar(1, 1); % Bound to feasibility / relaxation parameter
    ga = sdpvar(1, 1); % Upper bound on the energy gain
    X  = sdpvar(lcl);  % Lyapunov certificate
    K  = sdpvar(lxc+act, lxc+mea, 'full');      % Controller parameters
    G  = sdpvar(3*lcl+2*err, lcl+err, 'full');  % Slack variable
    F  = sdpvar(3*lcl+2*err, lcl+err, 'full');  % Slack variable

    % *Constraints*
    % The first constraint results from applying elimination twice
    % to generate slack variables G (first appl.) and F (second appl.)
    Te = [-eye(lcl), Ac(K), zeros(lcl, err), Bc(K); ...
          zeros(err, lcl), Cc(K), -eye(err), Dc(K)];
    if sys.Ts == 0 % Continuous time
        Te = [blkdiag([zeros(lcl), X; X, -r*eye(lcl)], eye(err), ...
              -ga*eye(dis)), Te'; Te, zeros(lcl+err)];
    else % Discrete time
        Te = [blkdiag(X, -X - r*eye(lcl), eye(err), -ga*eye(dis)),...
              Te'; Te, zeros(lcl+err)];
    end
    L   = blkdiag(eye(2*lcl+err), [zeros(dis, lcl+err); eye(lcl+err)]);
    Te  = Te + (L * G) * (L * F)' + (L * F) * (L * G)'; 
    Con = [Te <= -op.eps * eye(3*lcl + 2*err + dis); ...
           X  >=  op.eps * eye(lcl)];
    % General affine constraints on the controller matrices
    if ~isempty(op.A) && ~isempty(op.b)
        Con = [Con; op.A * K(:) <= op.b];
    end
    if ~isempty(op.Aeq) && ~isempty(op.beq)
        Con = [Con; op.Aeq * K(:) == op.beq];
    end
    
    % *Optimization problems*
    % For the first phase we use a very small penalty on the upper 
    % bound ga in order to avoid starting the second phase with a
    % gigantic energy gain.
    YOP1 = optimizer(Con, r+0.001*ga, op.opt, G, {r, ga, X, K, F});
    YOP2 = optimizer(Con, ga, op.opt, {r, G}, {ga, X, K, F});
    
    %% Phase 1 (Trying to find a stabilizing controller)
    
    displ('**Starting phase 1**');
    
    % *Iteratively solve LMIs*
    Gv  = [eye(lcl), zeros(lcl, err); ...
           zeros(lcl, lcl+err);...
           zeros(err, lcl), eye(err); ...
           -eye(lcl+err)]; % Initial slack variable
    ind = 1;
    while ind <= op.max_ite_ph1
        % Solve an LMI problem with fixed slack variable G
        sol = YOP1(Gv);
        rv  = sol{1};
        Gv  = sol{5}; % = F, this is possible due to the Hermitian part 

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

    % *Compute energy gain and get a reasonable slack variable G for 
    %  continuing*    
    if isempty(op.con) % Obtained or provied controller
        Kv = sol{4};
    else
        Kv = [op.con.a, op.con.b; op.con.c, op.con.d];
    end
    
    % Part of the slack variable
    G  = sdpvar(2*lcl+err, lcl+err, 'full'); % Slack variable
    % Constraints (corresponds to applying elimination only once)
    L  = [G; zeros(dis, lcl+err)];
    R  = [-eye(lcl), Ac(Kv), zeros(lcl, err), Bc(Kv); ...
          zeros(err, lcl), Cc(Kv), -eye(err), Dc(Kv)];
    if sys.Ts == 0 % Continuous time
        Te = blkdiag([zeros(lcl), X; X, zeros(lcl)], eye(err), ...
                     -ga * eye(dis));
    else % Discrete time
        Te = blkdiag(X, -X, eye(err), -ga * eye(dis));
    end
    Con = [Te + (L * R) + (L * R)' <= -op.eps * eye(2*lcl+err+dis); ...
           X  >=  op.eps * eye(lcl)];
    % Optimize
    t = optimize(Con, ga, op.opt);

    % Determined solutions
    gav = sqrt(value(ga));
    Gv  = [value(G); -eye(lcl+err)]; % The new slack variable used as
                                     % initialization for phase 2 

    if t.problem == 0
        displ(sprintf(['Found a stabilizing controller satisfying the ',...
               'constraints \nwith closed-loop energy-gain: ', ...
               num2str(gav)]));
    else
        warning(['Numerical issues or infeasibility have occurred, ', ...
                 'I continue, but do not have expectations']);
    end
    %% Phase 2 (Minimize the closed-loop energy gain)
    
    displ('**Starting phase 2**');
  
    % *Iteratively solve LMIs*
    ind = 1;
    gav = NaN * ones(op.max_ite_ph2, 1);
    while ind <= op.max_ite_ph2
        % Solve an LMI problem with fixed slack variable G
        sol      = YOP2({0, Gv});
        gav(ind) = sqrt(sol{1});
        Gv       = sol{4}; % = F, possible due to the Hermitian part

        % Display progress
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
    con = mat2ss(sol{3}, lxc, sys.Ts);
    
    % Perform a posteriori analysis
    displ(['The closed-loop energy gain is ', ...
           num2str(hinfnorm(lft(sys, con)))]);   
end
