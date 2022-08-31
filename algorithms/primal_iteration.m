%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : primal_iteration.m                                            %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 02                                                            %
% Date    : 08.08.2022                                                    %
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
% Instead of searching for standard static controllers one can also search 
% controllers satisfy arbitrary affine inequality or equality constraints:
%            A * K(:) <= b     and/or     Aeq * K(:) = beq.
%
% The design procedure is based on iteratively solving LMIs that involve
% a slack variable [1] and on alternatingly fixing this slack variable and
% the describing matrices of the controller. If the controller is fixed,
% then one essentialyl solves an analysis problem in which one could
% eliminate the slack variable. However, the reconstruction doesn't seem to
% produce nice slack variables for the subsequent LMI problems.
%
% [1] Ebihara, Y. et al., S-Variable Approach to LMI-Based Robust Control, 
%     2015.
%
% See also pk_iteration, svar_iteration, fo_iteration, dual_iteration
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
function [con, ga] = primal_iteration(sys, mea, act, op)
    % Some sanity checks
    arguments
        sys {mustBeA(sys, "ss")}
        mea (1, 1) {mustBeInteger, mustBeNonnegative}
        act (1, 1) {mustBeInteger, mustBeNonnegative}
        op.max_ite_ph1 (1, 1) {mustBeInteger, mustBeNonnegative} = 200
        op.max_ite_ph2 (1, 1) {mustBeInteger, mustBeNonnegative} = 100
        op.stp_slw_prg (1, 1) {mustBeInteger, mustBeNonnegative} = 10
        op.opt = sdpsettings('solver', 'mosek', 'verbose', 0);
        op.disp {mustBeA(op.disp, "logical")} = false
        op.eps (1, 1) {mustBeNumeric} = 1e-6
        op.Aeq = []
        op.beq = []
        op.A   = []
        op.b   = []
        op.con {mustBeA(op.con, "ss")} = ss([])
    end
    
    % Override display function of no display is desired
    if op.disp
        displ = @(x) disp(x);
    else
        displ = @(x) 0; 
    end

    if sys.Ts == 0 % Continuous time
        sblk = @(X) blkodiag(X, X);
    else % Discrete time
        sblk = @(X) blkdiag(X, -X);
    end
    
    % Skip first phase if a stabilizing controller is provided
    if ~isempty(op.con)
        op.max_ite_ph1 = 0; 
        Kv = op.con.d;
    end

    %% Abbreviations
    
    % *Some lengths*
    lx  = size(sys.a, 1);     % Number of system states
    err = size(sys, 1) - mea; % Size of error signal
    dis = size(sys, 2) - act; % Size of generalized disturbance
    inp = [dis, act];         % All input dimensions
    out = [err, mea];         % All output dimensions
    
    % *Outer factor*
    [A, B, C, D] = sssdata(sys, [dis, act], [err, mea]); % System data   
    
    OX = [A, B{1}, B{2}; eye(lx), zeros(lx, dis+act); ...
          C{1}, D{1,1}, D{1,2}; zeros(dis, lx), eye(dis), zeros(dis, act)];

    % Related dimension
    lox = size(OX, 2);
    
    %% Preparations for both phases
    
    % *All appearing variables*
    r  = sdpvar(1);  % Bound to feasibility / relaxation parameter
    ga = sdpvar(1);  % Upper bound on the energy gain
    X  = sdpvar(lx); % (Primal) Lyapunov certificate
    Y  = sdpvar(lx); % (Dual) Lyapunov certificate
    F  = sdpvar(act, lx+dis+act, 'full'); % Slack variable
    K  = sdpvar(act, mea, 'full');        % Controller
    
    % *Inner term and hermitian part*
    IX = blkdiag(sblk(X), eye(err), -ga*eye(dis));
    H  = F' * [K * C{2}, K * D{2, 1}, -eye(act)];

    % *Constraints for the problems with fixed controller*
    ConK = [OX' * IX * OX + H + H' <= (r-op.eps) * eye(lox); ...
            X >= op.eps * eye(lx)];
    % *Constraints for the problems with fixed slack variable*
    ConF = ConK;
    % General affine constraints on the controller matrices
    if ~isempty(op.A) && ~isempty(op.b)
        ConF = [ConF; op.A * K(:) <= op.b];
    end
    if ~isempty(op.Aeq) && ~isempty(op.beq)
        ConF = [ConF; op.Aeq * K(:) == op.beq];
    end

    % *Optimizers*
    YOPF1 = optimizer(ConF, r + 0.001*ga, op.opt, F, {r, ga, X, K});
    YOPK1 = optimizer(ConK, r + 0.001*ga, op.opt, K, {r, ga, X, F});
    YOPF  = optimizer(ConF, ga, op.opt, {r, F}, {ga, X, K});
    YOPK  = optimizer(ConK, ga, op.opt, {r, K}, {ga, X, F});

    %% Phase 1a (Build a full information gain for a suitable slack var.)

    displ('**Starting phase 1**');
  
    if isempty(op.con)
        % Dual outer factor
        OY = [eye(lx), zeros(lx, err); -A', -C{1}'; ...
              zeros(err, lx), eye(err); -B{1}', -D{1, 1}'] * ...
             null([B{2}', D{1, 2}']);
        % Splitting for direct optimization over the upper bound ga
        OY1 = OY(1:end-dis, :);
        OY2 = OY(end-dis+1:end, :);
        
        % Related dimension
        loy = size(OY, 2);
        
        % Inner term
        IY0 = blkdiag(sblk(Y), eye(err));
        
        % *Constraints*
        % Dual system LMI (after Schur)
        Con = [[OY1' * IY0 * OY1, OY2'; OY2, ga*eye(dis)] >= ...
                                                    op.eps * eye(loy+dis)];
        % Positivity of Lyapunov certificate
        Con = [Con; Y >= op.eps * eye(lx)];

        % *Solve full information design LMIs*
        optimize(Con, ga, op.opt);
    
        % *Reconstruct full information gain*
        Y = inv(chol(value(Y)));
        Y = Y * Y';
        
        % Involved middle matrix
        P = blkdiag(sblk(Y), eye(err), -value(ga)*eye(dis));
        P = middlematrix_permutation(P, dis, err, lx); 
        
        % Matrices W, U, V corresponding to the full information design 
        % problem
        [W, U, V] = WUV(ss2fiss(sys, inp, out), act, lx+sum(inp)-act);

        % Apply the elimination lemma
        Fv = elimi(P, U, V, W); % Get full information gain
        Fv = -[Fv, -eye(act)];  % Adjustments for a proper slack variable
    end
    
    %% Phase 1b (Trying to find a stabilizing controller)
        
    if isempty(op.con)
        % *Iteratively solve LMIs*
        ind = 1;
        while ind <= op.max_ite_ph1
            % Solve an LMI problem with fixed slack variable
            sol = YOPF1(Fv);
            rv  = sol{1};
            Kv  = sol{4};
    
            % Display progress
            displ(sprintf('Step: %i, Bound to feasibility: %0.5f',ind,rv));
            ind = ind + 1;
    
            % Stop if we have found a stabilizing controller
            if rv < 0; break; end
    
            % Solve an LMI problem with fixed controller matrices
            sol = YOPK1(Kv);
            rv  = sol{1};
            Fv  = sol{4};
    
            % Display progress
            displ(sprintf('Step: %i, Bound to feasibility: %0.5f',ind,rv));
            ind = ind + 1;
    
            % Stop if we have found a stabilizing controller
            if rv < 0; break; end
    
            if ind >= op.max_ite_ph1 % Then we were unsuccessfull
                error(['No stabilizing controller satisfying the ', ...
                       'given constraints was found. Try to weaken ', ...
                       'those or to increase max_ite_ph1']);
            end
        end
    end

    % *Compute energy gain and get a reasonable slack variable*   
    sol = YOPK({0, Kv});
    gav = sqrt(sol{1});
    Fv  = sol{3};

    displ(sprintf(['Found a stabilizing controller satisfying the ', ...
       'constraints and \nwith closed-loop energy-gain: ',...
       num2str(gav)]));


    %% Phase 2 (Minimize the closed-loop energy gain)

    displ('**Starting phase 2**');
    
    % *Iteratively solve LMIs*
    ind = 1;
    gav = NaN * ones(op.max_ite_ph2, 1);
    while ind <= op.max_ite_ph2-1
        % Solve design problem with fixed full information gain
        sol      = YOPF({0, Fv});
        gav(ind) = sqrt(sol{1});
        Kv       = sol{3};

        % Display progress
        displ(sprintf('gamma(%i) = %0.4f', ind, gav(ind)));
        ind = ind + 1;
        

        % Solve Design problem with fixed controller matrices
        sol      = YOPK({0, Kv});
        gav(ind) = sqrt(sol{1});
        Fv       = sol{3};

        % Display progress
        displ(sprintf('gamma(%i) = %0.4f', ind, gav(ind)));
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

    %% Design a corresponding static output-feedback controller

    % Final upper bound
    ga = gav(ind - 1);

    % Build final controller
    con = ss(Kv);
    
    % Perform a posteriori analysis
    displ(['The closed-loop energy gain is ', ...
           num2str(hinfnorm(lft(sys, con)))]);
end



