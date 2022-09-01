%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : primal_dual_iteration.m                                       %
%                                                                         %
% Author  : Tobias Holicki                                                %
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
% The design procedure is based on the dual iteration proposed in [1] and
% which as been described in detail in the paper [2]. The phase of finding
% a stabilizing gain is essentially from [3].
% I took the freedom to denote this algorithm as primal dual iteration
% since, indeed, both primal and dual controller design problems are
% solved and to contrast this with other algorithms that only solve primal
% problems.
%
% [1] Iwasaki, T., The dual iteration for fixed-order control, 1999.
% [2] Holicki, T. and Scherer, C. W., Revisiting the dual iteration for 
%     static and robust output-feedback, 2021.
% [3] El Ghaoui, L. and Oustry, F. and AitRami, M., A cone complementarity
%     linearization algorithm for static output-feedback and related
%     problems, 1997.
%
% See also pk_iteration, svar_iteration, fo_iteration, primal_iteration
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
%   con         - Some a priori known stabilizing controller.
% ----- Output ---------------------------------------------------------
%   con         - Computed controller
%   ga          - Computed upper bound on the closed-loop energy gain.
%  
function [con, ga] = primal_dual_iteration(sys, mea, act, op)
    % Some sanity checks
    arguments
        sys {mustBeA(sys, "ss")}
        mea (1, 1) {mustBeInteger, mustBeNonnegative}
        act (1, 1) {mustBeInteger, mustBeNonnegative}
        op.max_ite_ph1 (1, 1) {mustBeInteger, mustBeNonnegative} = 10
        op.max_ite_ph2 (1, 1) {mustBeInteger, mustBeNonnegative} = 20
        op.stp_slw_prg (1, 1) {mustBeInteger, mustBeNonnegative} = 5
        op.opt = sdpsettings('solver', 'sdpt3', 'verbose', 0);
        op.disp {mustBeA(op.disp, "logical")} = false
        op.eps (1, 1) {mustBeNumeric} = 1e-6
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
    
    %% Abbreviations
    
    % *Some lengths*
    lx  = size(sys.a, 1);     % Number of system states
    err = size(sys, 1) - mea; % Size of error signal
    dis = size(sys, 2) - act; % Size of generalized disturbance
    inp = [dis, act];         % All input dimensions
    out = [err, mea];         % All output dimensions
    
    % *Outer factors*
    [A, B, C, D] = sssdata(sys, [dis, act], [err, mea]); % System data
    OX = [A, B{1}; eye(lx), zeros(lx, dis); ...
          C{1}, D{1, 1}; zeros(dis, lx), eye(dis)];
    OY = [eye(lx), zeros(lx, err); -A', -C{1}'; ...
          zeros(err, lx), eye(err); -B{1}', -D{1, 1}'];
    % Factors involving the full information/actuation gains
    OXF = @(F) OX + [B{2}; zeros(lx, act); D{1, 2}; zeros(dis, act)] * F;
    OYE = @(E) OY - [zeros(mea, lx), C{2}, zeros(mea, err), D{2, 1}]' * E';
    % Factors with annihilators
    OX = OX * null([C{2}, D{2, 1}]);
    OY = OY * null([B{2}', D{1, 2}']);
    
    % Splitting for direct optimization over the upper bound ga
    OY1 = OY(1:end-dis, :);
    OY2 = OY(end-dis+1:end, :);
    
    % Related dimensions
    lox = size(OX, 2);
    loy = size(OY, 2);
    
    %% Preparations for both phases
    
    % *All appearing variables*
    r  = sdpvar(1);  % Bound to feasibility / relaxation parameter
    ga = sdpvar(1);  % Upper bound on the energy gain
    X  = sdpvar(lx); % (Primal) Lyapunov certificate
    Y  = sdpvar(lx); % (Dual) Lyapunov certificate
    X0 = sdpvar(lx); % Iterated primal Lyapunov certificate for phase 1
    Y0 = sdpvar(lx); % Iterated dual Lyapunof certificate for phase 1
    F  = sdpvar(act, lx+dis, 'full'); % Full information gain
    E  = sdpvar(lx+err, mea, 'full'); % Full actuation gain
    
    % *Inner terms*
    IX  = blkdiag(sblk(X), eye(err), -ga*eye(dis));
    IY  = blkdiag(sblk(Y), ga*eye(err), -eye(dis));
    IY0 = blkdiag(sblk(Y), eye(err));
    
    % *Constraints*
    % Primal problem involving the full information gain F
    ConF = [OX' * IX * OX <= -op.eps * eye(lox); ...
            OXF(F)' * IX * OXF(F) <= -op.eps * eye(lx+dis)];
    % Dual problem involving the full actuation gain E
    ConE = [OY' * IY * OY >= op.eps * eye(loy); ...
            OYE(E)' * IY * OYE(E) >= op.eps * eye(lx+err)];
    % Dynamic output feedback design problem appearing in phase 1
    Con = [OX' * IX * OX <= - op.eps * eye(lox); ... % Primal
           [OY1' * IY0 * OY1, OY2'; OY2, ga*eye(dis)] >= ...
                            op.eps * eye(loy+dis); ... % Dual (after Schur)
           [X, eye(lx); eye(lx), Y] >= op.eps * eye(2*lx)];
    
    % *Cost*
    % There is a small penalty on ga to avoid starting the second phase 
    % with huge bounds on the energy gain. 
    Cos =  trace(X * Y0 + X0 * Y) + 0.001*ga;

    % *Optimizers*
    YOPF = optimizer(ConF, ga, op.opt, F, {ga, X});
    YOPE = optimizer(ConE, ga, op.opt, E, {ga, Y});
    YOP0 = optimizer(Con, Cos, op.opt, {X0, Y0}, {ga, X, Y});


    %% Phase 1a (Trying to find a stabilizing controller)
    % This is a heuristic that often works quite well and begins with the
    % LMIs for dynamic output feedback design. 
    
    if isempty(op.con)
        displ('**Starting phase 1**');
        
        % Initial certificates
        X0i = eye(lx);
        Y0i = eye(lx);

        % *Iteratively solve LMIs*
        for i = 1 : op.max_ite_ph1
            % *Solve dynamic output-feedback design problem*
            sol = YOP0({X0i, Y0i});

            % Display progress (by omitting ga in the cost, res should be
            % monotonically decreasing)
            res = trace(X0i * sol{3} + Y0i * sol{2});
            displ(sprintf(['Push heuristic : tr(X_%i*Y_%i + Y_%i*X%i)', ...
                                        ' = %0.4f'], i-1, i, i-1, i, res));
            % Update iterates
            X0i = sol{2};
            Y0i = sol{3};

            % Potentially stop earlier (res is bounded from below by 2*lx)
            if res <= 2*lx + 1e-1
                break;
            end
        end

        % *Build a corresponding stabilizing full actuation gain*
        ga = sol{1};
        Ei = build_fa_gain(sqrt(ga), X0i);
    end   
    
    %% Phase 1b (Still trying to find a stabilizing controller)
    % X might not be a suitable Lyapunov matrix for continuing since its
    % inverse might not be suitable enough. We consider another problem
    % with some relaxation parameter r.
    
    if isempty(op.con)
        % *Inner term and constraints*
        IY  = blkdiag(sblk(Y), value(ga)*eye(err), -eye(dis));
        Con = [OY' * IY * OY >= (op.eps-r) * eye(loy); ...
               OYE(Ei)' * IY * OYE(Ei) >= op.eps * eye(lx+err)];

        % *Solve the dual design problem*
        t = optimize(Con, r, op.opt);
        
        if t.problem ~= 0
            warning(['Numerical problems occured in the initial dual ', ...
                     'design problem']);
        end
        if value(r) > 0 
            % In this case one could also try to iteratively reduce r, but
            % this is not implemented here.
            error(['The heuristic for finding a stabilizing controller',...
                   'failed. Please provide a stabilizing static gain.']);
        else
            displ(['Found a static stabilizing controller achieving ', ...
                   'the closed-loop energy-gain: ', ...
                   num2str(sqrt(value(ga)))]);

            % *Build a corresponding full information gain*
            Fi = build_fi_gain(sqrt(value(ga)), value(Y));
        end
    else
        displ(['The initially provided static stabilizing controller ', ...
               'achieves the closed-loop energy gain: ', ...
               num2str(hinfnorm(lft(sys, op.con)))]);
        
        % Get corresponding full information gain
        Fi = op.con.d * [C{2}, D{2, 1}];  
    end     

    
    %% Phase 2 (Minimize the closed-loop energy gain)

    displ('**Starting phase 2**');
    
    % *Iteratively solve LMIs*
    ind = 1;
    gav = NaN * ones(op.max_ite_ph2, 1);
    while ind <= op.max_ite_ph2-1
        % Solve the primal design problem for a given full information gain
        sol      = YOPF(Fi);
        gav(ind) = sqrt(sol{1});

        % Display progress
        displ(sprintf('gamma(%i) = %0.4f', ind, gav(ind)));
        ind = ind + 1;
        
        % Build full actuation gain
        Ei = build_fa_gain(gav(ind-1), sol{2});

        % Solve the dual design problem for a given full actuation gain
        sol      = YOPE(Ei);
        gav(ind) = sqrt(sol{1});
        
        % Display progress
        displ(sprintf('gamma(%i) = %0.4f', ind, gav(ind)));
        ind = ind + 1;
        
        % Build full information gain
        Fi = build_fi_gain(gav(ind-1), sol{2});
        
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

    % Final design
    sol = YOPF(Fi);
    gav = sol{1};
    Xv  = sol{2};

    % Middle matrix
    P = blkdiag(sblk(Xv), eye(err), -gav*eye(dis));
    P = middlematrix_permutation(P, dis, err, lx); 

    % U, V, W...
    [W, U, V] = WUV(sys, act, mea);

    % Build final controller by the elimination lemma
    con = ss(elimi(P, U, V, W));

    % Final upper bound
    ga = sqrt(gav);
    
    % Perform a posteriori analysis
    displ(['The closed-loop energy gain is ', ...
           num2str(hinfnorm(lft(sys, con)))]);

    
    %% Auxilliary functions
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function: build_fa_gain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function constructs the full actuation gain corresponding to the
    % upper bound ga and the Lyapunov certificate X based on elimination.
    function E = build_fa_gain(ga, X)
        % Involved middle matrix
        P = blkdiag(sblk(X), eye(err), -ga^2 * eye(dis));
        P = middlematrix_permutation(P, dis, err, lx); 
        
        % Matrices W,U,V corresponding to the full actuation design problem
        [W, U, V] = WUV(ss2fass(sys, inp, out), lx+sum(out)-mea, mea);
  
        % Apply the elimination lemma
        E = elimi(P, U, V, W); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function: build_fi_gain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function constructs the full information gain corresponding to
    % the upper bound ga and the Lyapunov certificate Y based on 
    % elimination.
    function F = build_fi_gain(ga, Y)
        % Y is a dual certificate and it is somewhat mroe convenient to
        % work with the primal
        Y = inv(chol(Y));
        Y = Y * Y';
        
        % Involved middle matrix
        P = blkdiag(sblk(Y), eye(err)/ga^2, -eye(dis));
        P = middlematrix_permutation(P, dis, err, lx); 
        
        % Matrices W,U,V corresponding to the full information design 
        % problem
        [W, U, V] = WUV(ss2fiss(sys, inp, out), act, lx+sum(inp)-act);

        % Apply the elimination lemma
        F = elimi(P, U, V, W);
    end
end



