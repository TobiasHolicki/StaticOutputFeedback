%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : dual_iteration.m                                              %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 04                                                            %
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
% The design procedure is based on the dual iteration proposed in [1] and
% which as been described in detail in the paper [2]. The phase of finding
% a stabilizing gain is essentially from [3].
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
function [con, ga] = dual_iteration(sys, mea, act, op)
    % Some sanity checks
    arguments
        sys {mustBeA(sys, "ss")}
        mea (1, 1) {mustBeInteger, mustBeNonnegative}
        act (1, 1) {mustBeInteger, mustBeNonnegative}
        op.max_ite_ph1 (1, 1) {mustBeInteger, mustBeNonnegative} = 10
        op.max_ite_ph2 (1, 1) {mustBeInteger, mustBeNonnegative} = 20
        op.stp_slw_prg (1, 1) {mustBeInteger, mustBeNonnegative} = 5
        op.opt = sdpsettings('solver', 'sdpt3', 'verbose', 1);
        op.disp {mustBeA(op.disp, "logical")} = false
        op.eps (1, 1) {mustBeNumeric} = 1e-4
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
    ga = sdpvar(1);  % Upper bound on the energy gain
    X  = sdpvar(lx); % (Primal) Lyapunov certificate
    Y  = sdpvar(lx); % (Dual) Lyapunov certificate

    % *Inner terms*
    IX  = blkdiag(sblk(X), eye(err), -ga*eye(dis));
    IY0 = blkdiag(sblk(Y), eye(err));
    
    % *Constraints*
    % Dynamic output feedback design problem appearing in phase 1
    Con = [OX' * IX * OX <= - op.eps * eye(lox); ... % Primal
           [OY1' * IY0 * OY1, OY2'; OY2, ga*eye(dis)] >= ...
                            op.eps * eye(loy+dis); ... % Dual (after Schur)
           [X, eye(lx); eye(lx), Y] >= op.eps * eye(2*lx);
           X <= 1e3 * eye(lx); Y <= 1e3 * eye(lx)];
    
    % *Cost*
    % There is a small penalty on ga to avoid starting the second phase 
    % with huge bounds on the energy gain. 
    Cos =  ga;

    % *Optimizers*
    t = optimize(Con, Cos, op.opt)

    


    X = value(X);
    Y = value(Y);
    ga = value(ga);
    
    
 

    U = inv(Y) - X;
    Xe = [X, U; U, -U];

    % Middle matrix
    P = blkdiag(blkodiag(Xe, Xe), eye(err), -ga*eye(dis));
    P = middlematrix_permutation(P, dis, err, 2*lx);
    
    % Preparations for applying elimination
    [W, U, V] = WUV(sys, act, mea, lx);
    con         = mat2ss(elimi(P, U, V, W), lx);
    ga = sqrt(ga);
end



