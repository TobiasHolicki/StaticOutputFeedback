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
        op.opt = sdpsettings('solver', 'mosek', 'verbose', 0);
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
    t = optimize(Con, Cos);

    X = value(X);
    Y = value(Y);
    ga = value(ga);
    
    
    X = [  528.7977   14.0677   10.8903  -82.2705   25.4280   64.2921   20.0852
   14.0677   29.0870   31.8307  -58.6602    7.5963   48.9625   19.9646
   10.8903   31.8307  136.4501    8.5842   18.3011  102.3715   78.8410
  -82.2705  -58.6602    8.5842  777.9732  -88.6888 -222.9135    9.5882
   25.4280    7.5963   18.3011  -88.6888  529.1399   56.5574   19.4665
   64.2921   48.9625  102.3715 -222.9135   56.5574  560.0435  231.0430
   20.0852   19.9646   78.8410    9.5882   19.4665  231.0430  138.6695];
    Y = [    0.9553    0.0021   -1.0203   -0.0205   -0.0012    0.2411    0.1627
    0.0021    0.6392   -0.9921   -0.3176    0.0854   16.3288    4.8086
   -1.0203   -0.9921    4.8649    1.4795    0.0494   -4.9662   12.3043
   -0.0205   -0.3176    1.4795    0.9317   -0.0116   -1.4140   -1.4154
   -0.0012    0.0854    0.0494   -0.0116    6.3866    5.2009    3.3501
    0.2411   16.3288   -4.9662   -1.4140    5.2009  842.5047  202.9977
    0.1627    4.8086   12.3043   -1.4154    3.3501  202.9977  737.5330];
    ga = 11.8455;

U = inv(Y) - X;
    Xe = [X, U; U, -U];

    % Middle matrix
    P = blkdiag(blkodiag(Xe, Xe), eye(err), -ga*eye(dis));
    P = middlematrix_permutation(P, dis, err, 2*lx);
    
    % Preparations for applying elimination
    [W, U, V] = WUV(sys, act, mea, lx);
    con         = mat2ss(elimi(P, U, V, W), lx);

   
    
con = 0;
ga = 0;
end



