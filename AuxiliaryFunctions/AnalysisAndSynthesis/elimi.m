%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : elimi.m                                                       %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 05                                                            %
% Date    : 19.07.2022                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements the reconstruction corresponding to the two
% elimination lemmas from [1]. 
%
% The more general one states that, for a given symmetric lxl matrix P and 
% matrices U in R^{mxl}, V in R^{nxk}, W in R^{lxk}, there exists a 
% solution Z in R^{mxn} of the inequality
% 			[W + U' * Z * V]' * P * [W + U' * Z * V] < 0
% if and only if
% 			V_' * W' * P * W * V_ < 0  
% holds and 
%  			(W U') * P * (W U')  has at least k negative eigenvalues.
% Here V_ is a basismatrix of the kernel of V.
%
% The more common one states that, for a given symmetrix (k+l)x(k+l) matrix
% P with exactly k negative and l positive eigenvalues and for matrices
% U in R^{mxl}, V in R^{nxk}, W in R^{lxk}, there exists a solution 
% Z in R^{mxn} of the inequality
%         [I_k; W + U' * Z * V]' * P * [I_k; W + U' * Z * V] < 0
% if and only if
%               V_' * [I_k; W]' * P * [I_k; W] * V_ < 0
% and
%        U_' * [-W'; I_l]' * P^{-1} * [-W'; I_l] * U_' > 0
% hold, where U_ and V_ are basismatrices of the kernels of U and V,
% respectively.
%
% [1] A. Helmersson, IQC synthesis based on inertia constraints, 1999
%
% ----- Input ---------------------------------------------------------- 
%   P, U, V, W - The matrices as above
% ----- Output ---------------------------------------------------------
%            Z - A solution of the LMI of interest
function [ Z ] = elimi(P, U, V, W)
     % Some sanity checks
    if size(W, 1) ~= size(U', 1)
        error('elimi::U^T and W must have the same number of rows')
    elseif size(W, 2) ~= size(V , 2)
        error('elimi::V and W must have the same number of columns')
    elseif norm(P - P', 2) > 1e-10
        error('elimi::P must be symmetric')
    end
    
    % Try to find the variant that is supposed to be used
    [l, k] = size(W);
    if size(P, 1) == k+l % Second variant
        % We simply redefine W and U
        W = [eye(k); W];
        U = [zeros(size(U, 1), k), U];
    elseif size(P, 1) ~= l % 
        error('elimi::Please check the dimensions of P and W')
    else % i.e., size(P, 1) == l
        % Then we do not have to do anything
    end 
    
    % Check the assumptions
    V_ = null(V); % Annihilator 
    
    if max(eig((W * V_)' * P * (W * V_))) > 0
        warning('elimi::The primal LMI seems not to be satisfied');
    elseif find(sort(eig([W, U']' * P * [W, U'])) < 0, 1, 'last') < k
        warning(['elimi:: [W, U^T]^T P [W U^T] seems not to have the '...
                 'right number of negative eigenvalues'])
    end

    % Transform V such that the annihilator V_ looks nice
    [Tv, sv, Wv] = svd(V);

    % Stuff for later
    k1 = rank(sv);
    sv = sv(1:k1, 1:k1);
    
    % Transformed W
    Wh = W * Wv;
    
    % Build R and S as in the proof
    R = [Wh(:, 1:k1), U'];
    S = Wh(:, k1+1:end);

    % Build Pt. If the LMI has a solution there exists Zt with
    % [eye(n1); Zt]' * Pt * [eye(n1); Zt] < 0. 
    h1 = chol(-S' * P * S)';
    h2 = h1 \ (S' * P * R);
    Pt = R' * P * R + h2' * h2;

    % Get eigenvalues and orthogonal eigenvectors of Pt and sort them 
    % beginning with the smalles eigenvalue.
    [u, t]  = schur(Pt);
    [us, ~] = ordschur(u, t, 'lhp');
    us      = us(:, 1:k1);

    % Now we should have us' * Pt * us < 0 and can construct Zt as above.
    Z2 = us(k1 + 1:end, :);
    Z1 = us(1:k1, :);

    % Small disturbance for invertability...
     if cond(Z1) > 1e11 || cond(Z1) < 1e-11
         Z1 = disturb(Z1, 1e-10, 0.001);
     end
    Zt = Z2 / Z1;

    % Transform back
    Z11 = Zt / sv;
    
    % This zeros in the following can be replaced by anything else
    Zh = [Z11, zeros(size(U, 1), size(V, 1) - k1)];

    % Transform back
    Z = Zh * Tv';

    % Final check
    O = W + U' * Z * V;
    e = max(eig(O' * P * O));
    if e > 0
       warning(['elimi::The reconstruction failed somehow. Maximum ', ...
                'eigenvalue in the desired LMI is ', num2str(e)]) 
    end

end

