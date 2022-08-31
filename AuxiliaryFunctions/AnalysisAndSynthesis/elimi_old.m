%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : elimi.m                                                       %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 03                                                            %
% Date    : 18.03.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements the reconstruction corresponding to the 
% elimination lemma. The latter lemma states that for a given invertible 
% and symmetric P with exactly k negative eigenvalues and matrices U, V, W 
% there exists a solution Z of
% [eye(k); U' * Z * V + W]' * P * [eye(k); U' * Z * V + W] < 0
% if and only if
% V_' * [eye(k); W]' * P * [eye(k); W] * V_ < 0 and
% U_' * [W'; -eye(l)]' * P^(-1) * [W'; -eye(l)] * U_ > 0.
% Here V_ is a basismatrix of the kernel of V and U_ is a basismatrix of
% the kernel of U.
%
% ----- Input ---------------------------------------------------------- 
%   P, U, V, W - The matrices as above
% ----- Output ---------------------------------------------------------
%            Z - A solution of the LMI of interest
%
function [ Z ] = elimi(P, U, V, W)

% Get relevant dimension
k  = size(W, 2);

% Checking assumptions
V_ = null(V); % Annihilator
U_ = null(U); % Annihilator
OV = [eye(k); W] * V_;            % Primal outer factor
OU = [W'; -eye(size(W, 1))] * U_; % Dual outer factor

if max(eig(OV' * P * OV)) > 0
    warning('elimi::The primal LMI seems not to be satisfied');
elseif min(eig(OU' * (P \ OU))) < 0
    warning('elimi::The dual LMI seems not to be satisfied');
elseif find(sort(eig(P)) < 0, 1, 'last') ~= k
    warning(['elimi::P seems not to have the right number of ', ...
            'negative eigenvalues'])
end

% Transform U and V such that the annihilators U_ and V_ look nice
[Tu, su, Wu] = svd(U);
[Tv, sv, Wv] = svd(V);

% Stuff for later
l1 = rank(su);
l2 = size(W, 1) - l1;
k1 = rank(sv);
k2 = k - k1;
su = su(1:l1, 1:l1);
sv = sv(1:k1, 1:k1);

% Transformed W and transformed P
Wh = Wu' * (W * Wv);
Ph = blkdiag(Wv, Wu)' * P * blkdiag(Wv, Wu);

% Build R and S as in the proof
R = [[eye(k); Wh] * [eye(k1); zeros(k2, k1)], ...
     [zeros(k, l1); eye(l1); zeros(l2, l1)]];
S = [eye(k); Wh] * [zeros(k1, k2); eye(k2)];

% Build Pt. If the LMI has a solution there exists Zt with
% [eye(k1); Zt]' * Pt * [eye(k1); Zt] < 0. 
h1 = chol(-S' * Ph * S)';
h2 = h1 \ (S' * Ph * R);
Pt = R' * Ph * R + h2' * h2;

% Get eigenvalues and orthogonal eigenvectors of Pt and sort them beginning
% with the smalles eigenvalue
[u, t]  = schur(Pt);
[us, ~] = ordschur(u, t, 'lhp');
us      = us(:, 1:k1);

% Now we should have X(1:k1, :)' * Pt * X(1:k1, :) < 0 and can construct Zt
% as above. (Test: k1 == sum(ei < 0))
Z2 = us(k1 + 1:end, :);
Z1 = us(1:k1, :);

% Small disturbance for invertability...
 if cond(Z1) > 1e11 || cond(Z1) < 1e-11
     Z1 = disturb(Z1, 1e-10, 0.001);
 end
Zt = Z2 / Z1;

% Transform back
Z11 = su \ (Zt / sv);

% This choice is arbitrary
Zh = zeros(size(U, 1), size(V, 1));

% Relevant part of the transformed solution
Zh(1:l1, 1:k1) = Z11;

% Transform back
Z = Tu * Zh * Tv';

% Final check
O = [eye(size(W, 2)); W + U' * Z * V];
e = max(real(eig(O' * P * O)))
if e > 0
   warning(['elimi::The reconstruction failed somehow. Maximum ', ...
            'eigenvalue in the desired LMI is ', num2str(e)]) 
end

end

