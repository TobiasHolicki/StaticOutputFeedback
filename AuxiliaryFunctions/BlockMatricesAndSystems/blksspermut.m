%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : blksspermut.m                                                 %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 03                                                            %
% Date    : 18.03.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Let G be an LTI system with state-space realization [A, B, C, D] where
% A is partitioned according bdim. This function returns the same LTI 
% system G where the states are reordered according to the permutation
% bper.
% Suppose for example that
%       A = [A11, A12, A13;
%            A21, A22, A23;
%            A31, A32, A33]
% with Akl of size bdim(k) times bdim(k) and that 
%       bper = [3, 1, 2].
% Then the output of this function is the system G with G.a
%           [A33, A31, A32;
%            A13, A11, A12;
%            A32, A21, A22]
% and G.b as well as G.c adjusted accordingly.
%
% ----- Input ---------------------------------------------------------- 
%   G    - Given LTI system
%   bdim - Vector of state dimensions
%   bper - Block permutation for the states
%
% ----- Output ---------------------------------------------------------
%   G    - The LTI system after state block permutations
%  
function [G] = blksspermut(G, bdim, bper)

% For sanity checks
la = size(G.a, 1); % Total number of states

% Some preparations
ls = sum(bdim); % Total number of states

% Some sanity checks
if la ~= ls
    error('blksspermut::Number of states and partitions do not match')
elseif length(bdim) ~= length(bper)
    error(['blksspermut::Length of state permutation and block ', ...
           'partition do not match'])
elseif min(bper) <= 0 || max(bper) > length(bdim)
    error('blksspermut::State permutation is not a valid permutation')
elseif length(unique(bper)) ~= length(bper)
    error('blksspermut::State permutation is not a valid permutation')
end
       
bs = mat2cell(1:ls, 1, bdim); % Indices of state partition 

bs = bs(bper); % Permute state blocks

% Permutation on original indices 
bs = cell2mat(bs);

% Permute the states of the given system G
A = G.a(bs, bs);
B = G.b(bs, :);
C = G.c(:, bs);
D = G.d;

G = ss(A, B, C, D);

end

