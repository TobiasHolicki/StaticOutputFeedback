%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : blkodiag.m                                                    %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 04                                                            %
% Date    : 18.03.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns for given objects A_1, ... A_N the object
% 
%      / 0      A_1 \
% A = |     ...      |
%      \A_N       0 /
% whenever blkdiag(A_1, ..., A_N) is well-defined.
%
% See also blkdiag, blklltriang, blklutriang, blkrltriang, blkrutriang
%
% ----- Input ---------------------------------------------------------- 
%   varargin - List of objects such as matrices, LTI systems, ...
%
% ----- Output ---------------------------------------------------------
%          A - An object as above
%
function [ A ] = blkodiag( varargin )

% Start with blkdiag (this is overloaded for anything useful)
A = blkdiag(varargin{:});

% Some abbreviations
N     = length(varargin); % Number of objects
idims = zeros(N, 1);      % Input dimensions

% Get relevant dimensions
for i = 1 : N
    idims(i) = size(varargin{i}, 2);
end
ln = sum(idims); % Total number of inputs

% Prepare (block) permutation
pn = mat2cell(1:ln, 1, idims);
pn = pn(N:-1:1);
pn = cell2mat(pn);

% Permute the block diagonal matrix to a block off diagonal matrix
A = A(:, pn);

end

