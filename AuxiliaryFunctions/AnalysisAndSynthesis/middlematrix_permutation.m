%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : middlematrix_permutation.m                                    %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 03                                                            %
% Date    : 06.10.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a block diagon matrix P as appearing in the standard synthesis LMIs,
% this function computes a permutation such that we can apply the 
% elimination lemma afterwards.
%
% ----- Input ---------------------------------------------------------- 
%         P - Middle matrix
%       inp - A vector containing the dimensions of size(Bi, 2) for all i
%       out - A vector containing the dimensions of size(Ci, 1) for all i
%  statedim - This should equal size(A, 1)
% ----- Output ---------------------------------------------------------
%         P - Permuted matrix
%
function [P] = middlematrix_permutation(P, inp, out, statedim)

% Incorporate state channel
inp = [statedim, inp(:)'];
out = [statedim, out(:)'];

% Some Abbreviations
n   = length(inp); % number of channels (should be the same as length(out))
io  = [inp, out];  % input and outpud dimensions
lio = sum(io);     % total size of inputs and outputs

% Initialize the permutation
p = 1 : lio;

% Initialize desired block permutation (left one)
pn = zeros(1, 2*n);

% Construct desired block permutation (left one)
for i = 1 : n
    pn(2*i-1) = i+n;
    pn(2*i)   = i;
end

% Construct desired permutation (right one)
pcr = mat2cell(p, 1, io(pn));
pnr = cell2mat(pcr([2:2:2*n, 1:2:2*n]));

P = P(pnr, pnr);

end

