%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : outerfactor_permutation.m                                     %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 03                                                            %
% Date    : 18.03.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes a permutation pnl such that for
% OU = [eye(...);
%       A, B1, ... BN; 
%       C1, D11, ... D1N;
%       :    :         :
%       CN, DN1, ... DNN];
% the output of OU(pnl, :) equals
% [A, B1, ..., BN;
%  I,  0, ....,  0;
%  C1, D11, ..., D1N;
%  0,  I,  ...,  0;
%   :   :    :    :
%  CN, DN1, ..., DNN;
%   0,  0, ....,  I].
% The permutation pnl corresponds to a multiplication with a permutation
% matrix Q from the left. This function also computes the permutation pnr
% that corresponds to the multiplication of Q from the right.
%
% ----- Input ---------------------------------------------------------- 
%   inp     - A vector containing the dimensions of size(Bi, 2) for all i
%   out     - A vector containing the dimensions of size(Ci, 1) for all i
%  statedim - This should equal size(A, 1)
% ----- Output ---------------------------------------------------------
%   pnl     - The permutation described above
%   pnr     - The permutation described above
%
function [pnl, pnr] = outerfactor_permutation(inp, out, statedim)

% Incorporate state channel
inp = [statedim, inp(1, :)];
out = [statedim, out(1, :)];

% Some Abbreviations
n   = length(inp); % number of channels (should be the same as length(out))
io  = [inp, out];  % input and outpud dimensions
lio = sum(io);     % total size of inputs and outputs

% Initialize the permutation
p = 1 : lio;

% transform to cell
pcl = mat2cell(p, 1, io);

% Initialize desired block permutation (left one)
pn = zeros(1, 2*n);

% Construct desired block permutation (left one)
for i = 1 : n
    pn(2*i-1) = i+n;
    pn(2*i)   = i;
end

% Construct desired permutation (left one)
pnl = cell2mat(pcl(pn));

% Construct desired permutation (right one)
pcr = mat2cell(p, 1, io(pn));
pnr = cell2mat(pcr([2:2:2*n, 1:2:2*n]));


end

