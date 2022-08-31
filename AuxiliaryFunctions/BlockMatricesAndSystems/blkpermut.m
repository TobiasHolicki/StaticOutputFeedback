%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : blkpermut.m                                                   %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 01                                                            %
% Date    : 12.03.2019                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a given matrix or LTI system G with block partition according to inp
% and out this functions returns a matrix or LTI where the blocks are
% permuted according to the input permutations bpi and output permutations
% bpo.
% Suppose for example that
%       G = [G11, G12, G13, G14;
%            G21, G22, G23, G24;
%            G31, G32, G33, G34]
% with Gkl of size out(k) times inp(l) and that 
%     bpi = [3, 1, 2, 4] as well as bpo = [3, 2, 1].
% Then the output of this function is
%           [G33, G31, G32, G34;
%            G23, G21, G22, G24;
%            G13, G11, G12, G14].
% Inputs can be abbreviated with ':'. This is for example useful if only 
% rows need to be permuted while the columns stay fixed.  
%
% ----- Input ---------------------------------------------------------- 
%   G   - Given matrix or LTI system
%   inp - Vector of input dimensions
%   out - Vector of output dimensions
%   bpi - Block permutation for inputs
%   bpo - Block permutation for outputs
%
% ----- Output ---------------------------------------------------------
%   G   - The fiven matrix or LTI system after block permutations
%  
function [G] = blkpermut(G, inp, out, bpi, bpo)

% Dimensions of the given G
[ou, in] = size(G);

% Handle some abbreviations
if inp == ':'
    inp = in;
    bpi = 1; 
elseif out == ':'
    out = ou;
    bpo = 1;
end
if bpi == ':'
    bpi = 1 : length(inp);
elseif bpo == ':'
    bpo = 1 : length(out);
end

% Some preparations
lin = sum(inp); % Total number of inputs
lou = sum(out); % Total number of outputs

% Some sanity checks
if ou ~= lou
    error('blkpermut::Output dimensions and partitions do not match')
elseif in ~= lin
    error('blkpermut::Input dimensions and partitions do not match')
elseif length(inp) ~= length(bpi)
    error(['blkpermut::Length of input permutation and input block ', ...
           'partition do not match'])
elseif length(out) ~= length(bpo)
    error(['blkpermut::Length of output permutation and output block ', ...
           'partition do not match'])
elseif min(bpi) <= 0 || max(bpi) > length(inp)
    error('blkpermut::Input permutation is not a valid permutation')
elseif min(bpo) <= 0 || max(bpo) > length(out)
    error('blkpermut::Output permutation is not a valid permutation')
end
       
bin = mat2cell(1:lin, 1, inp); % Indices of input partition in cells
bou = mat2cell(1:lou, 1, out); % Indices of output partition in cells

bin = bin(bpi); % Permute input blocks
bou = bou(bpo); % Permute output blocks

% Get permuted indices
bin = cell2mat(bin);
bou = cell2mat(bou);

% Permute the given G
G = G(bou, bin);

end

