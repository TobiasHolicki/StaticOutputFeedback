%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : zez.m                                                         %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 01                                                            %
% Date    : 23.06.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For numbers k, l, m this functions returns the matrix
% [0_{k x l};
%  I_{l x l};
%  0_{m x l}].
%
% ----- Input ---------------------------------------------------------- 
%  k, l, m - The above numbers. If m is omitted it is set to zero.
% ----- Output ---------------------------------------------------------
%       M - The matrix as above.
%
function [M] = zez(k, l, m)
    % Some sanity checks
    arguments
        k (1, 1) {mustBeInteger, mustBeNonnegative}
        l (1, 1) {mustBeInteger, mustBeNonnegative}
        m (1, 1) {mustBeInteger, mustBeNonnegative} = 0
    end

    M = [zeros(k, l); eye(l); zeros(m, l)];
end

