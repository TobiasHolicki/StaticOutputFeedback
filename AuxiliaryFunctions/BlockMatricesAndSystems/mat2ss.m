%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : mat2ss.m                                                      %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 01                                                            %
% Date    : 06.10.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a given partitioned matrix M = [A, B; C, D] with A being of dimension
% n x n, this function returns the system with realization [A, B, C, D].
%
function [sys] = mat2ss(M, n, varargin)

    A = M(    1:n,     1:n);
    B = M(    1:n, n+1:end);
    C = M(n+1:end,     1:n);
    D = M(n+1:end, n+1:end);

    lv = length(varargin);
    if lv == 0
        sys = ss(A, B, C, D);
    elseif lv == 1
        sys = ss(A, B, C, D, varargin{1});
    else
        error('Wrong number of input arguments')
    end

end

