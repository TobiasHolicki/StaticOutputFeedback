%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : sssdata.m                                                     %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 04                                                            %
% Date    : 02.07.2022                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a given system G as well as row and column partitions out, inp this
% functions returns the partitioned system matrices A, B, C, D. This is
% done by using mat2cell. 
%
% ----- Input ---------------------------------------------------------- 
%   G                    - System in ss representation
%   inp = [c1, ..., cm]  - Input/column partition
%   out = [r1, ..., rn]  - Output/row parition
% ----- Output ---------------------------------------------------------
%   A, B, C, D           - Partitioned matrices as cell, 
%                          e.g. D{i,j} is a ci times rj matrix
%
function [A, B, C, D] = sssdata(G, inp, out)
    % Some sanity checks
    arguments
        G {mustBeA(G, "ss")}
        inp (1, :) {mustBeInteger, mustBeNonnegative}
        out (1, :) {mustBeInteger, mustBeNonnegative}
    end
    if size(G, 2) ~= sum(inp)
      error('sssdata::Number of inputs and input partition do not match')
    elseif size(G, 1) ~= sum(out)
      error('sssdata::Number of outputs and output partition do not match')
    end
    
    % Number of states
    la = size(G.a, 1);

    % Partitioned system matrices 
    A = G.a;
    B = mat2cell(G.b, la, inp);
    C = mat2cell(G.c, out, la);
    D = mat2cell(G.d, out, inp);
end

