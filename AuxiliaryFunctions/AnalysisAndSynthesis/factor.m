%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : factor.m                                                      %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 02                                                            %
% Date    : 02.07.2022                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
function [ fact ] = factor( sys, opt )
    arguments
        sys {mustBeA(sys, "ss")}
        opt.primaldual {mustBeMember(opt.primaldual,["primal","dual"])} ...
                        = 'primal'
    end
    
    % Get system matrices
    [a, b, c, d] = ssdata(sys);

    % Build outer factor
    if strcmp(opt.primaldual, 'primal')
        fact = [a, b; eye(size([a, b])); c, d];
    else
        fact = [eye(size([a', c'])); -a', -c'; -b', -d'];
    end
end