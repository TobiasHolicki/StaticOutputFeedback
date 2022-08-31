%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : ss2fiss.m                                                     %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 01                                                            %
% Date    : 06.10.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function changes the state-space description of a system sys such
% that it has a structure appearing in a full information controller design
% problem.
%
% ----- Input ---------------------------------------------------------- 
%  sys - The given system as ss object
%  inp - Vector of input dimensions
%  out - Vector of output dimensions
% ----- Output ---------------------------------------------------------
%  sys - Modified system
%
function [sys] = ss2fiss(sys, inp, out)

% Abbreviations
mea = out(end);       % Number of measurements
act = inp(end);       % Number of actuators
la  = size(sys.a, 1); % Number of states

% Get system data
[A, B, C, D] = ssdata(sys);

% Adjust C and D together
CD = [C(1:end-mea, :), D(1:end-mea, :); ...
      eye(la+sum(inp)-act, la+sum(inp))];

sys = ss(A, B, CD(:, 1:la), CD(:, la+1:end));

end

