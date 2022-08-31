%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : ss2fass.m                                                     %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 01                                                            %
% Date    : 06.10.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function changes the state-space description of a system sys such
% that it has a structure appearing in a full actuation controller design
% problem.
%
% ----- Input ---------------------------------------------------------- 
%  sys - The given system as ss object
%  inp - Vector of input dimensions
%  out - Vector of output dimensions
% ----- Output ---------------------------------------------------------
%  sys - Modified system
%
function [sys] = ss2fass(sys, inp, out)

% Abbreviations
mea = out(end);       % Number of measurements
act = inp(end);       % Number of actuators
la  = size(sys.a, 1); % Number of states

% Get system data
[A, B, C, D] = ssdata(sys);

% Adjust B and D together
BD = [[B(:, 1:end-act); D(:, 1:end-act)], ...
      eye(la+sum(out)-mea, la+sum(out))'];

sys = ss(A, BD(1:la, :), C, BD(la+1:end, :));

end

