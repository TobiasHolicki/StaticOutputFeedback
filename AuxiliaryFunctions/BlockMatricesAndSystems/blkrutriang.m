%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : blkrutriang.m                                                 %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 04                                                            %
% Date    : 30.07.2022                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns for given matrices or systems A_1, ..., A_anz a 
% block right upper triangular matrix or system with blocks A_1, ...,A_anz.
% Therefore there must exist some N with anz = N * (N + 1) / 2. 
% As an example the output for anz = 6 is
%      / A_1 A_2 A_3 \
% A = |   0  A_4 A_5  |.
%      \  0   0  A_6 /
% Moreover, if A_i is not on the diagonal of the resulting matrix, 
% A_i = 'z' results in a zero matrix of correct size and 
% A_i = 'e' results in an identity matrix of correct size.
%
% See also blkdiag, blkodiag, blklutriang, blkrltriang, blklltriang
%
% ----- Input ---------------------------------------------------------- 
%   varargin - List of objects such as matrices, LTI systems, ...
% ----- Output ---------------------------------------------------------
%          A - Object as described above
%
function [ A ] = blkrutriang( varargin )

anz = length(varargin);            % Number of inputs
n   = -0.5 + sqrt(0.25 + 2 * anz); % Number of block columns

% Sanity check
if n - floor(n) ~= 0
    error('blkrutriang::Wrong number of input arguments');
end

% Some initializations
r     = zeros(n, 1);
c     = zeros(n, 1);
sr    = ones(n + 1, 1);
sc    = ones(n + 1, 1);
temp  = 1;


% Get the dimensions of the block matrices and sum them up by considering
% only the objects on the diagonal.
for i = 1 : n
    [r(i), c(i)] = size(varargin{temp});
    sr(i + 1)    = sr(i) + r(i);
    sc(i + 1)    = sc(i) + c(i);
    temp         = temp + n + 1 - i;
end

% Check whether an ss object is contained in the input.
issys = 0;
for i = 1:anz
   if isa(varargin{i}, 'ss')
      issys = 1;
      break;
   end
end

% Construct the block left lower triangular matrix.
index = 1;
A     = zeros(sr(end) - 1, sc(end) - 1);
if issys
   A = ss(A); % Convert to ss if an ss object is involved
end
for i = 1 : n
    rows = sr(i):sr(i+1)-1;
    for j = 1 : n
       if i <= j
          columns = sc(j):sc(j+1)-1; 
          % Handle sortcuts
          if     isa(varargin{index}, 'char') && varargin{index} == 'z'
              A(rows, columns) = zeros(r(i), c(j));
          elseif isa(varargin{index}, 'char') && varargin{index} == 'e'
              A(rows, columns) = eye(r(i), c(j));
          else
              [tr, tc] = size(varargin{index});
              % Another sanity check
              if length(rows) ~= tr || length(columns) ~= tc
                error('blkrutriang::check dimensions please')
              end
              A(rows, columns) = varargin{index};
          end
          index = index + 1;
       end
    end
end


end

