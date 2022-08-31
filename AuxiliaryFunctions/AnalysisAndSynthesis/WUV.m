%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : WUV.m                                                         %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 03                                                            %
% Date    : 06.10.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a given LTI system sys
%  dot x =  A x +  B1 p +  B2 u
%      q = C1 x + D11 p + D12 u
%      y = C2 x + D21 p
% this function returns matrices W, U, V such that the closed-loop system 
% matrices of the above system with a controller of the form
%  dot xc =  Ac xc +  Bc1 wc +  Bc2 y
%      zc = Cc1 xc + Dc11 wc + Dc12 y
%       u = Cc2 xc + Dc21 wc + Dc22 y
% are exactly given by W + U^T K V. Here, 
%     | Ac  Bc1  Bc2|
% K = |Cc1 Dc11 Dc12|
%     |Cc2 Dc21 Dc22| 
% and the signals of the closed-loop system are ordered as follows:
%  | dot x|                 | x|
%  |dot xc|                 |xc|
%  |     q| = (W + U^T K V) | p|.
%  |    zc|                 |wc|
% Such a decomposition is required for applying the elimination lemma.
%
% Note that for actually applying the elimination lemma afterwards, it
% might be required to perform a permutation of the matrices W, U, V as 
% typically q = [z; e] and p = [w; d] (e.g. for gain-scheduling synthesis).
% One has to make sure that the singals are correctly ordered!
%
% ----- Input ---------------------------------------------------------- 
%   sys - Given system with description as above
%   nu  - Number of actuators
%   ny  - Number of measurements
%   nxc - Degree of a to-be-designed controller 
%   nwc - Input dimension of the schedulung channel of the controller
%   nzc - Output dimension of the scheduling channel of the controller
% ----- Output ---------------------------------------------------------
%   W, U, V - See above
%
function [W, U, V] = WUV(sys, nu, ny, varargin)

if isempty(varargin)
    nxc = 0;
    nwc = 0;
    nzc = 0;
elseif length(varargin) == 1
    nxc = varargin{1};
    nwc = 0;
    nzc = 0;
elseif length(varargin) == 2
    nxc = varargin{1};
    nwc = varargin{2};
    nzc = nwc;
elseif length(varargin) == 3
    nxc = varargin{1};
    nwc = varargin{2};
    nzc = varargin{3};
else
    error('WUV::Wrong number of input arguments');
end

% System data
[A, B1, C1, D11] = ssdata(sys(1:end-ny, 1:end-nu));
[~, B2,  ~, D12] = ssdata(sys(1:end-ny, end-nu+1:end));
[~,  ~, C2, D21] = ssdata(sys(end-ny+1:end, 1:end-nu));

% Some lengths
[nq, np] = size(D11);
nx       = size(A, 1);

W = [blkdiag(A, zeros(nxc)), blkdiag(B1, zeros(nxc, nwc)); ...
     blkdiag(C1, zeros(nzc, nxc)), blkdiag(D11, zeros(nzc, nwc))];
U = [zeros(nx, nxc + nzc), B2; ...
     eye(nxc), zeros(nxc, nzc + nu); ...
     zeros(nq, nxc + nzc), D12; ...
     zeros(nzc, nxc), eye(nzc), zeros(nzc, nu)]';
V = [zeros(nxc, nx), eye(nxc), zeros(nxc, np + nwc); ...
     zeros(nwc, nx + nxc + np), eye(nwc); ...
     C2, zeros(ny, nxc), D21, zeros(ny, nwc)];
end

