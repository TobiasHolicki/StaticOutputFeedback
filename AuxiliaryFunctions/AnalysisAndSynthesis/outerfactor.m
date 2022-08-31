%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : outerfactor.m                                                 %
%                                                                         %
% Author  : Tobias Holicki                                                %
% Version : 04                                                            %
% Date    : 06.10.2020                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input to this function is a system sys with realization
%       |A  |  B1, ...  BN | 
%       |------------------|
%       |C1 | D11, ... D1N |
%       |:  |  :         : |
%       |CN | DN1, ... DNN |
% and partition of the input and outputs described by inp and out,
% respectively. 
% If varargin = 'synthesis', 'syn' or 's' this functions returns the outer 
% factors as appearing in LMI based controller design. These are of the 
% form
% [   A,    B1, ...,    BN-1;         [    I,        0, ...,         0;
%     I,     0, ...,       0;            -A',     -C1', ...,    -CN-1';
%    C1,   D11, ...,   D1N-1;              0,        I, ...,         0;
%     0,     I, ...,       0;           -B1',    -D11', ...,   -DN-11';
%     :      :    :        :                :        :    :          :
%  CN-1, DN-11, ..., DN-1N-1;               0,       0, ...,         I;
%     0,     0, ...,       I] * V,     -BN-1', -D1N-1', ..., -DN-1N-1'] * U
% with V and U beinf the annihilators of
%           [CN, DN1, ..., DNN-1]  and  [BN', D1N', ..., DN-1N],
% respectively.
% If varargin = 'analysis', 'ana' or 'a' this functions returns the outer
% factors as appearing in an LMI based analysis. These are of the form
% [ A,  B1, ...,  BN;     [   I,     0, ...,     0;
%   I,   0, ...,   0;       -A',  -C1', ...,  -CN';
%  C1, D11, ..., D1N;         0,     I, ...,     0;
%   0,   I, ...,   0;      -B1', -D11', ..., -DN1';
%   :    :    :    :          :      :    :      :
%  CN, DN1, ..., DNN;         0,     0, ...,     I;
%   0,   0, ...,   I],     -BN', -D1N', ..., -DNN'].
%
% ----- Input ---------------------------------------------------------- 
%      sys - An LTI system as ss object
%      inp - Partition of input channels
%      out - Partition of output channels
% varargin - Specification for synthesis or analysis mode (default: synth)
% ----- Output ---------------------------------------------------------
%   Oprim - Outer factor corresponding to the primal LMI
%   Odual - Outer factor corresponding to the dual LMI
%
function [Oprim, Odual] = outerfactor(sys, inp, out, varargin)

% Some sanity checks
if length(inp) ~= length(out)
    error(['outerfactor::The number of input and output channels ', ...
           'should be the same']);
end
if sum(inp) ~= size(sys, 2)
    error(['outerfactor::The number of inputs and the input ', ...
           'partition do not match']);
end
if sum(out) ~= size(sys, 1)
    error(['outerfactor::The number of outputs and the output ', ...
           'partition do not match']);
end
if length(varargin) > 1
    error('outerfactor::Wrong number of input arguments');
elseif isempty(varargin)
    mode = 'synthesis';
elseif strcmp(varargin{1}, 'synthesis') || ...
       strcmp(varargin{1}, 'syn') || strcmp(varargin{1}, 's')
    mode = 'synthesis';
elseif strcmp(varargin{1}, 'analysis') || ...
       strcmp(varargin{1}, 'ana') || strcmp(varargin{1}, 'a')
    mode = 'analysis';
else
    error('outerfactor::Please correctly specify the desired mode');
end

% Ensure that inp and out are row vectors
inp = inp(:)';
out = out(:)';


if strcmp(mode, 'synthesis')
    % Abbreviations
    n   = length(inp);         % Number of channels 
    mea = out(end);            % Number of measurements
    act = inp(end);            % Number of actuators
    la  = size(sys.a, 1);      % Number of states

    % Get relevant system data
    [A,  B,  C,  D] = ssdata(sys(1 : end - mea, 1 : end - act));
    [~,  ~, Cp, Dp] = ssdata(sys(end - mea + 1 : end, 1 : end - act));
    [~, Bd,  ~, Dd] = ssdata(sys(1 : end - mea, end - act + 1 : end));

    % Prepare permutation
    % IO dimensions plus 2*states, but excluding the control channel
    io  = [la, inp(1:end-1), la, out(1:end-1)];  
    lio = sum(io); % Dimension of all of those

    % Group the channels
    p = mat2cell(1 : lio, 1, io); 

    % Construct relevant permutation
    pn = zeros(1, 2*n);
    for i = 1 : n
        pn(2*i-1) = i+n;
        pn(2*i)   = i;
    end
    p = cell2mat(p(pn));

    % Apply the permutation
    Oprim = [eye(la + sum(inp(1:end-1))); A, B; C, D];
    Oprim = Oprim(p, :) * null([Cp, Dp]);

    Odual = [-A', -C'; -B', -D'; eye(la + sum(out(1:end-1)))];
    Odual = Odual(p, :) * null([Bd', Dd']);
    
elseif strcmp(mode, 'analysis')
    % Abbreviations
    n   = length(inp)+1;       % Number of channels (including states)
    la  = size(sys.a, 1);      % Number of states

    % Get relevant system data
    [A, B, C, D] = ssdata(sys);

    % Prepare permutation
    % Input and output dimensions plus 2*states
    io  = [la, inp, la, out];  
    lio = sum(io); % Dimension of all of those

    % Group the channels
    p = mat2cell(1 : lio, 1, io); 

    % Construct relevant permutation
    pn = zeros(1, 2*n);
    for i = 1 : n
        pn(2*i-1) = i+n;
        pn(2*i)   = i;
    end
    p = cell2mat(p(pn));

    % Apply the permutation
    Oprim = [eye(la + sum(inp)); A, B; C, D];
    Oprim = Oprim(p, :);

    Odual = [-A', -C'; -B', -D'; eye(la + sum(out))];
    Odual = Odual(p, :);
end

end

