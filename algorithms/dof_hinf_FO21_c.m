function out = dof_hinf_FO21_c(sys, varargin)

%% Dimensões
n = size(sys.A,1);
mu = size(sys.Bu,2);
mw = size(sys.Bw,2);
pz = size(sys.Cz,1);
py = size(sys.Cy,1);


if nargin > 1
    if nargin == 2
        options = varargin{1};
    else
        options = struct(varargin{:});
    end
else
    options = [];
end
if ~isfield(options,'xi')
    options.xi = 1;
end

if ~isfield(options,'gamma')
    options.gamma = [];
end
if ~isfield(options,'solver')
    options.solver = 'mosek';
end
if ~isfield(options,'echo')
    options.echo = 0;
end
if ~isfield(options,'maxIt')
    options.maxIt = 15;
end
if ~isfield(options,'tol')
    options.tol = 1e-7;
end
if ~isfield(options, 'nc')
    options.nc = 0;
end
if ~isfield(options, 'tolGamma')
    options.tolGamma = 1e-3;
end
if ~isfield(options,'B0')
    options.B0 = [eye(n+options.nc) options.xi*eye(n+options.nc) -eye(n+options.nc)];
end
if ~isfield(options,'X3X4')
    options.X3X4 = 'full';
end
if ~isfield(options,'X1X2X5')
    options.X1X2X5 = 'full';
end
if ~isfield(options,'earlyCheck')
    options.earlyCheck = 1;
end
if ~isfield(options,'Kini')
    options.Kini = [];
end
if ~isfield(options,'simTrans')
    options.simTrans = 0;
end


%% Sistema a ser controlado

A = sys.A;
Bu = sys.Bu;
Bw = sys.Bw;
Cz = sys.Cz;
Cy = sys.Cy;
Dzu = sys.Dzu;
Dzw = sys.Dzw;
Dyw = sys.Dyw;

[Acl,Bcl,Ccl,Dcl,L] = setupClosedLoop(A,Bu,Bw,Cz,Cy,Dzu,Dzw,Dyw,options);

out.feas = 0;
if isempty(options.Kini)    
	phase1 = phaseStabilization(Acl,L,n,options);
    if phase1.feas
        assign(L,phase1.L);
    else
        return;
    end
else
    assign(L,options.Kini);
    phase1.feas = 1;
end
if phase1.feas    
    gini=options.gamma;
    options.gamma=[];
    iniCond = evaluateFirstNorm(double(Acl),double(Bcl),double(Ccl),double(Dcl),n,mw,pz,options);    
    if iniCond.feas
        if ~isempty(gini)
            options.gamma = iniCond.hinf*gini;
            iniCond = evaluateFirstNorm(double(Acl),double(Bcl),double(Ccl),double(Dcl),n,mw,pz,options);
        end
        options.B0 = [iniCond.X1' iniCond.X2' iniCond.X3' iniCond.X4' -eye(n+options.nc)];
        options.gamma=[];
        out = phaseHinf(Acl,Bcl,Ccl,Dcl,L,n,mw,pz,options);        
    end
end

%__________________________________________________________________________
function [Acl,Bcl,Ccl,Dcl,L] = setupClosedLoop(A,Bu,Bw,Cz,Cy,Dzu,Dzw,Dyw,options)

[n,mu]=size(Bu);
[py,n]=size(Cy);
nc=options.nc;
%% Ganho por Realimentação de Saída

if options.nc == 0    
    L = sdpvar(mu, py,'full');        
    %% Sistema em Malha Fechada
    Acl = A+Bu*L*Cy;
    Bcl = Bu*L*Dyw+Bw;
    Ccl = Cz + Dzu*L*Cy;
    Dcl = Dzu*L*Dyw+Dzw;
else
    Ac =  sdpvar(nc, nc,'full');
    Bc =  sdpvar(nc, py,'full');
    Cc =  sdpvar(mu, nc,'full');
    Dc =  sdpvar(mu, py,'full');
        
    Acl = [A+Bu*Dc*Cy Bu*Cc;Bc*Cy Ac];
    Bcl = [Bw+Bu*Dc*Dyw;Bc*Dyw];
    Ccl = [Cz+Dzu*Dc*Cy Dzu*Cc];
    Dcl = Dzw+Dzu*Dc*Dyw;
    L = [Ac Bc;Cc Dc];
end

if options.simTrans
    T =[eye(n) zeros(n,nc);
         eye(nc) zeros(nc,n-nc) -eye(nc)];
    
    iT=inv(T);
    Acl = iT*Acl*T;
    Bcl = iT*Bcl;
    Ccl = Ccl*T;
end
%__________________________________________________________________________
function out = phaseHinf(Acl,Bcl,Ccl,Dcl,L,n,mw,pz,options)

LMIsGam = [];
rgamma = sdpvar;

P = sdpvar((n+options.nc), (n+options.nc), 'symmetric');

Q11 = zeros(n+options.nc);
Q21 = P;
Q22 = zeros(n+options.nc);
Q31 = zeros(mw,n+options.nc);
Q32 = zeros(mw,n+options.nc);
Q33 = -rgamma*eye(mw);
Q41 = Ccl;
Q42 = zeros(pz,n+options.nc);
Q43 = Dcl;
Q44 = -rgamma*eye(pz);
Ir=blkdiag(eye(n),eye(options.nc));
Q51 = Acl;
Q52 = -eye(n+options.nc);
Q53 = Bcl;
Q54 = zeros(n+options.nc,pz);
Q55 = zeros(n+options.nc);

Q = [Q11 Q21' Q31' Q41' Q51';
    Q21 Q22 Q32' Q42' Q52';
    Q31 Q32 Q33 Q43' Q53';
    Q41 Q42 Q43 Q44 Q54';
    Q51 Q52 Q53 Q54 Q55];

X11 = sdpvar((n+options.nc), (n+options.nc), options.X1X2X5);
X21 = sdpvar((n+options.nc), (n+options.nc), options.X1X2X5);
B11 = sdpvar((n+options.nc), (n+options.nc), options.X1X2X5);
B21 = sdpvar((n+options.nc), (n+options.nc), options.X1X2X5);
if ~isempty(options.X3X4)
    X31 = sdpvar(mw,     (n+options.nc), options.X3X4);
    X41 = sdpvar(pz,     (n+options.nc), options.X3X4);
    B31 = sdpvar(mw,     (n+options.nc), options.X3X4);
    B41 = sdpvar(pz,     (n+options.nc), options.X3X4);
else
    X31 = zeros(mw,     (n+options.nc));
    B31 = zeros(mw,     (n+options.nc));
    X41 = zeros(pz,     (n+options.nc));
    B41 = zeros(pz,     (n+options.nc));
end
X51 = sdpvar((n+options.nc), (n+options.nc), options.X1X2X5);
B51 = sdpvar((n+options.nc), (n+options.nc), options.X1X2X5);

X = [X11;X21;X31;X41;X51];
B = [B11;B21;B31;B41;B51]';

LMIs = [P>=0,Q+X*B+B'*X'<=0];
set = sdpsettings('verbose',0,'solver',options.solver, 'warning', 0);
HINF = optimizer(LMIs,rgamma, set, B, {rgamma,X, P,L});

out.it = 1;
out.feas = 1;
out.gammas = [];
obj = rgamma;
Bcal = options.B0;
out.cpusec_s =0;
while out.it <= options.maxIt
    t=clock;
    sol = HINF(Bcal);
    t=etime(clock,t);
    out.cpusec_s =out.cpusec_s + t;    
    out.L = sol{4};
    Bcal = sol{2}';      
    assign(L,out.L);
    out.gammas = [out.gammas sol{1}];
    if options.echo
        maxNorm=hinfnorm(ss(double(Acl),double(Bcl),double(Ccl),double(Dcl)),1e-3);
        fprintf('fase 2 it = %d -> gamma = %.4f (%.4f)\n',out.it,sol{1},maxNorm);
    end
    if size(out.gammas,2)> 2
        if abs(out.gammas(1,end-1) -  out.gammas(1,end)) < options.tolGamma
            break
        end
    end    
    out.it=out.it+1;    
end

if out.feas
    out.hinf = out.gammas(1,end);
    out.X = sol{2};
   
    out.maxNorm=hinfnorm(ss(double(Acl),double(Bcl),double(Ccl),double(Dcl)),1e-4);    
end

%__________________________________________________________________________
function out = evaluateFirstNorm(Acl,Bcl,Ccl,Dcl,n,mw,pz,options)

P=sdpvar(n+options.nc,n+options.nc,'symmetric');
LMIsGam = [];
rgamma=sdpvar;
if ~isempty(options.gamma)
    obj_gamma = [];
    LMIsGam = [rgamma<= options.gamma];
else
    obj_gamma=rgamma;
end


Q11 = zeros(n+options.nc);
Q21 = P;
Q22 = zeros(n+options.nc);
Q31 = zeros(mw,n+options.nc);
Q32 = zeros(mw,n+options.nc);
Q33 = -rgamma*eye(mw);
Q41 = Ccl;
Q42 = zeros(pz,n+options.nc);
Q43 = Dcl;
Q44 = -rgamma*eye(pz);

Q = [Q11 Q21' Q31' Q41';
    Q21 Q22 Q32' Q42';
    Q31 Q32 Q33 Q43';
    Q41 Q42 Q43 Q44];
Bcal = [Acl -eye(n+options.nc) Bcl zeros(n+options.nc,pz)];

X1 = sdpvar((n+options.nc), (n+options.nc), options.X1X2X5);
X2 = sdpvar((n+options.nc), (n+options.nc), options.X1X2X5);
if ~isempty(options.X3X4)
    X3 = sdpvar(mw,     (n+options.nc), options.X3X4);
    X4 = sdpvar(pz,     (n+options.nc), options.X3X4);
else
    X3 = zeros(mw,     (n+options.nc));
    X4 = zeros(pz,     (n+options.nc));
end

X =[X1;X2;X3;X4];
LMIs = [LMIsGam,P>=0,Q+X*Bcal+Bcal'*X'<=0];
sol = solvesdp(LMIs,obj_gamma,sdpsettings('verbose',0,'solver','mosek'));

out.p=min(checkset(LMIs));
out.X1=double(X1);
out.X2=double(X2);
if ~isempty(options.X3X4)
    out.X3=double(X3);
    out.X4=double(X4);
else
    out.X3=X3;
    out.X4=X4;
end
out.hinf=double(rgamma);
out.feas = (out.p>-1e-6);

%__________________________________________________________________________
function out = phaseStabilization(Acl,L,n,options)

rho = sdpvar;

P = sdpvar((n+options.nc), (n+options.nc), 'symmetric');

Q11 = zeros(n+options.nc);
Q21 = P;
Q22 = zeros(n+options.nc);
Q31 = Acl - rho*eye(n+options.nc);
Q32 = -eye(n+options.nc);
Q33 = zeros(n+options.nc);

Q = [Q11 Q21' Q31';
     Q21 Q22 Q32';
     Q31 Q32 Q33];

X1 = sdpvar((n+options.nc), (n+options.nc),  options.X1X2X5);
X2 = sdpvar((n+options.nc), (n+options.nc),  options.X1X2X5);
X3 = sdpvar((n+options.nc), (n+options.nc),  options.X1X2X5);

B1 = sdpvar((n+options.nc), (n+options.nc),  options.X1X2X5);
B2 = sdpvar((n+options.nc), (n+options.nc),  options.X1X2X5);
B3 = sdpvar((n+options.nc), (n+options.nc),  options.X1X2X5);

X = [X1;X2;X3];
B = [B1;B2;B3]';

LMIs = [P >= 0, Q + X*B + B'*X' <= 0];
set = sdpsettings('verbose',0,'solver',options.solver, 'warning', 0);
STAB = optimizer(LMIs,rho, set, B, {rho,X, P,L});
out.it = 1;
out.feas = 0;
out.r = [];
out.cpusec_s =0;
Bcal = options.B0;
while out.it <= options.maxIt
    t=clock;
    sol = STAB(Bcal);
    t=etime(clock,t);    
    out.cpusec_s =out.cpusec_s + t;    
    Bcal = sol{2}';
    out.L = sol{4};
    assign(L,out.L);
    out.it=out.it+1;
    out.r = [out.r sol{1}];
    muOpt = max(real(eig(double(Acl))));
    if options.echo == 1
        fprintf('fase 1 - it %d -> r %.4f (muOpt = %.4f) \n',out.it, sol{1},muOpt);
    end
    if  sol{4} <= 0                
        out.feas = 1;
        break;        
    end
    if options.earlyCheck & muOpt < 0         
        out.feas = 1;
        break;
    end    
    
    %if out.perturb == 0 & out.it > 10 & out.r(1,end-1)-out.r(1,end) < 0.01
    if  out.it > 2 & abs((out.r(end-1)-out.r(end))/out.r(end-1)) < 0.05
        break;
    end
end


