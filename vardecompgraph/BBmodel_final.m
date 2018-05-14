clear all; clc; close all;

% IN THIS VERSION, I LINEARIZE VARIABLES THAT CANNOT BE LOGLINEARIZED. 
% This applies to Debt, Trade balances, interest rates (R)
% WE ALSO LINEARIZE TRADE BALANCE OVER OUTPUT
% IN THE OTHER VERSION, I CONVERTED TO DEVIATION FROM STEADY STATE
% COMPARED TO GDP

run('BBparameter.m')
load BBparameter


% Index of State variables
Lambda  = 1;
C       = 2;
N1      = 3;
N2      = 4;
K1      = 5;
K2      = 6;
Mtilde  = 7;
Y       = 8;
Ytilde  = 9;
K       = 10;
R1      = 11;
R2      = 12;
W1      = 13;
W2      = 14;
Ptilde  = 15;
D       = 16;
I       = 17; % INVestment
R       = 18;
TB      = 19;
TBtilde = 20;
Ygdp    = 21;
TBtotal = 22;

A       = 23;
Atilde  = 24;
G       = 25;
S       = 26;
Nu      = 27;
Mu      = 28;
Plag    = 29;

ExpLambda = 30;
ExpNu     = 31;
ExpR1     = 32;
ExpR2     = 33;
ExpK      = 34;
ExpG      = 35;

Klag = 36;

% OBSERVABLES
Ygrowth = 37;
Cgrowth = 38;
Igrowth = 39;
TBoverY = 40;


% Index of expectational errors
% Index of exogenous shocks
EpsA      = 1;
EpsAtilde = 2;
EpsG      = 3;
EpsS      = 4;
EpsNu     = 5;
EpsMu     = 6;
EpsP      = 7;

NY   = 40;
NETA = 6; % # of expectational errors
NEPS = 7; % # of exogenous shocks

g1  = zeros(NY,NY);
g0  = zeros(NY,NY);
Pi  = zeros(NY,NETA);
Psi = zeros(NY,NEPS);
const   = zeros(NY,1);

% Eq1 Marginal Utility Def 
g0(1,Lambda) = -(1/gamma)*(lambda^(-(1/gamma)));
g0(1,C)      = -c;
g0(1,N1)     = (theta/omega)*(n1^omega);
g0(1,N2)     = (theta/omegatilde)*(n2^omegatilde);

% Eq2 Labor Supply1
g0(2,N1) = omega-1;
g0(2,W1) = -1;

% Eq3 Labor Supply2 
g0(3,N2) = omegatilde-1;
g0(3,W2) = -1;

% Eq4 Euler Equation Capital1 
g0(4,ExpNu)     = 1;
g0(4,ExpLambda) = 1;
g0(4,ExpR1)     = b*(g^(-gamma))*r1;
g0(4,ExpK)      = b*(g^(2-gamma))*phi;
g0(4,ExpG)      = b*(g^(2-gamma))*phi;
g0(4,K)         = -phi*g*(1+b*g^(1-gamma));
g0(4,Nu)        = -1;
g0(4,Lambda)    = -1;
g0(4,G)         = -(phi*g+gamma);

g1(4,K)         = -phi*g;


% Eq5 Euler Euqation Capital2
g0(5,ExpR1)    = 1;
g0(5,ExpR2)    = -1;

% Eq6 Euler Equation Debt 

g0(6,ExpLambda) = 1;
g0(6,ExpNu)     = 1;
g0(6,R)      = 1/(1+r); % LINEARIZED
g0(6,Nu)     = -1;
g0(6,Lambda) = -1;
g0(6,G)      = -gamma;

% Eq7 Budget Constraint COME BACK

g0(7,Y)         = y;
g0(7,Ytilde)    = ptilde*ytilde;
g0(7,Ptilde)    = ptilde*(ytilde-mtilde);  
g0(7,Mtilde)    = -ptilde*mtilde;
g0(7,D)         = g/(1+r);  % LINEARIZED
g0(7,G)         = (d*g)/(1+r);
g0(7,R)         = -(d*g)/(1+r)^2; % LINEARIZED
g0(7,C)         = -c;
g0(7,I)         = -INV;
g0(7,S)         = -s;
g1(7,D)         = 1; % LINEARIZED

% Eq8 Interest rate process
g0(8,R)      = 1;
g0(8,D)      = -psi; % LINEARIZED
g0(8,Ptilde) = -xi;
g0(8,Mu)     = -1;

% Eq9 R1 MPK1
g0(9,R1)    = 1;
g0(9,Y)     = -1;
g1(9,K1)    = -1;

% Eq10 R2 MPK2
g0(10,R2)       = 1;
g0(10,Ytilde)   = -1;
g0(10,Ptilde)   = -1;
g1(10,K2)      = -1;

% Eq11 W1 MPN1
g0(11,W1)   = 1;
g0(11,Y)    = -1;
g0(11,N1)   = 1;

% Eq12 W2 MPN2
g0(12,W2)       = 1;
g0(12,Ytilde)   = -1;
g0(12,Ptilde)   = -1;
g0(12,N2)       = 1;

% Eq13 MP of commoditiy
g0(13,Ptilde)   = 1;
g0(13,Y)        = -1;
g0(13,Mtilde)   = 1;

% Eq14 Production fcn of numeraire goods
g0(14,Y)        = 1;
g0(14,A)        = -1;
g0(14,Mtilde)   = -alpm;
g0(14,N1)       = -(1-alpk-alpm);
g0(14,G)        = -(1-alpk-alpm);
g1(14,K1)       = alpk;

% Eq15 Prod fcn of commodity goods
g0(15,Ytilde)   = 1;
g0(15,Atilde)   = -1;
g0(15,N2)       = -(1-alpktilde);
g0(15,G)        = -(1-alpktilde);
g1(15,K2)       = alpktilde;

% Eq16 TB Definition
% We defined TB as TB to Ygdp before. 
% Now just report TB as the level. 
% linearize: TB(t)-TB(ss)

g0(16,TB)   = 1; % LINEARIZED
g0(16,Y)    = -y;
g0(16,C)    = c;
g0(16,I)    = INV;
g0(16,S)    = s;

% Eq17 TBtilde Def

g0(17,TBtilde)  = 1; % LINEARIZED
g0(17,Ptilde)   = -ptilde*(ytilde-mtilde);
g0(17,Ytilde)   = -ptilde*ytilde;
g0(17,Mtilde)   = ptilde*mtilde;

% Eq18 Ygdp Def

g0(18,Ygdp)     = ygdp;
g0(18,Y)        = -y;
g0(18,TBtilde)  = -1; % LINEARIZED

% Eq19 TBtotal Def

g0(19,TBtotal)  = 1; % LINEARIZED
g0(19,TB)       = -1;% LINEARIZED
g0(19,TBtilde)  = -1;% LINEARIZED

% Eq20 K definition K=K1+K2

g0(20,K)    = k;
g0(20,K1)   = -k1;
g0(20,K2)   = -k2;

% Eq21 K lag Definition: Klag(t+1)=K(t)

g0(21,Klag) = 1;
g1(21,K)    = 1;

% Eq22 Ptilde lag definition Plag(t)=P(t-1)

g0(22,Plag)     = 1;
g1(22,Ptilde)   = 1;

% Stochastic process
RHO=diag([rhoa rhoatilde rhog rhos rhonu rhomu]);
SIG=diag([siga sigatilde sigg sigs signu sigmu]);

% Eq23 A
% Eq24 Atilde
% Eq25 G
% Eq26 S
% Eq27 Nu
% Eq28 Mu

g0(23:28,A:Mu)          = eye(6);
g1(23:28,A:Mu)          = RHO;
Psi(23:28,EpsA:EpsMu)   = SIG;

% Eq29 Ptilde

g0(29,Ptilde)   = 1;
g1(29,Ptilde)   = rhop1;
g1(29,Plag)     = rhop2;
Psi(29,EpsP)    = sigp;

% Eq30 INVestment Definition

g0(30,I) = INV;
g0(30,K) = -k*g;
g0(30,G) = -k*g;
g1(30,K) = -(1-delta)*k;

% Expectational errors

g0(31,Nu)=1;
g1(31,ExpNu) = 1;
Pi(31,1) =1;

g0(32,Lambda)=1;
g1(32,ExpLambda) = 1;
Pi(32,2) =1;

g0(33,R1)=1;
g1(33,ExpR1) = 1;
Pi(33,3) =1;

g0(34,K)=1;
g1(34,ExpK) = 1;
Pi(34,4) =1;

g0(35,G)=1;
g1(35,ExpG) = 1;
Pi(35,5) =1;

g0(36,R2)=1;
g1(36,ExpR2) = 1;
Pi(36,6) =1;

% Definition of OBSERVABLES
% I should add constant but I skipped

% Y growth
g0(37,Ygrowth) = 1;
g0(37,Ygdp)    = -1;
g1(37, Ygdp)   = -1;
g1(37, G)      = 1;
const(37) = log(g);

% Cgrowth 
g0(38, Cgrowth) = 1;
g0(38, C)       = -1;
g1(38, C)       = -1;
g1(38, G)       = 1;
const(38) = log(g);

% Igrowth
g0(39, Igrowth) = 1;
g0(39, I)       = -1;
g1(39, I)       = -1;
g1(39, G)       = 1;
const(39) = log(g);

% TB/Ygdp
% TB/YGDP(t)-tb/ygdp = (TB(t)-tb)/ygdp -tb/ygdp *Ygdphat(t) 
g0(40, TBoverY) = 1;
g0(40, TBtotal) = -1/ygdp;
g0(40, Ygdp)    = tbtotal/ygdp;


[G1, C0, G0, fmat, fwt, ywt, gev, eu] = gensys(g0, g1, const, Psi, Pi);
if eu ~= ones(2,1)
    Loss = inf;
    return
end

%% Kalman Filter & Smoother

MeanS = zeros(NY,1);   % Unconditional mean
VarS = zeros(NY,NY); % Sigma^s (1|0)=unconditional variance
error = 1; count = 0; maxit = 1e5; tol = 1e-5;
while (error>tol && count<maxit)
    VarS1 = G1*VarS*(G1.')+ G0*(G0.');
    error = max(max(abs(VarS1-VarS))); count = count+1;
    VarS = VarS1;
end

% Data
DATA = xlsread('DataVAR.xlsx');
DATA = DATA(:,2:5);
diff = DATA(2:end,1:3)-DATA(1:end-1,1:3);
DATA = [diff DATA(2:end,4)];
clear diff
DATA = DATA'; % DATA(:,t);
T = size(DATA,2);

%{ 
Kalman Filter equations:
y(t) = M1*s(t)+u(t), u(t)~N(0,R)
s(t) = G1*s(t-1)+G0*e(t), e(t)~N(0,Q)

Here, we let u(t) = 0. No measurement error. 
%}
% Measurement equations
M1=zeros(4,NY);
M1(1,Ygrowth) = 1;
M1(2,Cgrowth) = 1;
M1(3,Igrowth) = 1;
M1(4,TBoverY) = 1;

% Store the signal extraction part
KFmean = zeros(NY,T);
KFvar  = zeros(NY,NY,T);

% Store the forecasted shock distribution
FCmean = zeros(NY,T);
% Apply Kalman Filter to Data(:,t)
for t = 1:T
    % Initiate the Kalman Filter with s(t|Info till t-1)
    % Forecasting y using the measurement eqn. 
    % Y1 = E(Y(t)|Info till t-1)
    % Y2 = V(Y(t)|Info till t-1)
    Y1 = M1*MeanS;       % Conditional Mean, (NYx1)
    Y2 = M1*VarS*(M1.'); % Conditional Variance, (NYxNY)
    if det(Y2)==0
        disp('error')
        return
    end
    invY2 = (eye(length(Y2))/Y2);
    % Signal extraction using s(t|Into till t-1) and y(t|Info till t-1)
    S1 = MeanS + VarS*(M1.')*invY2*(DATA(:,t)-Y1); % mean, (NSx1)
    S2 = VarS - VarS*(M1.')*invY2*M1*VarS; % variance, (NSxNS)
    % Store results
    KFmean(:,t)  = S1;
    KFvar(:,:,t) = S2;
    % Forecasting s using the transition eqn.
    MeanS = G1*S1; 
    VarS  = G1*S2*(G1.')+G0*(G0.');
end

% Kalman smoother
% Store results
KSmean = zeros(NY,T);
KSvar  = zeros(NY,NY,T);
% Initiate from the last period
KSmean(:,T)  = KFmean(:,T);
KSvar(:,:,T) = KFvar(:,:,T);

for t=T-1:-1:1
    mean_t  = KFmean(:,t);
    var_t   = KFvar(:,:,t);
    % Forecast one period ahead
    mean_t1 = G1*mean_t;
    var_t1  = G1*var_t*(G1')+G0*(G0');
    % Kalman smoothing
    invt1=inversePD(A);
    %invt1 = eye(length(var_t1))/var_t1;
    J = (var_t)*(G1')*invt1;
    KSmean(:,t)  = (mean_t)+J*(KSmean(:,t+1)-mean_t1);
    KSvar(:,:,t) = (var_t)+J*(KSvar(:,:,t+1)-var_t1)*(J');
end

%% Take out Innovations
% s(t)=G1*s(t-1)+G0*e(t)
% s(0) = G0*e(0) => e(0)= inv(G0'*G0)*(G0'*s(0))
% e(t) = inv(G0'*G0)*(G0'*s(t)-G0'*G1*s(t-1))

%Store shocks
Exog = zeros(NEPS,T);
invgg = eye(length(G0'*G0))/(G0'*G0);
Exog(:,1) = invgg*(G0'*KSmean(:,1));
for t=2:T
    Exog(:,t)=invgg*(G0'*KSmean(:,t)-G0'*G1*KSmean(:,t-1));
end

%% Historical decomposition of shocks

d = zeros(1,NY);
d(Ygrowth) = 1;
HistD = zeros(NEPS,T);
HistD(:,1) = (d*G0).*(Exog(:,1)');

for t=2:T
    histdc = 0;
    for j=t-1:-1:0
       histdc = histdc+(d*G1^j*G0).*(Exog(:,t-j)');
    end
    HistD(:,t)=histdc;
end

HD = zeros(4,T);
HD(1,:) = HistD(EpsA,:)+HistD(EpsAtilde,:);
HD(2,:) = HistD(EpsG,:);
HD(3,:) = HistD(EpsP,:);
HD(4,:) = HistD(EpsMu,:)+HistD(EpsS,:)+HistD(EpsNu,:);

h = sum(HD,1);
plot(h); hold on
plot(DATA(1,:))

dgdp = DATA(1,:);
save vardecomp.mat HD dgdp 

function A_=inversePD(A)
%A:positive definite matrix
M=size(A,1);
[R b] = chol(A);
if b~=0
    return
end
R_ = R \ eye(M);
A_ = R_ * R_';
end

 
