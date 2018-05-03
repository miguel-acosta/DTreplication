clear all; clc; close all;
% This file stores the parameters used in Drechsel and Tenreyro (2018).


% Parameter
ptilde     = 0.5244;
dstar      = -0.001;
s          = 0.0189;
xi         = -0.199;
g          = 1.0117204;
psi        = 2.8;
alpk       = 0.32;
alpm       = 0.05;
alpktilde  = 0.32;
delta      = 0.1255;
phi        = 6;
%b          = 0.93;
b          = 0.9224;

gamma      = 2;
theta      = 1.6;
omega      = 1.6;
omegatilde = 1.6;

% Steady State
r = (1/b)*(g^(gamma))-1;
d = dstar;
rstar = r;
r1 = (1/b)*(g^gamma)-(1-delta);
r2 = r1;

ytilde_k2 = r2/(alpktilde*ptilde);
k2_n2     = g*(ytilde_k2)^(-1/(1-alpktilde));
ytilde_n2 = g^(1-alpktilde)*(k2_n2)^alpktilde;
w2        = (1-alpktilde)*(ytilde_n2)*ptilde;
n2        = (w2/theta)^(1/(omegatilde-1));
k2        = (k2_n2)*n2;
ytilde    = (ytilde_n2)*n2;

y_k1     = r1/alpk;
y_mtilde = ptilde/alpm;
syms k1_n1 mtilde_n1
eqns  = [g^(1-alpk-alpm)*(k1_n1)^(alpk-1)*(mtilde_n1)^(alpm) == y_k1,...
         g^(1-alpk-alpm)*(k1_n1)^(alpk)*(mtilde_n1)^(alpm-1) == y_mtilde]; 
S     = vpasolve(eqns, [k1_n1 mtilde_n1]);
k1_n1     = real(S.k1_n1);
mtilde_n1 = real(S.mtilde_n1);

y_n1    = g^(1-alpk-alpm)*(k1_n1)^(alpk)*(mtilde_n1)^(alpm);
w1      = (1-alpk-alpm)*(y_n1);
n1      = (w1/theta)^(1/(omega-1));
k1      = (k1_n1)*n1;
mtilde  = (mtilde_n1)*n1;
y       = (y_n1)*n1;

k       = k1+k2;
inv     = k*g-(1-delta)*k;
c       = y+ptilde*ytilde-ptilde*mtilde-inv+d*(g/(1+r)-1)-s;
lambda  = (c-(theta/omega)*n1^omega-(theta/omegatilde)*n2^omegatilde)^(-gamma);
tb      = y-c-k*g+(1-delta)*k-s;
tbtilde = ptilde*(ytilde-mtilde);
ygdp    = y+tbtilde;
tbtotal = tb+tbtilde;

% Stochastic process

rhoa      = 0.9;
rhoatilde = 0.9;
rhog      = 0.9;
rhos      = 0.9;
rhonu     = 0.9;
rhomu     = 0.9;

siga      = 0.1;
sigatilde = 0.1;
sigg      = 0.1;
sigs      = 0.1;
signu     = 0.1;
sigmu     = 0.1;

%{
rhoa = 0.8277;
rhoatilde = 0.5887;
rhog = 0.5244;
rhos = 0.6440;
rhonu = 0.8687;
rhomu = 0.9199;

siga = 0.0295;
sigatilde = 0.0525;
sigg = 0.0261;
sigs = 0.1876;
signu = 0.4582;
sigmu = 0.0547;
%}

rhop1 = 0.95;
rhop2 = -0.13;
sigp  = 0.1064;

save BBparameter