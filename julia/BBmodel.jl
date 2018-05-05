using Gensys
function BBmodel(PAR::Dict{String,Float64})    
    @unpack p_til, dstar, s, ξ, g, Ψ, αk, αm, αk_til,
    δ, ϕ,b, Γ, θ, ω, ω_til,
    ρ_a, ρ_a_til, ρ_g, ρ_s, ρ_ν, ρ_μ,
    σ_a, σ_a_til, σ_g, σ_s, σ_ν, σ_μ, ρ_p1, ρ_p2, σ_p,
    d, r, rstar, r1, r2, y_til_k2, k2_n2, y_til_n2, w2, n2,
    k2, y_til, y_k1, y_m_til, k1_n1, m_til_n1, y_n1, w1,
    n1, k1, m_til, y, k, INV, c, λ,
    tb, tb_til, ygdp, ygdpf, tbtotal = PAR
    


# Structure to hold the index of each variable


# Index of State variables
indEndo = Dict("λ"        => 1,
               "C"        => 2,
               "N1"       => 3,
               "N2"       => 4,
               "K1"       => 5,
               "K2"       => 6,
               "M_til"    => 7,
               "Y"        => 8,
               "Y_til"    => 9,
               "K"        => 10,
               "R1"       => 11,
               "R2"       => 12,
               "W1"       => 13,
               "W2"       => 14,
               "P_til"    => 15,
               "D"        => 16,
               "I"        => 17,
               "R"        => 18,
               "TB"       => 19,
               "TB_til"   => 20,
               "Ygdp"     => 21,
               "TBtotal"  => 22,
               "A"        => 23,
               "A_til"    => 24,
               "G"        => 25,
               "S"        => 26,
               "ν"        => 27,
               "μ"        => 28,
               "Plag"     => 29,
               "TBYobs"   => 36, 
               "ΔGDPobs"  => 37,
               "ΔCobs"    => 38,
               "ΔINVobs"  => 39)

# Index of expectational errors
indExpE = Dict("Exp_λ"    => 30,
               "Exp_ν"    => 31,
               "Exp_R1"   => 32,
               "Exp_G"    => 33,
               "Exp_R2"   => 34,
               "Exp_K"    => 35)



# Index of exogenous shocks
indExog = Dict("ϵ_A"     => 1,
               "ϵ_A_til" => 2,
               "ϵ_G"     => 3,
               "ϵ_S"     => 4,
               "ϵ_ν"     => 5,
               "ϵ_μ"     => 6,
               "ϵ_P"     => 7)

ind  = merge(indEndo, merge(indExpE, indExog))
NETA = length(keys(indExpE)) # number of expectational errors
NY   = length(keys(indEndo)) +NETA # Number of endogenous variables
NEPS = length(keys(indExog)) # number of exogenous shocks

g1       = zeros(NY,NY);
g0       = zeros(NY,NY);
Pi       = zeros(NY,NETA);
Psi      = zeros(NY,NEPS);
constant = zeros(NY,1);

# Eq1 Marginal Utility Def  #checked.
g0[1,ind["λ"]]      = -1/Γ*(λ^(-1/Γ));
g0[1,ind["C"]]      = -c;
g0[1,ind["N1"]]     = θ*(n1^ω); 
g0[1,ind["N2"]]     = θ*(n2^ω_til);

# Eq2 Labor Supply1 #checked.
g0[2,ind["N1"]]     = ω-1;
g0[2,ind["W1"]]     = -1;

# Eq3 Labor Supply2  #checked.
g0[3,ind["N2"]]     = ω_til-1;
g0[3,ind["W2"]]     = -1;


# Eq4 Euler Equation Capital1
g0[4,ind["Exp_ν"]]     = 1;
g0[4,ind["Exp_λ"]]     = 1;
g0[4,ind["Exp_R1"]]    = b*(g^(-Γ))*r1;
g0[4,ind["Exp_K"]]     = b*(g^(2-Γ))*ϕ;
g0[4,ind["Exp_G"]]     = b*(g^(2-Γ))*ϕ;
g0[4,ind["ν"]]         = -1;
g0[4,ind["λ"]]         = -1;
g0[4,ind["K"]]         = -ϕ*g*(1+b*g^(1-Γ));
g0[4,ind["G"]]         = -ϕ*g+Γ;

g1[4,ind["K"]]         = -ϕ*g;


# Eq5 Euler Euqation Capital2
g0[5,ind["Exp_R1"]]    = 1;
g0[5,ind["Exp_R2"]]    = -1;

# Eq6 Euler Equation Debt

g0[6,ind["Exp_λ"]] = 1;
g0[6,ind["Exp_ν"]]     = 1;
g0[6,ind["R"]]         = 1/(1+r);
g0[6,ind["ν"]]         = -1;
g0[6,ind["λ"]]         = -1;
g0[6,ind["G"]]         = -Γ;

# Eq7 Budget Constraint

g0[7,ind["Y"]]         = y;
g0[7,ind["Y_til"]]     = p_til*y_til;
g0[7,ind["P_til"]]     = p_til*(y_til-m_til);
g0[7,ind["M_til"]]     = -p_til*m_til;
g0[7,ind["D"]]         = (d*g)/(1+r);
g0[7,ind["G"]]         = (d*g)/(1+r);
g0[7,ind["R"]]         = -(d*g)/(1+r)^2;
g0[7,ind["C"]]         = -c;
g0[7,ind["I"]]         = -INV;
g0[7,ind["S"]]         = -s;

g1[7,ind["D"]]         = d;

# Eq8 Interest rate process
g0[8,ind["R"]]         = 1;
g0[8,ind["D"]]         = -Ψ*ygdp; ## We fixed this from d to ygdp because debt can't be log linearized
g0[8,ind["P_til"]]     = -ξ;
g0[8,ind["μ"]]         = -1;

# Eq9 R1 MPK1
g0[9,ind["R1"]]        = 1;
g0[9,ind["Y"]]         = -1;
g1[9,ind["K1"]]        = -1;

# Eq10 R2 MPK2
g0[10,ind["R2"]]       = 1;
g0[10,ind["Y_til"]]    = -1;
g0[10,ind["P_til"]]    = -1;
g1[10,ind["K2"]]       = -1;

# Eq11 W1 MPN1
g0[11,ind["W1"]]       = 1;
g0[11,ind["Y"]]        = -1;
g0[11,ind["N1"]]       = 1;

# Eq12 W2 MPN2
g0[12,ind["W2"]]       = 1;
g0[12,ind["Y_til"]]    = -1;
g0[12,ind["P_til"]]    = -1;
g0[12,ind["N2"]]       = 1;

# Eq13 MP of commoditiy
g0[13,ind["P_til"]]   = 1;
g0[13,ind["Y"]]       = -1;
g0[13,ind["M_til"]]   = 1;

# Eq14 Production fcn of numeraire goods
g0[14,ind["Y"]]        = 1;
g0[14,ind["A"]]        = -1;
g0[14,ind["M_til"]]    = -αm;
g0[14,ind["N1"]]       = -(1-αk-αm);
g0[14,ind["G"]]        = -(1-αk-αm);

g1[14,ind["K1"]]       = αk;

# Eq15 Prod fcn of commodity goods
g0[15,ind["Y_til"]]    = 1;
g0[15,ind["A_til"]]    = -1;
g0[15,ind["N2"]]       = -(1-αk_til);
g0[15,ind["G"]]        = -(1-αk_til);

g1[15,ind["K2"]]       = αk_til;

# Eq16 TB Definition
g0[16,ind["TB"]]       = ygdp;
g0[16,ind["Y"]]        = -y;
g0[16,ind["C"]]        = c;
g0[16,ind["I"]]        = INV;
g0[16,ind["S"]]        = s;
g0[16,ind["Ygdp"]]     = tb;

# Eq17 TB_til Def
g0[17,ind["TB_til"]]   = ygdp;
g0[17,ind["P_til"]]    = -p_til*(y_til-m_til);
g0[17,ind["Y_til"]]    = -p_til*y_til;
g0[17,ind["M_til"]]    = p_til*m_til;
g0[17,ind["Ygdp"]]     = tb_til;

# Eq18 Ygdp Def
g0[18,ind["Ygdp"]]     = ygdp;
g0[18,ind["Y"]]        = -y;
g0[18,ind["P_til"]]    = p_til*(m_til-y_til);
g0[18,ind["Y_til"]]    = -p_til*y_til;
g0[18,ind["M_til"]]    = p_til*m_til;

# Eq19 TBtotal Def
g0[19,ind["TBtotal"]]  = 1;
g0[19,ind["TB"]]       = -1;
g0[19,ind["TB_til"]]   = -1;

# Eq20 K definition K=K1+K2
g0[20,ind["K"]]        = k;
g0[20,ind["K1"]]       = -k1;
g0[20,ind["K2"]]       = -k2;


# Eq21 investment Definition
g0[21,ind["I"]]        = INV;
g0[21,ind["K"]]        = -k*g;
g0[21,ind["G"]]        = -k*g;
g1[21,ind["K"]]        = -(1-δ)*k;

# Eq22 P_til lag definition Plag(t)=P(t-1)
g0[22,ind["Plag"]]     = 1;
g1[22,ind["P_til"]]    = 1;

# Stochastic process 
RHO=diagm([ρ_a, ρ_a_til, ρ_g, ρ_s, ρ_ν, ρ_μ]);
SIG=diagm([σ_a, σ_a_til, σ_g, σ_s, σ_ν, σ_μ]);

# Eq23 A
# Eq24 A_til
# Eq25 G
# Eq26 S
# Eq27 ν
# Eq28 μ
g0[23:28,ind["A"]:ind["μ"]]          = eye(6);
g1[23:28,ind["A"]:ind["μ"]]          = RHO;
Psi[23:28,ind["ϵ_A"]:ind["ϵ_μ"]]     = SIG;

# Eq29 P_til
g0[29,ind["P_til"]]     = 1;
g1[29,ind["P_til"]]     = ρ_p1;
g1[29,ind["Plag"]]      = -ρ_p2;
Psi[29,ind["ϵ_P"]]      = σ_p;

# Definition with expectational errors
g0[30,ind["λ"]]         = 1;
g1[30,ind["Exp_λ"]]     = 1;
Pi[30,1]                = 1;

g0[31,ind["ν"]]         = 1;
g1[31,ind["Exp_ν"]]     = 1;
Pi[31,2]                = 1;

g0[32,ind["R1"]]        = 1;
g1[32,ind["Exp_R1"]]    = 1;
Pi[32,3]                = 1;

g0[33,ind["R2"]]        = 1;
g1[33,ind["Exp_R2"]]    = 1;
Pi[33,4]                = 1;

g0[34,ind["K"]]         = 1;
g1[34,ind["Exp_K"]]     = 1;
Pi[34,5]                = 1;

g0[35,ind["G"]]         = 1;
g1[35,ind["Exp_G"]]     = 1;
Pi[35,6]                = 1;

## Observables 
# GDP observable
g0[36,ind["ΔGDPobs"]]   = 1;
g0[36,ind["Ygdp"]]      = -1;
g1[36,ind["Ygdp"]]      = -1;
g1[36,ind["G"]]         = 1;
constant[36]            = log(g);

# Consumption observable
g0[37,ind["ΔCobs"]]     = 1;
g0[37,ind["C"]]         = -1;
g1[37,ind["C"]]         = -1;
g1[37,ind["G"]]         = 1;
constant[37]            = log(g);

# Investment obervable
g0[38,ind["ΔINVobs"]]   = 1;
g0[38,ind["I"]]       = -1;
g1[38,ind["I"]]       = -1;
g1[38,ind["G"]]         = 1;
constant[38]            = log(g);

#tby_obs = tbaggout;                               % this is the empirical ratio
g0[39,ind["TBYobs"]]    = 1;
g0[39,ind["TBtotal"]]   = -1;
g0[39,ind["Ygdp"]]      = 1; 



G1, C0, G0, fmat, fwt, ywt, gev, eu = Gensys.gensysdt(g0, g1, constant, Psi, Pi);
#if eu ~= ones(2,1)
#    return
#end

# Observable equations for estimation
Nobs = 4
M1   = zeros(Nobs, NY);
M1[1,ind["ΔGDPobs"]] = 1
M1[2,ind["ΔCobs"]]   = 1
M1[3,ind["ΔINVobs"]] = 1
M1[4,ind["TBYobs"]]  = 1




return(G1, C0, G0, fmat, fwt, ywt, gev, eu, ind, NY, NEPS, M1)

end
