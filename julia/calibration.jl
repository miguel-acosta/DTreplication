include("BBsteadystate.jl")
include("BBmodel.jl")
include("IRFplot.jl")

# Structural parameters
par = Dict("p_til"   => 0.5244,
           "dstar"   => -0.001,
           "s"       => 0.0189,
           "ξ"       => -0.199,
           "g"       => 1.0117204,
           "Ψ"       => 2.8,
           "αk"      => 0.32,
           "αm"      => 0.05,
           "αk_til"  => 0.32,
           "δ"       => 0.1255,
           "ϕ"       => 6.0,
           "b"       => 0.9224,
           "Γ"       => 2.0,
           "θ"       => 1.6,
           "ω"       => 1.6,
           "ω_til"   => 1.6,
           "ρ_a"     => 0.9,
           "ρ_a_til" => 0.9,
           "ρ_g"     => 0.9,
           "ρ_s"     => 0.9,
           "ρ_ν"     => 0.9,
           "ρ_μ"     => 0.9,
           "σ_a"     => 0.1,
           "σ_a_til" => 0.1,
           "σ_g"     => 0.1,
           "σ_s"     => 0.1,
           "σ_ν"     => 0.1,
           "σ_μ"     => 0.1,
           "ρ_p1"    => 0.95,
           "ρ_p2"    => 0.13,
           "σ_p"     => 0.1064)

# Get the steady state
ss = BBsteadystate(par);
#dt = load('Crosscheck/Solution files/DTsteady.mat');

# Get the solution 
G1, C0, G0, fmat, fwt, ywt, gev, eu, ind, NY, NEPS = BBmodel(ss)


 # Calculate Impulses Responses 
T = 10;
 
IRF = zeros(NY,NEPS,T);
for t=1:T
    IRF[:,:,t]=G1^(t-1)*G0;
end
 
 # Figure 4
shock     = ind["ϵ_P"]
titles    = ["GDP", "Consumption", "Investment", "Trade balance/GDP"]
variables = [ind["Ygdp"], ind["C"], ind["I"], ind["TBYobs"]]

plotIRF(IRF[:,shock,:]*100, variables, shock, titles, "figures/figure4")
