using CSV, Distributions, Stats
include("BBsteadystate.jl")
include("BBmodel.jl")
include("IRFplot.jl")
indir = "../posterior/fullSample/"


##----------------------------------------------------------------------------##
## Auxiliary Functions and Ordering ov variables 
##----------------------------------------------------------------------------##
## Function to merge estimating parameters with full parameter list
function paramVec(params,vectorOfParams,order)
    params_ret = copy(params)
    for vv in 1:length(vectorOfParams)
        params_ret[order[vv]] = vectorOfParams[vv]
    end
    return(params_ret)
end
#param =  ["ξ", "ψ", "ρ_p1", "ρ_p2", "σ_p",
#           "ρ_a", "ρ_a_til", "ρ_g", "ρ_s", "ρ_nu", "ρ_mu",
#           "σ_a", "σ_a_til", "σ_g", "σ_s", "σ_nu", "σ_mu"]

## Order of variables in vector of estimated parameters
order = Dict(1 => "ξ",    2 => "Ψ",
             3 => "ρ_p1", 4 => "ρ_p2",      5 => "σ_p",
             6 => "ρ_a",  7 => "ρ_a_til",   8 => "ρ_g",  9 => "ρ_s", 10 => "ρ_ν", 11 => "ρ_μ",             
             12 => "σ_a", 13 => "σ_a_til", 14 => "σ_g", 15 => "σ_s", 16 => "σ_ν", 17 => "σ_μ")

## Calibrated parameters
calibrated = Dict("p_til"   => 0.5244,
                  "dstar"   => -0.001,
                  "s"       => 0.0189,
                  "g"       => 1.0117204,
                  "αk"      => 0.32,
                  "αm"      => 0.05,
                  "αk_til"  => 0.32,
                  "δ"       => 0.1255,
                  "ϕ"       => 6.0,
                  "b"       => 0.9224,
                  "Γ"       => 2.0,
                  "θ"       => 1.6,
                  "ω"       => 1.6,
                  "ω_til"   => 1.6)

##----------------------------------------------------------------------------##
## Read in and merge posterior draws
##----------------------------------------------------------------------------##
postFiles = readdir("posterior")

##postFiles = ["post9.csv" , "post8.csv"]
draws = zeros(length(keys(order)),1)
for postFile in postFiles 
    tempFile = convert(Array,CSV.read(string(indir, postFile),header=false))
    draws    = [draws tempFile]
end
draws = draws[:,2:end]
ndraws = size(draws)[2]

##----------------------------------------------------------------------------##
## CI for IRF
##----------------------------------------------------------------------------##
# Solve model once, to get the dimensions and modal response
postMode = [mode(draws[ii,:]) for ii in 1:size(draws)[1]]
#[-0.199, 2.8, 0.95, 0.13, 0.1064,
#            0.9,0.9,0.9,0.9,0.9,0.9,
#            0.1,0.1,0.1,0.1,0.1,0.1]            
 #
T   = 10
ss  = BBsteadystate(paramVec(calibrated,postMode,order));
G1mode, C0_, G0mode, fmat_, fwt_, ywt_, gev_, eu_, ind_, NY_, NEPS_ = BBmodel(ss)
IRFmode = zeros(NY,NEPS,T);
for t=1:T
    IRFmode[:,:,t]=G1mode^(t-1)*G0mode;
end



#IRF = zeros(NY,NEPS,T,ndraws);
#for pp = 1:ndraws
IRF = SharedArray{Float64}(NY,NEPS,T,ndraws);
@parallel for pp = 1:ndraws
    par = paramVec(calibrated,draws[:,pp],order)
    ss  = BBsteadystate(par);
    G1, C0, G0, fmat, fwt, ywt, gev, eu, ind, NY, NEPS = BBmodel(ss)

    for t=1:T
        IRF[:,:,t,pp]=G1^(t-1)*G0;
    end
end
IRF = convert(Array, IRF)

shock     = ind["ϵ_P"]

## GDP
titles    = ["GDP"]
variables = [ind["Ygdp"]]
plotIRF(IRFmode*100, variables, shock, titles, "figures/figure4"; CI = true, IRFCI = IRF*100);
plotIRF(IRFmode*100, variables, shock, titles, "figures/figure4"; CI = true, IRFCI = IRF*100);
## Consumption

## Investment

## Trade Balance


