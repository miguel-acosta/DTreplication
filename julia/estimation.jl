using Distributions
using Optim
using IterableTables, DataFrames, ExcelReaders, PyPlot, Parameters, CSV
include("kalmanLikelihood.jl")
include("BBmodel.jl")
include("BBsteadystate.jl")
##----------------------------------------------------------------------------##
## Auxiliary Functions
##----------------------------------------------------------------------------##

## Function for stability of AR(2) coefficients:
function stableAR2(ϕ1::Float64, ϕ2::Float64)
    return( (ϕ2<1+ϕ1) & (ϕ2<1-ϕ1) & (ϕ2>-1) )
end

## Function to parameterize beta pdf with mean and variance
## (thanks stack exchange for the formulas)
function BetaMV(μ::Float64,σ::Float64)
    α = ((1-μ)/(σ^2) - 1/μ) * μ^2
    β = α*(1/μ - 1)
    return(Beta(α,β))
end

## Function to parameterize inverse gamma  pdf with mean and variance
## (thanks wikipedia and a pencil)
function InverseGammaMV(μ::Float64,σ::Float64)
    α = (μ^2)/(σ^2) + 2
    β = μ*(α-1)
    return(InverseGamma(α,β))
end

## Function to merge estimating parameters with full parameter list
function paramVec(params::Dict{String,Float64},vectorOfParams::Array{Float64,1},order::Dict{Int64,String})
    params_ret = copy(params)
    for vv in 1:length(vectorOfParams)
        params_ret[order[vv]] = vectorOfParams[vv]
    end
    return(params_ret)
end

## Function to turn dictionary of parameters into vector
## based on some pre-defined order 
function estParamVec(dictOfParams::Dict{String,Float64},order::Dict{Int64,String})
    np  = length(keys(dictOfParams))
    VEC = zeros(np,1)
    orderRev = map(reverse,order)
    for pp in keys(dictOfParams)
        VEC[orderRev[pp]] = dictOfParams[pp]
    end
    return(vec(VEC))
end

function gradientFD(f,x;step = 1e-5,index=0)
    if index > 0
        x_copy        = copy(x)
        x_copy[index] = x_copy[index] + step
        grad          = (f(x_copy) - f(x))/step
    else 
        grad = zeros(length(x),1)
        for ii = 1:length(x)
            x_copy     = copy(x)
            x_copy[ii] = x_copy[ii] + step
            grad[ii]   = (f(x_copy) - f(x))/step
        end
    end
    return(grad)
end

## Inverse hessian via finite differencing
function hessianFD(f,x;step= 1e-5)
    N    = length(x)
    hess = zeros(N,N)
    gx = gradientFD(f,x)
    for ii = 1:N
        for jj = 1:N
            x_move      = copy(x)
            x_move[jj]  = x_move[jj] + step
            hess[ii,jj] = (gradientFD(f,x_move;index=ii)-gx[ii])/step
        end
    end
    return(hess)
end

##----------------------------------------------------------------------------##
## Set the prior 
##----------------------------------------------------------------------------##
function logprior(PARAMS::Dict{String,Float64})
    # This function evaluates the log prior given a vector of parameter
    # values by using the predetermined distributions.

    # Input : a (17,1) vector of 17 model parameters 
    # Output: logP = log joint prior of 17 parameters
    # NOTE  : logP = -inf indicates the support error. 


    # Check if each parameber is NOT in its support
    #     (-1,1) for AR(1) coefficients
    #     [0,Inf) for standard deviations
    #     stable AR(2) coefficients according to the magic triangle
    unstableAR1 = any([(PARAMS[pp] < -1) | (PARAMS[pp] > 1) for pp in ["ρ_a", "ρ_a_til", "ρ_g", "ρ_s", "ρ_ν", "ρ_μ"]])
    negStdev    = any([PARAMS[pp] < 0                  for pp in ["σ_a", "σ_a_til", "σ_g", "σ_s", "σ_ν", "σ_μ", "σ_p"]])
    unstableAR2 = !stableAR2(PARAMS["ρ_p1"], PARAMS["ρ_p2"])
    
    if unstableAR1 | negStdev | unstableAR2
        logP = -Inf
    else
        Prob = [pdf(Normal(-0.199,0.045)  , PARAMS["ξ"]),       
                pdf(Normal(2.8,0.5)       , PARAMS["Ψ"]),       
                pdf(BetaMV(0.8,0.2)       , PARAMS["ρ_p1"]),    
                pdf(BetaMV(0.15,0.1)      , PARAMS["ρ_p2"]),              
                pdf(InverseGammaMV(0.05,2.0), PARAMS["σ_p"]),     
                pdf(BetaMV(0.5,0.2)       , PARAMS["ρ_a"]),     
                pdf(BetaMV(0.5,0.2)       , PARAMS["ρ_a_til"]), 
                pdf(BetaMV(0.5,0.2)       , PARAMS["ρ_g"]),     
                pdf(BetaMV(0.5,0.2)       , PARAMS["ρ_s"]),     
                pdf(BetaMV(0.5,0.2)       , PARAMS["ρ_ν"]),     
                pdf(BetaMV(0.5,0.2)       , PARAMS["ρ_μ"]),              
                pdf(InverseGammaMV(0.05,2.0), PARAMS["σ_a"]),     
                pdf(InverseGammaMV(0.05,2.0), PARAMS["σ_a_til"]), 
                pdf(InverseGammaMV(0.05,2.0), PARAMS["σ_g"]),     
                pdf(InverseGammaMV(0.05,2.0), PARAMS["σ_s"]),     
                pdf(InverseGammaMV(0.05,2.0), PARAMS["σ_ν"]),     
                pdf(InverseGammaMV(0.05,2.0), PARAMS["σ_μ"])]
        logP = sum(log.(Prob));
    end
    return(logP)
end

# Derive the posterior mode and save for future uses.
function logposterior(calibrated::Dict{String,Float64},p::Array{Float64,1},order::Dict{Int64,String},DATA::Array{Float64,2})
    prior = logprior(paramVec(calibrated,p,order))
    if prior == -Inf
        return(-Inf)
    else
        return(prior + kalmanLikelihood(paramVec(calibrated,p,order),DATA))
    end
end

##----------------------------------------------------------------------------##
## Setup: set directories and import data 
##----------------------------------------------------------------------------##
function estimate(short::Bool=false)
    srand(4)
    if contains(pwd(), "jma2241") # Tells us whether we're on the cluster
        outdirMCMC = "/rigel/sscc/users/jma2241/DT/"
    else
        outdirMCMC = "../posterior/"
    end

    ## Import data
    if !short 
        startDat = 1
        outdirMCMC = string(outdirMCMC, "posterior_full/")
        modeSuffix = ""
    else
        startDat = 50
        outdirMCMC = string(outdirMCMC, "posterior_short/")
        #modeSuffix = "_short"
        modeSuffix = ""
    end
    dat     = convert(Array{Float64,2},readxlsheet("../VAR/DataVAR.xlsx", "Sheet1", skipstartrows=1))
    datVAR  = dat[:,2:5]
    DATA    = [datVAR[2:end,1:3] - datVAR[1:(end-1),1:3] datVAR[2:end,4]].' #DATA(:,t) refers to period t's observations.
    DATA    = DATA[:,startDat:end] # short sample
    
    ##----------------------------------------------------------------------------##
    ## Calibrate parameters and establish order of parameters in vector
    ##----------------------------------------------------------------------------##
    
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
    ## Use minimizer to find (or, hopefully, get closer to) posterior mode.
    ##----------------------------------------------------------------------------##
    
    neglogposterior(p::Array{Float64,1}) = -logposterior(calibrated,p,order,DATA)::Float64
    logposterior_calib(p::Array{Float64,1}) = logposterior(calibrated,p,order,DATA)::Float64
    # Prior Mean of Theta
    θ0 = Dict("ξ"        => -0.199,       
              "Ψ"        => 2.8,       
              "ρ_p1"     => 0.8,    
              "ρ_p2"     => 0.15,    
              "σ_p"      => 0.05,     
              "ρ_a"      => 0.5,     
              "ρ_a_til"  => 0.5, 
              "ρ_g"      => 0.5,     
              "ρ_s"      => 0.5,     
              "ρ_ν"      => 0.5,     
              "ρ_μ"      => 0.5,     
              "σ_a"      => 0.05,     
              "σ_a_til"  => 0.05, 
              "σ_g"      => 0.05,     
              "σ_s"      => 0.05,     
              "σ_ν"      => 0.05,     
              "σ_μ"      => 0.05)
    θ0vec = estParamVec(θ0,order)
    
    
    ## Initialize the MCMC.
    NP     = length(θ0);     # number of parameters to estimate
    findMode = false
    if findMode == true
        result = optimize(neglogposterior,θ0vec,LBFGS(), Optim.Options(iterations=200,show_trace=true,x_tol=1e-16,f_tol=1e-16,allow_f_increases=true, g_tol=1e-16))
        Theta  = Optim.minimizer(result)
        CSV.write(string("mode",modeSuffix,".csv"), DataFrame(Theta))
        H = inv(hessianFD(neglogposterior,Theta))
        CSV.write(string("H",modeSuffix,".csv"), DataFrame(H))
    end
    Theta  = vec(convert(Array, CSV.read(string("mode",modeSuffix,".csv"),allowmissing=:none)))
    H      = Symmetric(convert(Array, CSV.read(string("H",modeSuffix,".csv"),allowmissing=:none)))
     
    ##----------------------------------------------------------------------------##
    ## MCMC!
    ##----------------------------------------------------------------------------##
    scale = 0.2;
    VarR  = scale*H;     ## Variance of the random walk. 
    N     = 500_000;     ## number of iterations 
    Nsave = 10_000;
    θ0    = vec(rand(MvNormal(Theta, VarR),1))
    
    while logposterior_calib(θ0) == -Inf
        θ0   = vec(rand(MvNormal(Theta, VarR),1))
    end
    
    POST = zeros(NP, N) ## POST(;,t) refers to the posterior of time t
    naccept = 0
    #tic()    
    t = 1
    while t <= N
        θ1::Array{Float64,1} = vec(θ0+ rand(MvNormal(vec(zeros(NP,1)), VarR),1))
        post1::Float64       = logposterior_calib(θ1)
        if post1 == -Inf
        else
            post0 = logposterior_calib(θ0)
            alpha = minimum([1.0, exp(post1-post0)])
            if rand(Uniform()) < alpha
                θ0        = copy(θ1)
                naccept += 1
            end
            POST[:,t] = copy(θ0)
            t += 1
            if t % Nsave == 0
                CSV.write(string(outdirMCMC, "post", Int(t/Nsave) , ".csv"), DataFrame(POST[:,(t-Nsave+1):t]), header = false)
                println(string("You are at iteration ", t, " of ", N, ". Acceptance rate is: ", naccept/t))
                #toc()
                #tic()
            end
        end

    end

end

estimate(true) ## True for short sample, false for full sample
