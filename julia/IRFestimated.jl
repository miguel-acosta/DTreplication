using CSV, Distributions, Stats, DataFrames
include("BBsteadystate.jl")
include("BBmodel.jl")
include("IRFplot.jl")

##----------------------------------------------------------------------------##
## Options
##----------------------------------------------------------------------------##

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

function modeBin(vv::Array{Float64}; nbin::Int64=1000)
    bins::Array{Float64}   = (minimum(vv)-1e-8):((maximum(vv)-minimum(vv))/nbin):maximum(vv)
    nbin::Int64            = length(bins)
    binMed::Array{Float64} = zeros(nbin-1,1)
    binVal::Array{Float64} = zeros(nbin-1,1)
    J::Int64 = length(vv)
    for ii = 2:nbin
        binMed[ii-1] = (bins[ii]+bins[ii-1])/2
        binVal[ii-1] = sum((vv .> bins[ii-1]) .* (vv .<= bins[ii]))
#        binVal[ii-1] = sum([(vv[jj] > bins[ii-1]) & (vv[jj] <= bins[ii-1]) for jj in 1:J])
    end
    mxval,mxind = findmax(binVal)
    return(binMed[mxind[1]])
end
#param =  ["ξ", "ψ", "ρ_p1", "ρ_p2", "σ_p",
#           "ρ_a", "ρ_a_til", "ρ_g", "ρ_s", "ρ_nu", "ρ_mu",
#           "σ_a", "σ_a_til", "σ_g", "σ_s", "σ_nu", "σ_mu"]


##----------------------------------------------------------------------------##
## Read in and merge posterior draws
##----------------------------------------------------------------------------##
function getDraws(DIR, drawsPerFile, nfiles, order, burnin; mergeDraws = true)
    if mergeDraws
        ndraws = drawsPerFile * nfiles; 
        draws = SharedArray{Float64}(length(keys(order)),ndraws)
        @time @sync@parallel for pp in 1:nfiles
            tempFile = readcsv(string(DIR, "post", pp, ".csv"))
            draws[:,((pp-1)*drawsPerFile+1):(pp*drawsPerFile)] = tempFile
        end
        draws = convert(Array,draws[:,Int(ndraws*burnin):ndraws])
        writecsv(string(DIR, "postAll.csv"), draws)
    else
        draws = readcsv(string(DIR, "postAll.csv"))
    end
    ndraws = size(draws)[2]
    return(draws, ndraws)
end

##----------------------------------------------------------------------------##
## CI for IRF
##----------------------------------------------------------------------------##
function getIRFatMode(draws, order, calibrated, T)
    # Solve model once, to get the dimensions and modal response
    postMode = [modeBin(draws[ii,:]) for ii in 1:size(draws)[1]]
    calibratedParams = [-0.199, 2.8, 0.95, 0.13, 0.1064,
                        0.9,0.9,0.9,0.9,0.9,0.9,
                        0.1,0.1,0.1,0.1,0.1,0.1]

    ss  = BBsteadystate(paramVec(calibrated,postMode,order));
    G1mode, C0_, G0mode, fmat_, fwt_, ywt_, gev_, eu_, ind, NY, NEPS = BBmodel(ss)
    IRFmode = zeros(NY,NEPS,T);
    for t=1:T
        IRFmode[:,:,t]=G1mode^(t-1)*G0mode;
    end
    return(IRFmode, NY, NEPS, ind) 
end


function IRFdist(draws, NY, ind, calibrated, order, shock, T)
    ndraws    = size(draws)[2]
    shock     = ind["ϵ_P"]
    IRF       = SharedArray{Float64}(NY,T,ndraws) # zeros(NY,T,ndraws) #
    @time @sync @parallel for pp = 1:ndraws
    #for pp = 1:ndraws
        par::Dict{String,Float64} = paramVec(calibrated,draws[:,pp],order)
        ss::Dict{String,Float64}  = BBsteadystate(par);
        G1, C0, G0 = BBmodel(ss)

        for t=1:T
            irftt = G1^(t-1)*G0;
            IRF[:,t,pp] = irftt[:,shock]
        end
#        if pp % 1000 == 0
#            println(string("Draw ", pp, " of " ,ndraws))
#        end
    end
    IRF = convert(Array, IRF)
    return(IRF)
end

function getModalIRF(IRFall, inds, NY, T)
    modalIRF = zeros(NY,T);
    for tt = 1:T
        for ny = 1:length(inds)
            modalIRF[inds[ny],tt] = modeBin(IRFall[inds[ny],tt,:];nbin=100)
        end
    end
    return(modalIRF)
end


function makeIRFplots(pltdir; burnin  = 0.25, drawsPerFile = 10_000, nfiles = 0, T = 10)
    indir   = string("../posterior/", pltdir)
    if nfiles == 0
        nfiles  = length(readdir(indir))
        nfiles  = in("postAll.csv", readdir(indir)) ? nfiles -1 : nfiles
    end

    
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

    draws, ndraws = getDraws(indir, drawsPerFile, nfiles, order, burnin)
    IRFmode, NY, NEPS, ind = getIRFatMode(draws, order, calibrated, T)
    
    
    shock   = ind["ϵ_P"]
    indsAll = [ind["Ygdp"], ind["C"], ind["I"], ind["TBYobs"]]

    IRF = IRFdist(draws, NY, ind, calibrated, order, shock, T)

    modalIRF = getModalIRF(IRF, indsAll, NY, T)

    IRFforLine = copy(IRFmode[:,shock,:]) # or modalIRF

    ## GDP
    titles    = ["GDP"]
    variables = [ind["Ygdp"]]
    plotIRF(IRFforLine*100, variables, shock, titles, string("figures/",pltdir,"GDP"); CI = true, IRFCI = IRF*100)
    plotIRF(IRFforLine*100, variables, shock, titles, string("figures/",pltdir,"GDP"); CI = true, IRFCI = IRF*100, LINELABEL = "IRF at Posterior Mode", TITLE = "GDP")
    
    ## Consumption
    titles    = ["Consumption"]
    variables = [ind["C"]]
    plotIRF(IRFforLine*100, variables, shock, titles, string("figures/",pltdir,"C"); CI = true, IRFCI = IRF*100, LINELABEL = "IRF at Posterior Mode", TITLE = "Consumption")
    
    ## Investment
    titles    = ["Investment"]
    variables = [ind["I"]]
    plotIRF(IRFforLine*100, variables, shock, titles, string("figures/",pltdir,"I"); CI = true, IRFCI = IRF*100, LINELABEL = "IRF at Posterior Mode", TITLE = "Investment")
    
    ## Trade Balance
    titles    = ["Trade balance/GDP"]
    variables = [ind["TBYobs"]]
    plotIRF(IRFforLine*100, variables, shock, titles, string("figures/",pltdir,"TB"); CI = true, IRFCI = IRF*100, LINELABEL = "IRF at Posterior Mode", TITLE = "Trade balance/GDP")
    
end



makeIRFplots("posterior_short_20180511/") #; nfiles = 1)
