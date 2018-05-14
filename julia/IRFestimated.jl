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
## Write a nice table
##----------------------------------------------------------------------------##
function writeParamTable(fname, rows, centers, lows, highs;precision = 2)
    nlines = size(centers)[1]
    nspecs = size(centers)[2]
    open(string("figures/", fname, ".tex"),"w") do f
        for nn in 1:nlines
            rr = rows[nn]
            write(f,"$rr&")
            for ss in 1:nspecs
                cc = round(centers[nn,ss],2)
                write(f, "$cc")
                ss < nspecs ? write(f, "&") : write(f, "\\\\\n")
                
            end
            write(f,"&")
            for ss in 1:nspecs
                ll = isnan(lows[nn,ss]) ? "--" : round(lows[nn,ss],precision)
                hh = isnan(lows[nn,ss]) ? "--" : round(highs[nn,ss],precision)
                write(f, "[$ll, $hh]")
                ss < nspecs ? write(f, "&") : write(f, "\\\\\n")
            end
        end
    end
end
 

##----------------------------------------------------------------------------##
## Read in and merge posterior draws
##----------------------------------------------------------------------------##
function getDraws(DIR, drawsPerFile, nfiles, order, burnin, thin; mergeDraws = true)
    if mergeDraws
        ndraws = drawsPerFile * nfiles; 
        draws = SharedArray{Float64}(length(keys(order)),ndraws)
        @time @sync@parallel for pp in 1:nfiles
            tempFile = readcsv(string(DIR, "post", pp, ".csv"))
            draws[:,((pp-1)*drawsPerFile+1):(pp*drawsPerFile)] = tempFile[:,1:drawsPerFile]
        end
        draws = convert(Array,draws[:,Int(floor(ndraws*burnin)):thin:ndraws])
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
function getIRFat(meanORmode, draws, order, calibrated, T)
    # Solve model once, to get the dimensions and modal response
    if meanORmode == "mode"
        postMode = [modeBin(draws[ii,:]) for ii in 1:size(draws)[1]]
    elseif meanORmode == "mean"
        postMode = mean(draws,2)
    end
        
    calibratedParams = [-0.199, 2.8, 0.95, 0.13, 0.1064,
                        0.9,0.9,0.9,0.9,0.9,0.9,
                        0.1,0.1,0.1,0.1,0.1,0.1]

    ss  = BBsteadystate(paramVec(calibrated,postMode,order));
    G1mode, C0_, G0mode, fmat_, fwt_, ywt_, gev_, eu_, ind, NY, NEPS = BBmodel(ss)
    IRFmode = zeros(NY,NEPS,T);
    for t=1:T
        IRFmode[:,:,t]=G1mode^(t-1)*G0mode;
    end
    decomp = varDecomp(G0mode, G1mode, NY,ind)
    
    return(IRFmode, decomp, NY, NEPS, ind) 
end

function varDecomp(G0::Array{Float64,2}, G1::Array{Float64,2}, NY, ind)
    #vecΣy = inv(eye(NY^2) - kron(G1,G1)) * vec(G0 * G0.')
    #Σy    = reshape(vecΣy, NY,NY)
    Φ0::Array{Float64,2}    = copy(G0)
    Σy0::Array{Float64,2}   = Φ0 * (Φ0.')
    DC0::Array{Float64,2}   = Φ0.^2

    Φ1    = Array{Float64,2}
    DC1   = Array{Float64,2}
    Σy1   = Array{Float64,2}    
    
    diff::Float64  = 1.0;
    tol::Float64   = 1e-6;
    while diff > tol
        Φ1    = G1*Φ0
        Σy1   = Σy0+Φ1*(Φ1.')        
        DC1   = DC0+(Φ1.^2)
        diff  = maximum([maximum(abs.(DC1-DC0)), maximum(abs.(Σy1-Σy0))])
        #diff  = maximum(abs.(DC1-DC0))
        # Update for the next round of loop
        DC0   = copy(DC1)
        Φ0    = copy(Φ1)
        Σy0   = copy(Σy1)
    end
    # Final outcome for the unconditional variance
    DC::Array{Float64,2}   = copy(DC1)
    Σy::Array{Float64,2}   = copy(Σy1)
    
    # Show the percentage contribution of each shock

    # Forecast error variance
    # Each row of DC: the percentage contribution of that variable decomposed
    # into all present shoks
    for n=1:NY
        DC[n,:] = DC[n,:]/Σy[n,n]*100;
        #DC[n,:] = DC[n,:]/sum(DC[n,:])*100;
    end


    Table::Array{Float64,2} = zeros(4,6);
    
    indsOfInterest::Array{Int64,1} = [ind["ΔGDPobs"], ind["ΔCobs"],  ind["ΔINVobs"],  ind["TBYobs"]]
    
    obs = DC[indsOfInterest,:]; 
    Table[:,1] = obs[:,ind["ϵ_A"]]+obs[:,ind["ϵ_A_til"]]
    Table[:,2] = obs[:,ind["ϵ_G"]]
    Table[:,3] = obs[:,ind["ϵ_μ"]]
    Table[:,4] = obs[:,ind["ϵ_P"]]
    Table[:,5] = obs[:,ind["ϵ_S"]]
    Table[:,6] = obs[:,ind["ϵ_ν"]]
    return(Table)
end


function IRFdist(draws, NY, ind, calibrated, order, shock, T)
    ndraws    = size(draws)[2]
    shock     = ind["ϵ_P"]
    IRF       = SharedArray{Float64}(NY,T,ndraws) # zeros(NY,T,ndraws) #
    decomp    = SharedArray{Float64}(4,6,ndraws)  # zeros(NY,T,ndraws) #
    @time @sync @parallel for pp = 1:ndraws
    #for pp = 1:ndraws
        par::Dict{String,Float64} = paramVec(calibrated,draws[:,pp],order)
        ss::Dict{String,Float64}  = BBsteadystate(par);
        G1, C0, G0 = BBmodel(ss)

        for t=1:T
            irftt = G1^(t-1)*G0;
            IRF[:,t,pp] = irftt[:,shock]
        end
        decomp[:,:,pp] = varDecomp(G0, G1, NY,ind)
        #if pp % 1000 == 0
        #    println(string("Draw ", pp, " of " ,ndraws))
        #end
    end
    IRF = convert(Array, IRF)
    decomp = convert(Array, decomp)
    return(IRF,decomp)
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
function getMedianIRF(IRFall, inds, NY, T)
    medIRF = zeros(NY,T);
    for tt = 1:T
        for ny = 1:length(inds)
            medIRF[inds[ny],tt] = median(IRFall[inds[ny],tt,:])
        end
    end
    return(medIRF)
end

function makeIRFplots(pltdir; thinning = 10_000, burnin  = 0.25, drawsPerFile = 9_999, nfiles = 0, T = 10)
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

    draws, ndraws = getDraws(indir, drawsPerFile, nfiles, order, burnin, thinning)
    IRFmode, decompMode, NY, NEPS, ind = getIRFat("mode", draws, order, calibrated, T)
    IRFmean, decompMean, NY, NEPS, ind = getIRFat("mean", draws, order, calibrated, T)    
    
    shock   = ind["ϵ_P"]
    indsAll = [ind["Ygdp"], ind["C"], ind["I"], ind["TBYobs"]]

    IRF,decomp = IRFdist(draws, NY, ind, calibrated, order, shock, T)

    ## Variance Decomposition Table
    nobs = size(decompMode)[1]
    nexo = size(decompMode)[2]
    decompLOW  = zeros(nobs,nexo)
    decompHI   = zeros(nobs,nexo)
    Clevel = 90
    for ii = 1:nobs
        decompLOW[ii,:]  = [percentile(vec(decomp[ii,nn,:]),Int((100-Clevel)/2)) for nn in 1:nexo]
        decompHI[ii,:]   = [percentile(vec(decomp[ii,nn,:]),Int(100-(100-Clevel)/2)) for nn in 1:nexo]
    end

    writeParamTable(string(pltdir, "decomp"),
                    ["\$\\Delta \\text{GDP}_t\$",
                     "\$\\Delta C_t\$",
                     "\$\\Delta I_t\$",
                     "\$ \\text{TB}_t/\\text{Y}_t\$"],
                    decompMode, decompLOW, decompHI)
    writeParamTable(string(pltdir, "decomp_mean"),
                    ["\$\\Delta \\text{GDP}_t\$",
                     "\$\\Delta C_t\$",
                     "\$\\Delta I_t\$",
                     "\$ \\text{TB}_t/\\text{Y}_t\$"],
                    decompMean, decompLOW, decompHI)
    

    
    modalIRF   = getModalIRF(IRF, indsAll, NY, T)    
    PyPlot.plt[:hist](IRF[ind["Ygdp"], 7,:],100, color = "dodgerblue")
    savefig(string("figures/", pltdir, "dist_gdp_7.pdf"))

    ##--------------------------------------------- IRF at MODE -----------------------------------##
    linelabs = ["IRF at Mode", "Modal IRF", "IRF at Mean"]
    ## GDP
    plotGDP = [IRFmode[ind["Ygdp"],shock,:] modalIRF[ind["Ygdp"],:] IRFmean[ind["Ygdp"],shock,:]]
    titles    = ["GDP"]
    variables = [ind["Ygdp"]]
    plotIRF(plotGDP*100, variables, shock, titles, string("figures/",pltdir,"GDP");
            CI = true, IRFCI = IRF*100, LINELABEL = linelabs)
    plotIRF(plotGDP*100, variables, shock, titles, string("figures/",pltdir,"GDP");
            CI = true, IRFCI = IRF*100, LINELABEL = linelabs, TITLE = "GDP")
    
    ## Consumption
    plotC = [IRFmode[ind["C"],shock,:] modalIRF[ind["C"],:] IRFmean[ind["C"],shock,:]]    
    titles    = ["Consumption"]
    variables = [ind["C"]]
    plotIRF(plotC*100, variables, shock, titles, string("figures/",pltdir,"C");
            CI = true, IRFCI = IRF*100, LINELABEL = linelabs, TITLE = "Consumption")
    
    ## Investment
    plotI = [IRFmode[ind["I"],shock,:] modalIRF[ind["I"],:] IRFmean[ind["I"],shock,:]]        
    titles    = ["Investment"]
    variables = [ind["I"]]
    plotIRF(plotI*100, variables, shock, titles, string("figures/",pltdir,"I");
            CI = true, IRFCI = IRF*100, LINELABEL = linelabs, TITLE = "Investment")
    
    ## Trade Balance
    plotTBY   = [IRFmode[ind["TBYobs"],shock,:] modalIRF[ind["TBYobs"],:] IRFmean[ind["TBYobs"],shock,:]]        
    titles    = ["Trade balance/GDP"]
    variables = [ind["TBYobs"]]
    plotIRF(plotTBY*100, variables, shock, titles, string("figures/",pltdir,"TB");
            CI = true, IRFCI = IRF*100, LINELABEL = linelabs, TITLE = "Trade balance/GDP")

    
end
makeIRFplots("posterior_full_20180511/" ; thinning = 10, burnin = 0.25)
makeIRFplots("posterior_flat_20180511/" ; thinning = 10, burnin = 0.25)
makeIRFplots("posterior_short_20180511/"; thinning = 10, burnin = 0.25)

