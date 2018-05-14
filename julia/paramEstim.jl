using Stats

##----------------------------------------------------------------------------##
## Options
##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##
## Auxiliary Functions and Ordering ov variables 
##----------------------------------------------------------------------------##
function modeBin(vv::Array{Float64}; nbin::Int64=1000)
    bins::Array{Float64}   = (minimum(vv)-1e-8):((maximum(vv)-minimum(vv))/nbin):maximum(vv)
    nbin::Int64            = length(bins)
    binMed::Array{Float64} = zeros(nbin-1,1)
    binVal::Array{Float64} = zeros(nbin-1,1)
    J::Int64 = length(vv)
    for ii = 2:nbin
        binMed[ii-1] = (bins[ii]+bins[ii-1])/2
        binVal[ii-1] = sum((vv .> bins[ii-1]) .* (vv .<= bins[ii]))
    end
    mxval,mxind = findmax(binVal)
    return(binMed[mxind[1]])
end

##----------------------------------------------------------------------------##
## Read in and merge posterior draws
##----------------------------------------------------------------------------##
function getDraws(DIR, nparam; drawsPerFile = 9_999, burnin = 0.25, thin = 1)
    nfiles  = length(readdir(DIR))
    nfiles  = in("postAll.csv", readdir(DIR)) ? nfiles -1 : nfiles
    
    ndraws = drawsPerFile * nfiles; 
    draws = SharedArray{Float64}(nparam,ndraws)
    @time @sync@parallel for pp in 1:nfiles
        tempFile = readcsv(string(DIR, "post", pp, ".csv"))
        draws[:,((pp-1)*drawsPerFile+1):(pp*drawsPerFile)] = tempFile[:,1:drawsPerFile]
    end
        draws = convert(Array,draws[:,Int(floor(ndraws*burnin)):thin:ndraws])
    return(draws)
end

##----------------------------------------------------------------------------##
## CI for IRF
##----------------------------------------------------------------------------##
function getStats(draws; Clevel = 90)
    K = size(draws)[1]
    postMode = [modeBin(draws[ii,:]) for ii in 1:K]
    postMean = mean(draws,2)
    CIlow  = [percentile(vec(draws[kk,:]),Int((100-Clevel)/2)) for kk in 1:K]
    CIhigh = [percentile(vec(draws[kk,:]),Int(100-(100-Clevel)/2)) for kk in 1:K]
    return(postMode, postMean, CIlow, CIhigh)
end

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


#calib = [-0.199,2.8,0.95,0.13,0.1064,0.9,0.9,0.9,0.9,0.9,0.9,0.1,0.1,0.1,0.1,0.1,0.1]
prior = [-0.199,2.8,0.8,0.15,0.10,0.5,0.5,0.5,0.5,0.5,0.5,0.10,0.10,0.10,0.10,0.10,0.10]
priorLow  = repmat([NaN], length(prior))
priorHigh = repmat([NaN], length(prior))

paperMean = [-0.2212,3.2057,0.8060,0.1278,0.1765,0.8277,0.5887,0.5244,
             0.6440,0.8687,0.9199,0.0295,0.0525,0.0261,0.1876,0.4582,0.0547]
paperLow  = [-0.2876,2.5050,0.6840,0.0105,0.0876,0.7494,0.2827,0.3199,
            0.5075,0.8382,0.8743,0.0231,0.0242,0.0193,0.1659,0.4145,0.0410]
paperHigh = [-0.1550,3.8984,0.9388,0.2298,0.2652,0.9092,0.8980,0.7299,
             0.7832,0.8996,0.9693,0.0360,0.0810,0.0327,0.2089,0.5000,0.0683]

names = ["\$\\xi\$", "\$\\psi\$",
         "\$\\rho_{\\tilde p}^1\$", "\$-\\rho_{\\tilde p}^1\$", "\$\\sigma_{\\tilde p}\$",
         "\$\\rho_a\$","\$\\rho_{\\tilde a}\$","\$\\rho_g\$",
         "\$\\rho_s\$","\$\\rho_\\nu\$","\$\\rho_\\mu\$",
         "\$\\sigma_a\$","\$\\sigma_{\\tilde a}\$","\$\\sigma_g\$",
         "\$\\sigma_s\$","\$\\sigma_\\nu\$","\$\\sigma_\\mu\$"]
    

K          = length(paperMean)
drawsFull  = getDraws("../posterior/posterior_full_20180511/", K)
drawsShort = getDraws("../posterior/posterior_short_20180511/", K)
drawsFlat  = getDraws("../posterior/posterior_flat_20180511/", K)

modeFull,  meanFull, lowFull, highFull = getStats(drawsFull)
modeShort, meanShort, lowShort, highShort = getStats(drawsShort)
modeFlat,  meanFlat, lowFlat, highFlat = getStats(drawsFlat)

writeParamTable("params",names,
                [prior paperMean      meanShort meanFull meanFlat ],
                [priorLow paperLow    lowShort  lowFull  lowFlat  ],
                [priorHigh paperHigh  highShort highFull highFlat ])



writeParamTable("params_modes",names,
                [prior paperMean      modeShort modeFull modeFlat ],
                [priorLow paperLow    lowShort  lowFull  lowFlat  ],
                [priorHigh paperHigh  highShort highFull highFlat ])
