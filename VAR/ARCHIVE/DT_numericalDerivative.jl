include("VARfuncs.jl")
using IterableTables, DataFrames, ExcelReaders, PyPlot, Distributions, Stats

dat = convert(Array{Float64,2},readxlsheet("DataVAR.xlsx", "Sheet1", skipstartrows=1))

datVAR = dat[:,2:end]
datVAR = datVAR[:,end:-1:1]
TT     = size(dat)[1]
P      = 2
N      = size(datVAR)[2]
trend  = dat[:,1]

tIRF = 10
Ndraws = 9999
Clevel = 100
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##--## ## ESTIMATE THE VAR
Aolsct,μolsct,Σolsct,residsct = VARols(P, datVAR; cons=true,  Z = trend)

Pchol = chol(Σolsct).'
shock = Pchol*[1, 0, 0, 0, 0]
responses = IRF(Aolsct, 0*μolsct[:,1], tIRF-1, shock)'


##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Confidence Intervals

## Stack estimates
Π        = [μolsct reshape(Aolsct, size(Aolsct)[1], size(Aolsct)[2]*2)]
π        = vec(Π.')

## Build Q
dat_lags = lagmatrix(datVAR,2)
Z        = dat_lags[:,N+1:end]
X        = [ones(TT-P,1) trend[(P+1):end] Z] ## the i-th row of this is x_t from hamilton
Q        = (X.'*X)/(TT-P)

## Build G
K = N*P+2
G = zeros(N, N*K, tIRF)
Δ = 1e-4
for ee = 1:N*K
    ## Add a little to the parameter
    e     = zeros(length(π),1)
    Δeff  = ee % K != 2 ? copy(Δ) : 0.0 ## the trend is weird... 
    e[ee] = copy(Δeff)
    π_e   = π + e

    ## Reshape the parameter
    π_e_unpack     = reshape(π_e, N*P+2,N)
    μ_MC           = π_e_unpack[1:2,:].'
    A_MC           = reshape(π_e_unpack[3:end,:].', N,N,2)
    Π_MC           = [μ_MC reshape(A_MC, size(A_MC)[1], size(A_MC)[2]*2)]

    ## Re-estimate IRF
    A_MC,μ_MC,Σ_MC = VARols(2, datVAR; cons=true,  Z = trend, βfixed = Π_MC)
    Pchol_MC       = chol(Σ_MC).'
    shock_MC       = Pchol_MC*[1, 0, 0, 0, 0]
    println(shock_MC[1])
    responses_MC   = IRF(A_MC, 0*μolsct[:,1], tIRF-1, shock_MC)' 

    ## Fill in G
    for ss = 1:tIRF
        ψ          = vec(responses[ss,:])
        ψ_e        = vec(responses_MC[ss,:])
        G[:,ee,ss] = (ψ_e - ψ)./Δ
    end
end

## Calculate variance-covariance matrix
VCOV      = zeros(N,N,tIRF)
VCOVnosym = zeros(N,N,tIRF)
SD        = zeros(tIRF,N)
for ss = 1:tIRF
    VCOVnosym[:,:,ss] = G[:,:,ss] * kron(Σolsct,inv(Q))  * (G[:,:,ss].')/(TT-P)
    VCOV[:,:,ss]      = convert(Array{Float64,2},Symmetric(VCOVnosym[:,:,ss]))
    SD[ss,:]          = sqrt.(diag(VCOV[:,:,ss]))
end


## Confidence intervals
CIlow = zeros(tIRF,N)
CIhigh = zeros(tIRF,N)    
for nn = 1:N
    CIlow[:,nn]  = responses[:,nn]  -1.28*SD[:,nn]
    CIhigh[:,nn] = responses[:,nn] +1.28*SD[:,nn]
end


responses = responses * 100
CIlow = CIlow * 100
CIhigh = CIhigh * 100


plotIRF(responses[:,1], "Commodity Price (Deviation from mean)", "output/IRF_P_P";
        makeLegend=false  ,CIlow=CIlow[:,1], CIhigh=CIhigh[:,1],Clevel=80)
plotIRF(responses[:,2], "Trade Balance/GDP"                    , "output/IRF_TBY_P";
        makeLegend=false,CIlow=CIlow[:,2], CIhigh=CIhigh[:,2],Clevel=80,YTICKS=collect(-1.5:0.5:0.5))
plotIRF(responses[:,3], "Log Real Investment"                  , "output/IRF_I_P";
        makeLegend=false  ,CIlow=CIlow[:,3], CIhigh=CIhigh[:,3],Clevel=80,YTICKS=collect(-5:5:10))
plotIRF(responses[:,4], "Log Real Consumption"                 , "output/IRF_c_P";
        makeLegend=false  ,CIlow=CIlow[:,4], CIhigh=CIhigh[:,4],Clevel=80,YTICKS=collect(-1:1:3))
plotIRF(responses[:,5], "Log Real GDP  Capita"                 , "output/IRF_Y_P";
        makeLegend=false  ,CIlow=CIlow[:,5], CIhigh=CIhigh[:,5],Clevel=80,YTICKS=collect(-1:1:3))
