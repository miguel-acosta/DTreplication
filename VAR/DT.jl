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
Clevel = 80
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##--## ## ESTIMATE THE VAR
Aolsct,μolsct,Σolsct,residsct = VARols(P, datVAR; cons=true,  Z = trend)

Pchol = chol(Σolsct).'
shock = Pchol*[1, 0, 0, 0, 0]
responses = IRF(Aolsct, 0*μolsct[:,1], tIRF-1, shock)' * 100


##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## Confidence Intervals
Π = [μolsct reshape(Aolsct, size(Aolsct)[1], size(Aolsct)[2]*2)]
π = vec(Π.')
dat_lags = lagmatrix(datVAR,2)
Z   = dat_lags[:,N+1:end]
X  = [ones(TT-P,1) trend[(P+1):end] Z] ## the i-th row of this is x_t from hamilton
Q  = (X.'*X)/(TT-P)
VCOV = convert(Array{Float64,2},Symmetric(kron(Σolsct,inv(Q)))/(TT-P))
πdist = MvNormal(π, VCOV)
srand(8)
dall  = rand(πdist,Ndraws)
IRF_MC = zeros(tIRF, N, Ndraws)
for nn = 1:Ndraws
    d              = copy(dall[:,nn])
    d_unpack       = reshape(d, N*P+2,N)
    μ_MC           = d_unpack[1:2,:].'
    A_MC           = reshape(d_unpack[3:end,:].', N,N,2)
    Π_MC           = [μ_MC reshape(A_MC, size(A_MC)[1], size(A_MC)[2]*2)]
    Aols_MC,μolsMC,Σ_MC,resids_MC     = VARols(2, datVAR; cons=true,  Z = trend, βfixed = Π_MC)
    Pchol_MC       = chol(Σ_MC).'
    shock_MC       = Pchol_MC*[1, 0, 0, 0, 0]
    IRF_MC[:,:,nn] = IRF(A_MC, 0*μolsct[:,1], tIRF-1, shock_MC)' * 100
    if nn % 500    == 0
        println(string("iteration ", nn, " of ", Ndraws))
    end
        
end


CIlow = zeros(tIRF,N)
CIhigh = zeros(tIRF,N)
for nn = 1:N
    CIlow[:,nn] =  [percentile(vec(IRF_MC[tt,nn,:]),Int(    (100-Clevel)/2)) for tt in 1:tIRF]
    CIhigh[:,nn] = [percentile(vec(IRF_MC[tt,nn,:]),Int(100-(100-Clevel)/2)) for tt in 1:tIRF]
end
    
IRF_SD = std(IRF_MC,3)

crit = quantile(Normal(0,1),(100-Clevel)/200)
 for nn = 1:N
     CIlow[:,nn] = responses[:,nn]  -crit*IRF_SD[:,nn]
     CIhigh[:,nn] = responses[:,nn] +crit*IRF_SD[:,nn]
 end





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
