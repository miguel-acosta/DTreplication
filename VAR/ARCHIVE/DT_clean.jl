include("VARfuncs.jl")
using IterableTables, DataFrames, ExcelReaders, PyPlot #, HypothesisTests

dat = convert(Array{Float64,2},readxlsheet("DataVAR.xlsx", "Sheet1", skipstartrows=1))
#dat = readxlsheet(DataFrame, "DataVAR.xlsx", "Sheet1", header=true)
datVAR = dat[:,2:end]
datVAR = datVAR[:,end:-1:1]
TT  = size(dat)[1]
P   = 2
N   = size(datVAR)[2]


#function IRFdt(A, α, shock; TT = 20)
#    responses = zeros(TT, size(A)[2])
#    responses[1,:] = shock
#    responses[2,:] = A[:,:,1] * responses[1,:] + α
#    for tt = 3:TT
#        responses[tt,:] = A[:,:,1] * responses[tt-1,:] + A[:,:,2] * responses[tt-2,:] + α * (tt-1)
#    end
#    return(responses)
#end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##--## ## ESTIMATE THE VAR
#Aols,μols,Σols,resids         = VARols(2, datVAR; cons=false)
#Aolsc,μolsc,Σolsc,residsc     = VARols(2, datVAR)
#Aolst,μolst,Σolst,residst     = VARols(2, datVAR; cons=false, Z = collect(1:TT))
Aolsct,μolsct,Σolsct,residsct = VARols(2, datVAR; cons=true,  Z = collect(1:TT))

Pchol = chol(Σolsct)'
shock = Pchol*[1, 0, 0, 0, 0]
responses = IRF(Aolsct, 0*μolsct[:,1], 9, shock)' * 100
#responses = IRFdt(Aolsct, 0, shock;TT=10) * 100
plotIRF(responses[:,1], "Commodity Price (Deviation from mean)", "output/IRF_P_P"; makeLegend=false)
plotIRF(responses[:,2], "Trade Balance/GDP"                    , "output/IRF_TBY_P"; makeLegend=false)
plotIRF(responses[:,3], "Log Real Investment"                  , "output/IRF_I_P"; makeLegend=false)
plotIRF(responses[:,4], "Log Real Consumption"                 , "output/IRF_c_P"; makeLegend=false)
plotIRF(responses[:,5], "Log Real GDP  Capita"                 , "output/IRF_Y_P"; makeLegend=false)
