include("VARfuncs.jl")
using IterableTables, DataFrames, ExcelReaders, PyPlot #, HypothesisTests

dat = readxlsheet(DataFrame, "DataVAR.xlsx", "Sheet1", header=true)
datVAR = dat[:,2:end]
datVAR = datVAR[:,end:-1:1]
TT  = length(dat[1])
P   = 2
N   = size(datVAR)[2]


function IRFdt(A, α, shock; TT = 20)
    responses = zeros(TT, size(A)[2])
    responses[1,:] = shock
    responses[2,:] = A[:,:,1] * responses[1,:] + α
    for tt = 3:TT
        responses[tt,:] = A[:,:,1] * responses[tt-1,:] + A[:,:,2] * responses[tt-2,:] + α * (tt-1)
    end
    return(responses)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##--## ## ESTIMATE THE VAR
Aols,μols,Σols,resids         = VARols(2, datVAR; cons=false)
Aolsc,μolsc,Σolsc,residsc     = VARols(2, datVAR)
Aolst,μolst,Σolst,residst     = VARols(2, datVAR; cons=false, Z = collect(1:TT))
Aolsct,μolsct,Σolsct,residsct = VARols(2, datVAR; cons=true,  Z = collect(1:TT))

Pchol = chol(Σolsct)'
shock = Pchol*[1, 0, 0, 0, 0]
responses = IRF(Aolsct, 0*μolsct[:,1], 10, shock)' * 100
responsess = IRFdt(Aolsct, 0, shock;TT=10) * 100
plotIRF(responsess[:,1], "Commodity Price (Deviation from mean)", "output/IRF_P_P")
plotIRF(responsess[:,2], "Trade Balance/GDP"                    , "output/IRF_TBY_P")
plotIRF(responsess[:,3], "Log Real Investment"                  , "output/IRF_I_P")
plotIRF(responsess[:,4], "Log Real Consumption"                 , "output/IRF_c_P")
plotIRF(responsess[:,5], "Log Real GDP  Capita"                 , "output/IRF_Y_P")