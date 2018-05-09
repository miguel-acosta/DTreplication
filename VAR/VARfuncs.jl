using IterableTables, DataFrames
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function VARols(P, variables; chop=0, cons = true, Z= [], βfixed=[])
    if chop > 0
        variables = variables[(chop+1):end,:]
    end
    V        = size(variables)[2]
    TT       = size(variables)[1]
    A        = zeros(V,V,P)
    ## Vector to store coeff. on constant or any other non-dynamic variables Z
    if cons == true && !isempty(Z)
        Z = [ones(TT-P,1) Z[(P+1):TT,]]
        nz = size(Z)[2]
    elseif cons == true && isempty(Z)
        Z = ones(TT-P,1)
        nz = 1
    elseif cons == false && !isempty(Z)
        Z  = Z[(P+1):TT,]
        nz = length(size(Z)) > 1 ? size(Z)[2] : 1
    else
        nz = 0
    end
    μ        = zeros(V, nz)
    
    u        = zeros(TT-P,V)

    for depvar = 1:V
        ## Get y variable
        y = variables[(P+1):TT,depvar]

        ## Create x matrix
        X = zeros(TT-P, P*V+nz)
        if nz > 0
            X[:,1:nz] = Z
        end
        col = nz+1
        for pp = 1:P
            for vv = 1:V
                x = variables[:,vv]
                X[:,col] = x[(P-pp+1):(TT-pp)]
                col += 1
            end
        end

        if !isempty(βfixed)
            β = copy(βfixed[depvar,:])
        else
            β = (X.' * X) \ X.' * y
        end
        μ[depvar,:] = β[1:nz]
        u[:,depvar] = y - X*β

        for pp = 1:P
            inds = (V*(pp-1) + 1 + nz):(V*(pp) + nz)
            A[depvar,:,pp] = β[inds]
        end
    end
#    Σ = zeros(V,V)
#    for tt = 1:(TT-P)
#        Σ += u[tt,:] * u[tt,:].'
    #    end
#    print(mean(u))
#    print("\n")
    Σ = (u.'*u)/(size(u)[1])
    
#    Σ = cov(u[pmax:end,:],1,false)
    return(A, μ, Σ, u)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function moments(Y)
    LY = Y[1:(end-1),:]
    Y  = Y[2:end,:]
    TT = size(Y)[1]

    EY = repmat(mean(Y,1), TT,1)
    Γ0 = (Y-EY).' * (Y-EY) /TT
    Γ1 = (Y-EY).' * (LY-EY) /TT
    return EY[1,:], Γ0, Γ1
end


##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function companion(A,μ)
    N = length(μ)
    P = size(A)[3]
    NP = N*P
    ## Build constant
    M = repmat([0.0], NP,1)
    M[1:N] = μ
    ## Build companion
    C = repmat([0.0], NP,NP)
    C[1:N, :] = reshape(A, N, NP)
    C[N+1:NP,1:NP-N] = kron(eye(P-1), eye(N,N))
    ## Return
    return(M,C)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function VMA(A,μ,j)
    AA = companion(A,μ)[2]
    K = length(μ)
    P = size(A)[3]
    J = [eye(K) zeros(K,K*(P-1))]
    Φ = J * AA^j * J.'
    return(Φ)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function IRF(A, μ, nt, shock;G=eye(length(shock)))
    N = length(μ)
    P = size(A)[3]
    Y = zeros(N, nt+1)
    for jj = 0:nt
        Y[:,jj+1] = VMA(A,μ,jj) * G * shock
    end
    return(Y)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function plotts(y, name, snames;sub=true,
                units = repmat(["Percent"],length(snames),1),
                dims = (12,4), TITLE = "", LEGEND=false)
    N = size(y)[2]
    T = size(y)[1]    
    f = figure(figsize = dims)
    LINESTYLES = repmat(["-", "--", ":"], 1,3)
    COLORS = repmat(["black", "gray"], 1,3)    
    for ii = 1:length(snames)
        if sub
            subplot(100 + N*10 + ii)
        end
        plot(1:T, repmat([0],T,1), color = "black", linewidth = 1)
        plot(1:T, y[:,ii], linewidth = 2,label=snames[ii],color=COLORS[ii],
             linestyle=LINESTYLES[ii])
        xlim(1,T)
        xlabel("t")
        title(snames[ii])
        ylabel(units[ii])        
    end
    if length(TITLE) > 0
        title(TITLE)
    end
    if LEGEND
        legend()
    end
    
    savefig(string(name,".pdf"))
    close(f)
    close()
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function lagmatrix(y, P)
    if length(size(y)) > 1
        TT,K = size(y)
    else
        TT = size(y)[1]
        K = 1
    end
    L = zeros(TT-P, K * (P+1))
    for ll in 1:(P+1)
        L[:,(K*(ll-1)+1):(ll*K)] = Array(y[P-ll+2:(TT-ll+1),:])
    end
    return(L)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function simulateVAR(Y,k,p,companion, mean, resids)
    Y = Array(Y)
    LY = lagmatrix(Y,p)
    LAGS = LY[:,k+1:end]
    TT = size(LAGS)[1]
    Ysim = zeros(TT,k)
    for tt = 1:TT
        fcast = mean + companion * LAGS[tt,:] + [resids[tt,:]; zeros(k*(p-1))]
        Ysim[tt,:] = fcast[1:k]
        LY = lagmatrix([Y[1:p,:] ;Ysim],p)
        LAGS = LY[:,k+1:end]
    end
    return(Ysim)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function forecast(h,A,μ,errors)
    K    = length(μ)
    TT   = size(errors)[1]
    yth  = zeros(K,1)
    Φ(j) = VMA(A,μ,j)
    for jj = 1:TT
        yth += Φ(h+jj-1) * errors[TT-jj+1,:]
    end
    return(yth+μ)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
## https://www3.nd.edu/~esims1/time_series_notes.pdf
function FEV(A,μ,h)
    K        = size(μ)[1]
    fev      = zeros(K,1)
    fevShare = zeros(K,K)
    for kk1 = 1:K
        for hh = 0:h
            fev[kk1] += sum(VMA(A,μ,hh)[kk1,:].^2)
        end
        for kk2 = 1:K
            for hh = 0:h
                fevShare[kk1,kk2] += VMA(A,μ,hh)[kk1,kk2]^2
            end
            fevShare[kk1,kk2] = fevShare[kk1,kk2]/fev[kk1]
        end
    end
    return(fev, fevShare)
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function plotFEV(vname, fev)
    f=figure(figsize = (4,4))
    stacks = stackplot(1:12, fev.', colors = ["mediumspringgreen", "midnightblue", "m"])
    hatches=["/", "*", "o"]
    for ii in 1:length(stacks)
        stacks[ii][:set_hatch](hatches[ii])
    end
    legend(["GDP", "M3", "FF"],loc="lower right")
    PyPlot.title(string("Forecast Error Variance Decomposition: ", vname))
    savefig(string("FEV", vname, ".pdf"))
    close(f)
    close()
end

##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function runkle(Y, A, μ, resids, k, p, nbs, tirf;recursiveID=false,Clevel=95)
    clow = (100-Clevel)/200
    chigh = 1-clow
    TT    = size(resids)[1]
    M,C   = companion(A,μ)
    IRFbs = zeros(tirf,k,nbs,k)
    for ss = 1:nbs
        ## Draw sample
        bssample = wsample(1:TT, repmat([1/TT],TT), TT)
        residsBS = resids[bssample,:]
        Ybs      = simulateVAR(Y,k,p,C,M,residsBS)
        ## Re-estimate A
        Abs,μbs,Σbs,rbs  = VARols(p,DataFrame(Ybs))
        ## Structural?
        if recursiveID
            P=chol(Σbs).'
        else
            P=eye(length(μbs))
        end
        ## Compute IRF
        for kk = 1:k
            shock = zeros(k,1)
            shock[kk] = 1
            IRFbs[:,:,ss,kk] = IRF(Abs,μbs,tirf,P*shock).'
        end
    end
    IRFbs025 = zeros(tirf,k,k)
    IRFbs975 = zeros(tirf,k,k)    
    for kk1 = 1:k
        for kk2 = 1:k
            IRFbs025[:,kk1,kk2] = [quantile(IRFbs[tt,kk1,:,kk2], clow) for tt in 1:tirf]
            IRFbs975[:,kk1,kk2] = [quantile(IRFbs[tt,kk1,:,kk2], chigh) for tt in 1:tirf]            
        end
    end
    return(IRFbs025, IRFbs975)
end


##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function plotIRF(y, vname, fname, ;units="Percent", CIlow = [], CIhigh=[],Clevel=95,TITLE="", makeLegend=true, YTICKS = [])
    T = size(y)[1]    
    f = figure(figsize = (5,4))
    plot(1:T, repmat([0],T,1), color = "black", linewidth = 1)

    LINESTYLES = repmat(["-", "--", ":"], 1,10)
    COLORS = repmat(["black", "gray"], 1,10)    
    
    if length(size(y)) > 1
        for ss in 1:size(y)[2]
            plot(1:T, y[:,ss], linewidth = 2, label=vname[ss],
                 linestyle=LINESTYLES[ss], color=COLORS[ss]; myfonts)
        end
    else
        plot(1:T, y, linewidth = 2, label=vname,color="#ff3341")
        TITLE = vname
    end


    if length(CIlow)>0
        fill_between(1:T, CIlow, CIhigh, facecolor="#33eeff",
                     alpha=0.5,label=string(Clevel, "% CI"))
    end
    xlim(1,T)
    if length(YTICKS)>0
        yticks(YTICKS)
    end
    xlabel("t")
    title(TITLE)
    ylabel(units)
    if makeLegend
        legend()
    end
    savefig(string(fname,".pdf"))
    close(f)
    close()
end
