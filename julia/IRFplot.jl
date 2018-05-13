using PyPlot, Stats
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function plotIRF(IRF, idx, shock, vname, fname, ;units="Percent", CI = false,Clevel=90,LINELABEL = "", TITLE="",IRFCI=[])
    T = size(IRF)[2]    
    f = figure(figsize = (5,4))
    plot(1:T, repmat([0],T,1), color = "black", linewidth = 1)

    LINESTYLES = repmat(["-", "--", ":"], 1,10)
    COLORS = repmat(["black", "gray"], 1,10)

    y = IRF[idx, :].'
    if CI
        CIlow  = [percentile(vec(IRFCI[idx,tt,:]),Int((100-Clevel)/2)) for tt in 1:T]
        CIhigh = [percentile(vec(IRFCI[idx,tt,:]),Int(100-(100-Clevel)/2)) for tt in 1:T]
        println(CIlow)
        println(CIhigh)
    end

    if length(size(y)) > 1
        for ss in 1:size(y)[2]
            if (LINELABEL == "") & (typeof(vname) == String)
                LINELABEL = vname
            else
                LINELABEL = vname[ss]
            end
            
            plot(1:T, y[:,ss], linewidth = 2, label=LINELABEL,
                 linestyle=LINESTYLES[ss], color=COLORS[ss])            
        end
    else
        if LINELABEL == ""
            LINELABEL = vname
        end
        plot(1:T, y, linewidth = 2, label=LINELABEL)
        TITLE = vname
    end


    if length(IRFCI)>0
        fill_between(1:T, CIlow, CIhigh, facecolor="gray",
                     alpha=0.5,label=string(Clevel, "% CI"))
    end
    xlim(1,T)
    xlabel("t")
    title(TITLE)
    ylabel(units)
    legend()
    println(string(fname,".pdf"))
    savefig(string(fname,".pdf"))
    close(f)
    close()
end
