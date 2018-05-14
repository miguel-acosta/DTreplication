using PyPlot, Stats
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
function plotIRF(y, idx, shock, vname, fname, ;units="Percent", CI = false,Clevel=90,LINELABEL = [""], TITLE="",IRFCI=[])
    T = size(IRFCI)[2]    
    f = figure(figsize = (5,4))
    plot(1:T, repmat([0],T,1), color = "black", linewidth = 1)

    LINESTYLES = repmat(["-", "--", ":"], 1,10)
    COLORS = repmat(["#23a638", "#008bd2", "#f29fc5"], 1,10)

    
    if CI
        CIlow  = [percentile(vec(IRFCI[idx,tt,:]),Int((100-Clevel)/2)) for tt in 1:T]
        CIhigh = [percentile(vec(IRFCI[idx,tt,:]),Int(100-(100-Clevel)/2)) for tt in 1:T]
    end

    if length(size(y)) > 1
        for ss in 1:size(y)[2]
            plot(1:T, y[:,ss], linewidth = 2, label=LINELABEL[ss],
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
        fill_between(1:T, CIlow, CIhigh, facecolor="#d3f8a3",
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
