using NLsolve
using Parameters
function BBsteadystate(PAR::Dict{String,Float64})

    @unpack p_til, dstar, s, ξ, g, Ψ, αk, αm, αk_til,
            δ, ϕ,b, Γ, θ, ω, ω_til = PAR
    
    # Solve for the steady state
    r          = (1/b)*(g^(Γ))-1

    d          = dstar
    rstar      = r
    r1         = (1/b)*(g^Γ)-(1-δ)
    r2         = r1
  
    y_til_k2   = r2/(αk_til*p_til)
    k2_n2      = g*(y_til_k2)^(-1/(1-αk_til))
    y_til_n2   = g^(1-αk_til)*(k2_n2)^αk_til
    w2         = (1-αk_til)*(y_til_n2)*p_til
    n2         = (w2/θ)^(1/(ω_til-1))
    k2         = (k2_n2)*n2
    y_til      = (y_til_n2)*n2
  
    y_k1       = r1/αk
    y_m_til    = p_til/αm

    function f!(F,x)
        if x[1] <0 || x[2] < 0
            F[1] = 100000
            F[2] = 100000
        else
            F[1] = g^(1-αk-αm)*(x[1])^(αk-1)*(x[2])^(αm) - y_k1
            F[2] = g^(1-αk-αm)*(x[1])^(αk)*(x[2])^(αm-1) - y_m_til
        end
    end

#    opt        = Opt(:LD_MMA,2)
#    lower_bounds!(opt, [0.,0.])
#    min_objective!(opt,eqns)

#    (minf,minx,ret) = optimize(opt, [1.0,1.0])
    
    #    S          = optimize(eqns, [1.0,1.0], BFGS())
    S          = nlsolve(f!, [1.0,1.0])
    k1_n1      = S.zero[1]  #minx[1] 1.33570465442287390000 
    m_til_n1   = S.zero[2]  #minx[2] 0.09360122597253398000 
  
    y_n1       = g^(1-αk-αm)*(k1_n1)^(αk)*(m_til_n1)^(αm)
    w1         = (1-αk-αm)*(y_n1)
    n1         = (w1/θ)^(1/(ω-1))
    k1         = (k1_n1)*n1
    m_til      = (m_til_n1)*n1
    y          = (y_n1)*n1
               
    k          = k1+k2
    INV        = k*g-(1-δ)*k
    c          = y+p_til*y_til-p_til*m_til-INV+d*(g/(1+r)-1)-s
    λ     = (c-(θ/ω)*n1^ω-(θ/ω_til)*n2^ω_til)^(-Γ)
    tb         = y-c-k*g+(1-δ)*k-s
    tb_til     = p_til*(y_til-m_til)
    ygdp       = y+tb_til
    ygdpf      = y - p_til*m_til
    tbtotal    = tb+tb_til

  
    # Pack up the steady state
    ss = Dict("d"        => d        ,
              "r"        => r        ,
              "rstar"    => rstar    ,
              "r1"       => r1       ,
              "r2"       => r2       ,
              "y_til_k2" => y_til_k2 ,
              "k2_n2"    => k2_n2    ,
              "y_til_n2" => y_til_n2 ,
              "w2"       => w2       ,
              "n2"       => n2       ,
              "k2"       => k2       ,
              "y_til"    => y_til    ,
              "y_k1"     => y_k1     ,
              "y_m_til"  => y_m_til  ,
              "k1_n1"    => k1_n1    ,
              "m_til_n1" => m_til_n1 ,
              "y_n1"     => y_n1     ,
              "w1"       => w1       ,
              "n1"       => n1       ,
              "k1"       => k1       ,
              "m_til"    => m_til    ,
              "y"        => y        ,
              "k"        => k        ,
              "INV"      => INV      ,
              "c"        => c        ,
              "λ"        => λ        ,
              "tb"       => tb       ,
              "tb_til"   => tb_til   ,
              "ygdp"     => ygdp     ,
              "ygdpf"    => ygdpf    ,
              "tbtotal"  => tbtotal  )
    
    return(merge(PAR,ss))
end


