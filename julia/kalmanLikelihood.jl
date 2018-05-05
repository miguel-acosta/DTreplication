using Distributions


function kalmanLikelihood(parameters::Dict{String,Float64}, data::Array{Float64,2})
    # inputs: parameters, data 
    ss = BBsteadystate(parameters)
    G1, C0, G0, fmat, fwt, ywt, gev, eu, ind, NY, NEPS, M1 = BBmodel(ss)

    if sum(eu) < 2
        return(-Inf)
    end
    
    # Measurement Eqn: y(t) = M0+M1*s(t)+ U(t). U(t)~N(0, VarU)
    # NOTE: M0 = zeros(NY, 1), U = zeros(NY, 1)
    # Transition Eqn: s(t) = G1*s(t-1)+G0*EPS(t). EPS(t)~N(0, VarEps)

    NS = size(M1)[2]

    # Assign the initial expectation and var-covariance matrix.
    MeanS = zeros(NS,1);
    VarS  = zeros(NS,NS) # Sigma^s (1|0)=unconditional variance
    VarEps = eye(NEPS)
    error = 1; count = 0; maxit = 1e5; tol = 1e-5;
    while ((error>tol) & (count<maxit))
        VarS1 = G1*VarS*(G1.')+ G0*VarEps*(G0.');
        error = maximum(maximum(abs.(VarS1-VarS)));
        count += 1;
        VarS  = copy(VarS1);
    end
    
    # Derive the likelihood and combine it with the prior. 
    multi(MD,VAR) = -((NY/2)*log(2*pi)+log(det(VAR))/2+ (MD.'*inv(VAR)*MD)/2);
    # multi is the logged conditional multi-variate normal pdf
    # with input: MD(conditional mean deviation) and VAR(conditional variance)
    
    # Apply Kalman Filter and get the Likelihood LF.
    LF = 0;
    for t = 1:size(data)[2]
        Y1 = M1*MeanS # % mean, (NYx1)
        Y2 = M1*VarS*(M1.') # % variance, (NYxNY)
        if det(Y2)<=0
            return(-Inf)
        end
        LF = LF+multi(data[:,t]-Y1, Y2)
        S1 = MeanS + VarS*(M1.')*inv(Y2)*(data[:,t]-Y1) #% mean, (NSx1)
        S2 = VarS - VarS*(M1.')*inv(Y2)*M1*VarS # variance, (NSxNS)
        MeanS = G1*S1 
        VarS = G1*S2*(G1.') + G0*VarEps*(G0.')
    end
    
    return(LF[1])
end

