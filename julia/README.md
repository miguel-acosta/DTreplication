# Model Solution and Calibrated IRF

As noted in our writeup, we make use of Chris Sims' Gensys
routine to solve the model. This comes in a `Julia` package
provided by QuantEcon. The file `BBmodel.jl` contains our
model written in the form prescribed by Sims. 
The file `BBsteadystate.jl` computes the steady state of the
model. The file `calibration.jl` produces the impulse responses
at the calibration used by the authors, and makes use of
`IRFplot.jl`, which is just a plotting routine. 

# Model Estimation
As described in the main text, we estimate the model using
the random walk Metropolis-Hastings MCMC algorithm. The
implementation of this is in `estimation.jl`. This further
makes use of `kalmanLikelihood.jl` in order to evaluate the
likelihood of the model. The file `paramEstim.jl` prints
summary statistics of the estimated parameters.

# Estimated IRF and Variance Decomposition
The file `IRFestimated.jl` calculates the infinite
horizon forecast error variance decomposition of the
observable variables. It also computes the impulse
responses and associated 90% credible sets of the estimated models. 


# High Performance Computing
All files with the suffix `.sh` are the shell scripts
created to make use of Columbia's "Habanero" high performance
computing cluster. 
