set more off

use VARdata, clear

tsset Year
sort Year
gen t = _n
var com gdp con inv tb , lags(1/2)  exog(t) nocons // closest i can get to them
// var com gdp con inv tb , lags(1/2)  exog(t)  // what we do 
irf create model1, step(9) set(irf1, replace)
irf graph oirf, impulse(com) response(gdp con inv tb com) byopts(rescale rows(1)) level(80)


graph export irf.pdf, replace
