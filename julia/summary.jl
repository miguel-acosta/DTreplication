using CSV


param =  c("ξ", "ψ", "ρ_p1", "ρ_p2", "σ_p",
           "ρ_a", "ρ_a_til", "ρ_g", "ρ_s", "ρ_nu", "ρ_mu",
           "σ_a", "σ_a_til", "σ_g", "σ_s", "σ_nu", "σ_mu")

##postFiles <- readdir("posterior")
postFiles <- ["post9.csv" , "post8.csv"]
draws <- zeros(length(param),1)
for postFile in postFiles 
    tempFile <- CSV.read(string("posterior/", postFile))
    draws = [draws postFile]
end

