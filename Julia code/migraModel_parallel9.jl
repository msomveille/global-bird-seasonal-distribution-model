using DataFrames
using HDF5
using JLD

@everywhere include("migraModel_forCluster.jl")

run1 = @spawn migraModel(10, 0, 0/3, 1+(4*(76.5/9)), 0.01+(1/3*(0.24/9)))
run2 = @spawn migraModel(10, 0, 1/3, 1+(4*(76.5/9)), 0.01+(1/3*(0.24/9)))
run3 = @spawn migraModel(10, 0, 2/3, 1+(4*(76.5/9)), 0.01+(1/3*(0.24/9)))
run4 = @spawn migraModel(10, 0, 3/3, 1+(4*(76.5/9)), 0.01+(1/3*(0.24/9)))
run5 = @spawn migraModel(10, 0, 4/3, 1+(4*(76.5/9)), 0.01+(1/3*(0.24/9)))
#run6 = @spawn migraModel(10)
#run7 = @spawn migraModel(10)
#run8 = @spawn migraModel(10)
#run9 = @spawn migraModel(20)
#run10 = @spawn migraModel(20)
#run11 = @spawn migraModel(20)
#run12 = @spawn migraModel(20)
#run13 = @spawn migraModel(20)
#run14 = @spawn migraModel(20)
#run15 = @spawn migraModel(20)
#run16 = @spawn migraModel(20)

resmechh = hcat(fetch(run1),fetch(run2),fetch(run3),fetch(run4),fetch(run5)) #,fetch(run6),fetch(run7),fetch(run8)) #,fetch(run9),fetch(run10),fetch(run11),fetch(run12),fetch(run13),fetch(run14),fetch(run15),fetch(run16))

#abc = sample(collect(1:1:10000000),1)[1]
#pathh = string("resmec", abc, ".jld")
#save(pathh, "data", resmechh)
save("resmec9.jld", "data", resmechh)
