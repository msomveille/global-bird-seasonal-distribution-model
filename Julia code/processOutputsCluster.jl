using DataFrames
using HDF5
using JLD
using Gadfly
using Distances
using StatsBase

Gadfly.push_theme(:dark)

cd("$(homedir())/Desktop/PhD/Chapter 4 – mechanistic model/OutputsCluster/outputs_e0.0188")
pwd()
readdir(pwd())

## load and process results

resmec_BR = zeros(7352)
resmec_NB = zeros(7352)
resmec_res = zeros(7352)
resmec_mm = zeros(1)
resmec_rr = zeros(1)
resmec_tt = zeros(1)
resmec_ee = zeros(1)
for j in 1:length(readdir(pwd()))
  for i in 1:5
    resmec_BR = hcat(resmec_BR, load(readdir(pwd())[j])["data"][i][1][:,2:end][:,collect(1:3:(10*3))])
    resmec_NB = hcat(resmec_NB, load(readdir(pwd())[j])["data"][i][1][:,2:end][:,collect(2:3:(10*3))])
    resmec_res = hcat(resmec_res, load(readdir(pwd())[j])["data"][i][1][:,2:end][:,collect(3:3:(10*3))])
    resmec_mm = vcat(resmec_mm, load(readdir(pwd())[j])["data"][i][2][1,2:end])
    resmec_rr = vcat(resmec_rr, load(readdir(pwd())[j])["data"][i][2][2,2:end])
    resmec_tt = vcat(resmec_tt, load(readdir(pwd())[j])["data"][i][2][3,2:end])
    resmec_ee = vcat(resmec_ee, load(readdir(pwd())[j])["data"][i][2][4,2:end])
  end
end
resmec_BR = resmec_BR[:,2:end]
resmec_NB = resmec_NB[:,2:end]
resmec_res = resmec_res[:,2:end]
resmec_mm = resmec_mm[2:end]
resmec_rr = resmec_rr[2:end]
resmec_tt = resmec_tt[2:end]
resmec_ee = resmec_ee[2:end]
resmec_param = hcat(resmec_mm, resmec_rr, resmec_tt, resmec_ee)


writedlm("resmec_BR2.csv", resmec_BR, '\t')
writedlm("resmec_NB2.csv", resmec_NB, '\t')
writedlm("resmec_residents2.csv", resmec_res, '\t')
writedlm("resmec_param2.csv", resmec_param, '\t')


k = 2031
plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=resmec_BR[:,k], style(default_point_size=1.3pt))
plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=resmec_NB[:,k], style(default_point_size=1.3pt))
plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=(resmec_BR[:,k].+resmec_NB[:,k] .+ 1)./(resmec_BR[:,k].+resmec_NB[:,k].+resmec_res[:,k] .+ 1), style(default_point_size=1.3pt))

plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=resmec_BR[:,k], style(default_point_size=1.3pt))
plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=resmec_NB[:,k], style(default_point_size=1.3pt))
plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=(resmec_BR[:,k].-resmec_NB[:,k]), style(default_point_size=1.3pt))
plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=(resmec_BR[:,k].+resmec_NB[:,k])./(resmec_BR[:,k].+resmec_NB[:,k].+resmec_res[:,k]), style(default_point_size=1.3pt))

plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=vcat(breeding_WH,breeding_EH), style(default_point_size=1.3pt))
plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=vcat(nonbreeding_WH,nonbreeding_EH), style(default_point_size=1.3pt))
plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=vcat(prop_migr_obs_WH,prop_migr_obs_EH).+0.0001, style(default_point_size=1.3pt))


cor(resmec_BR[:,k], vcat(breeding_WH,breeding_EH))
cor(resmec_NB[:,k], vcat(nonbreeding_WH,nonbreeding_EH))
cor((resmec_BR[:,80].+resmec_NB[:,80])./(resmec_BR[:,80].+resmec_NB[:,80].+resmec_res[:,80]), vcat(prop_migr_obs_WH,prop_migr_obs_EH))
