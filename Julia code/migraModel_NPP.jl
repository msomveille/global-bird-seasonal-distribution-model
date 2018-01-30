using DataFrames
using Gadfly
using Distances
using StatsBase
using GLM

Gadfly.push_theme(:dark)

cd("$(homedir())/Desktop/PhD/Chapter 4 – mechanistic model/Loading-data")
pwd()

PresAbs_BR_NH = Array(readtable("PresAbs_BR_NH.csv")) ## dim = 7613, 1185
PresAbs_NB_NH = Array(readtable("PresAbs_NB_NH.csv"))
PresAbs_BR_SH = Array(readtable("PresAbs_BR_SH.csv")) ## dim = 7613, 256
PresAbs_NB_SH = Array(readtable("PresAbs_NB_SH.csv"))
PresAbs_res = Array(readtable("PresAbs_res.csv"))

# Get climate data
climateDataset = readtable("env_data_6months_npp.csv")

hexagonstokeep = (!isna(climateDataset[:,1])) | (!isna(climateDataset[:,2])) | (!isna(climateDataset[:,3])) | (!isna(climateDataset[:,4])) | (!isna(climateDataset[:,5])) | (!isna(climateDataset[:,6])) | (!isna(climateDataset[:,7]))
hexagonstokeepWH = hexagonstokeep & (climateDataset[:Longitude] .<= -30)
hexagonstokeepEH = hexagonstokeep & (climateDataset[:Longitude] .> -30)

climateData = climateDataset[hexagonstokeep,:]

## Subset the global hexagon grid (with species' presence/absence) into Western and Eastern Hemispheres
PresAbs_BR_NH_WH = PresAbs_BR_NH[hexagonstokeepWH, 2:end]
PresAbs_NB_NH_WH = PresAbs_NB_NH[hexagonstokeepWH, 2:end]
PresAbs_BR_SH_WH = PresAbs_BR_SH[hexagonstokeepWH, 2:end]
PresAbs_NB_SH_WH = PresAbs_NB_SH[hexagonstokeepWH, 2:end]
PresAbs_res_WH = PresAbs_res[hexagonstokeepWH, 2:end]
PresAbs_BR_NH_EH = PresAbs_BR_NH[hexagonstokeepEH, 2:end]
PresAbs_NB_NH_EH = PresAbs_NB_NH[hexagonstokeepEH, 2:end]
PresAbs_BR_SH_EH = PresAbs_BR_SH[hexagonstokeepEH, 2:end]
PresAbs_NB_SH_EH = PresAbs_NB_SH[hexagonstokeepEH, 2:end]
PresAbs_res_EH = PresAbs_res[hexagonstokeepEH, 2:end]

# Western-Northern Hemisphere
spprangeBR = zeros(1184)
spprangeNB = zeros(1184)
for i in 1:1184
  spprangeBR[i] = sum(PresAbs_BR_NH_WH[:,i])
  spprangeNB[i] = sum(PresAbs_NB_NH_WH[:,i])
end
speciestokeep = (spprangeBR .> 0) & (spprangeNB .> 0)
PresAbs_BR_NH_WH = PresAbs_BR_NH_WH[:,speciestokeep]
PresAbs_NB_NH_WH = PresAbs_NB_NH_WH[:,speciestokeep]
# Western-Southern Hemisphere
spprangeBR = zeros(255)
spprangeNB = zeros(255)
for i in 1:255
  spprangeBR[i] = sum(PresAbs_BR_SH_WH[:,i])
  spprangeNB[i] = sum(PresAbs_NB_SH_WH[:,i])
end
speciestokeep = (spprangeBR .> 0) & (spprangeNB .> 0)
PresAbs_BR_SH_WH = PresAbs_BR_SH_WH[:,speciestokeep]
PresAbs_NB_SH_WH = PresAbs_NB_SH_WH[:,speciestokeep]
# Eastern-Northern Hemisphere
spprangeBR = zeros(1184)
spprangeNB = zeros(1184)
for i in 1:1184
  spprangeBR[i] = sum(PresAbs_BR_NH_EH[:,i])
  spprangeNB[i] = sum(PresAbs_NB_NH_EH[:,i])
end
speciestokeep = (spprangeBR .> 0) & (spprangeNB .> 0)
PresAbs_BR_NH_EH = PresAbs_BR_NH_EH[:,speciestokeep]
PresAbs_NB_NH_EH = PresAbs_NB_NH_EH[:,speciestokeep]
# Eastern-Southern Hemisphere
spprangeBR = zeros(255)
spprangeNB = zeros(255)
for i in 1:255
  spprangeBR[i] = sum(PresAbs_BR_SH_EH[:,i])
  spprangeNB[i] = sum(PresAbs_NB_SH_EH[:,i])
end
speciestokeep = (spprangeBR .> 0) & (spprangeNB .> 0)
PresAbs_BR_SH_EH = PresAbs_BR_SH_EH[:,speciestokeep]
PresAbs_NB_SH_EH = PresAbs_NB_SH_EH[:,speciestokeep]

# Residents Western Hemisphere
spprangeRes = zeros(7819)
for i in 1:7819
  spprangeRes[i] = sum(PresAbs_res_WH[:,i])
end
speciestokeep = (spprangeRes .> 0)
PresAbs_res_WH = PresAbs_res_WH[:,speciestokeep]
# Residents Eastern Hemisphere
spprangeRes = zeros(7819)
for i in 1:7819
  spprangeRes[i] = sum(PresAbs_res_EH[:,i])
end
speciestokeep = (spprangeRes .> 0)
PresAbs_res_EH = PresAbs_res_EH[:,speciestokeep]


# From climate dataset, get climate and NDVI data for Northern Summer (NS; April to September) and Northern Winter (NW; November to March)
precNS = log(climateData[:precPresentNS]+1)
precNW = log(climateData[:precPresentNW]+1)
tempNS = climateData[:tempPresentNS]
tempNW = climateData[:tempPresentNW]
nppNS = climateData[:nppPresentNS]
nppNW = climateData[:nppPresentNW]
nppNS[nppNS .< 0] = 0
nppNS = log(nppNS.+1)
nppNW[nppNW .< 0] = 0
nppNW = log(nppNW.+1)
habitat = climateData[:Biome]

lonlat = climateData[:,[9,8]]
lonlatWH = lonlat[climateData[:Longitude] .< -30, :]
lonlatEH = lonlat[climateData[:Longitude] .>= -30, :]

#plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=vcat(energySupplyAnnualWH, energySupplyAnnualEH), style(default_point_size=1.3pt))
#plot(x=lonlat[:,1], y=lonlat[:,2], color=nppNS - nppNW, style(default_point_size=1.3pt))
#plot(x=lonlat[:,2], y=nppNW, style(default_point_size=1.3pt))


# z-tranform climate data to make axes comparable
sdTemp = std([tempNS; tempNW])
mnTemp = mean([tempNS; tempNW])
sdPrec = std([precNS; precNW])
mnPrec = mean([precNS; precNW])
tempNS_std = (tempNS - mnTemp) / sdTemp
tempNW_std = (tempNW - mnTemp) / sdTemp
precNS_std = (precNS - mnPrec) / sdPrec
precNW_std = (precNW - mnPrec) / sdPrec

# Subset climate data (and z-tranformed climate data) and NDVI data for western versus eastern hemispheres
tempNS_std_WH = tempNS_std[climateData[:Longitude] .< -30, :]
tempNW_std_WH = tempNW_std[climateData[:Longitude] .< -30, :]
precNS_std_WH = precNS_std[climateData[:Longitude] .< -30, :]
precNW_std_WH = precNW_std[climateData[:Longitude] .< -30, :]
tempNS_WH = tempNS[climateData[:Longitude] .< -30, :]
tempNW_WH = tempNW[climateData[:Longitude] .< -30, :]
precNS_WH = precNS[climateData[:Longitude] .< -30, :]
precNW_WH = precNW[climateData[:Longitude] .< -30, :]
nppNS_WH = nppNS[climateData[:Longitude] .< -30, :]
nppNW_WH = nppNW[climateData[:Longitude] .< -30, :]
habitat_WH = habitat[climateData[:Longitude] .< -30, :]
tempNS_std_EH = tempNS_std[climateData[:Longitude] .>= -30, :]
tempNW_std_EH = tempNW_std[climateData[:Longitude] .>= -30, :]
precNS_std_EH = precNS_std[climateData[:Longitude] .>= -30, :]
precNW_std_EH = precNW_std[climateData[:Longitude] .>= -30, :]
tempNS_EH = tempNS[climateData[:Longitude] .>= -30, :]
tempNW_EH = tempNW[climateData[:Longitude] .>= -30, :]
precNS_EH = precNS[climateData[:Longitude] .>= -30, :]
precNW_EH = precNW[climateData[:Longitude] .>= -30, :]
nppNS_EH = nppNS[climateData[:Longitude] .>= -30, :]
nppNW_EH = nppNW[climateData[:Longitude] .>= -30, :]
habitat_EH = habitat[climateData[:Longitude] .>= -30, :]

temp_meanAnnual_std_WH = (tempNS_std_WH .+ tempNW_std_WH) ./ 2
prec_meanAnnual_std_WH = (precNS_std_WH .+ precNW_std_WH) ./ 2
temp_meanAnnual_std_EH = (tempNS_std_EH .+ tempNW_std_EH) ./ 2
prec_meanAnnual_std_EH = (precNS_std_EH .+ precNW_std_EH) ./ 2



PresAbs_NS_WH = hcat(PresAbs_BR_NH_WH, PresAbs_NB_SH_WH, PresAbs_res_WH)
PresAbs_NW_WH = hcat(PresAbs_NB_NH_WH, PresAbs_BR_SH_WH, PresAbs_res_WH)
PresAbs_NS_EH = hcat(PresAbs_BR_NH_EH, PresAbs_NB_SH_EH, PresAbs_res_EH)
PresAbs_NW_EH = hcat(PresAbs_NB_NH_EH, PresAbs_BR_SH_EH, PresAbs_res_EH)

RichnessNS_WH = zeros(length(PresAbs_BR_NH_WH[:,1]))
RichnessNW_WH = zeros(length(PresAbs_BR_NH_WH[:,1]))
for i in 1:length(PresAbs_NS_WH[:,1])
  RichnessNS_WH[i] = sum(PresAbs_NS_WH[i,:])
  RichnessNW_WH[i] = sum(PresAbs_NW_WH[i,:])
end
RichnessNS_EH = zeros(length(PresAbs_BR_NH_EH[:,1]))
RichnessNW_EH = zeros(length(PresAbs_BR_NH_EH[:,1]))
for i in 1:length(PresAbs_NS_EH[:,1])
  RichnessNS_EH[i] = sum(PresAbs_NS_EH[i,:])
  RichnessNW_EH[i] = sum(PresAbs_NW_EH[i,:])
end



plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=nppNS_WH[:,1] - nppNW_WH[:,1], style(default_point_size=1.5pt))
plot(x=lonlatWH[:,1][(lonlatWH[:,2] .> 20)], y=lonlatWH[:,2][(lonlatWH[:,2] .> 20)], color=RichnessNW_WH[(lonlatWH[:,2] .> 20)], style(default_point_size=1.3pt))


aaa = DataFrame(NPP=nppNW_WH[:,1][lonlatWH[:,2] .> 0], Richness=sqrt(RichnessNW_WH)[lonlatWH[:,2] .> 0])
aaa2 = DataFrame(NPP=vcat(nppNW_WH[:,1][lonlatWH[:,2] .> 0], nppNW_EH[:,1][lonlatEH[:,2] .> 0]), Richness=vcat(RichnessNW_WH[lonlatWH[:,2] .> 0], RichnessNW_EH[lonlatEH[:,2] .> 0]))

plot(x=aaa2[:,1], y=aaa2[:,2], color=vcat(habitat_WH[(lonlatWH[:,2] .> 0)],habitat_EH[(lonlatEH[:,2] .> 0)]), style(default_point_size=1.5pt))

cor(aaa[:NPP], aaa[:Richness])

model1 = lm(Richness ~ 0 + NPP, aaa2)


ee = coef(model1)[1]

#energySupplyJulyWH = exp(ee * nppNS_WH[:,1]).-1
#energySupplyJanuaryWH = exp(ee * nppNW_WH[:,1]).-1
#energySupplyJulyEH = exp(ee * nppNS_EH[:,1]).-1
#energySupplyJanuaryEH = exp(ee * nppNW_EH[:,1]).-1

energySupplyJulyWH = ee * nppNS_WH[:,1]
energySupplyJanuaryWH = ee * nppNW_WH[:,1]
energySupplyJulyEH = ee * nppNS_EH[:,1]
energySupplyJanuaryEH = ee * nppNW_EH[:,1]


plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=energySupplyJanuaryWH, style(default_point_size=1.5pt))




function degtorad(deg)
  return(deg*pi/180)
end

function greatcircledistance(long1, lat1, long2, lat2)
  R = 6371 # Earth mean radius [km]
  deltaLong = long2 - long1
  deltaLat = lat2 - lat1
  a = sin(deltaLat/2)^2 + cos(lat1) * cos(lat2) * sin(deltaLong/2)^2
  c = 2 * asin(minimum([1 sqrt(a)]))
  d = R * c
  return(d) # Distance in km
end

pairwiseDistanceWH = zeros(length(lonlatWH[:,1]), length(lonlatWH[:,1]))
for i in 1:length(lonlatWH[:,1])
	for j in 1:length(lonlatWH[:,1])
		pairwiseDistanceWH[i,j] = greatcircledistance(degtorad(lonlatWH[:,1][i]), degtorad(lonlatWH[:,2][i]), degtorad(lonlatWH[:,1][j]), degtorad(lonlatWH[:,2][j]))
	end
end
pairwiseDistanceEH = zeros(length(lonlatEH[:,1]), length(lonlatEH[:,1]))
for i in 1:length(lonlatEH[:,1])
	for j in 1:length(lonlatEH[:,1])
		pairwiseDistanceEH[i,j] = greatcircledistance(degtorad(lonlatEH[:,1][i]), degtorad(lonlatEH[:,2][i]), degtorad(lonlatEH[:,1][j]), degtorad(lonlatEH[:,2][j]))
	end
end

#pairwiseDistanceWH = Array(readtable("pairwise_dist_WH_npp.csv", header = false))
pairwiseDistanceWH2 = copy(pairwiseDistanceWH)
pairwiseDistanceWH2[(pairwiseDistanceWH2 .< 250) & (pairwiseDistanceWH2 .> 0)] = 1
pairwiseDistanceWH2[(pairwiseDistanceWH2 .>= 250) | (pairwiseDistanceWH2 .<= 0)] = 0

#pairwiseDistanceEH = Array(readtable("pairwise_dist_EH_npp.csv", header = false))
pairwiseDistanceEH2 = copy(pairwiseDistanceEH)
pairwiseDistanceEH2[(pairwiseDistanceEH2 .< 250) & (pairwiseDistanceEH2 .> 0)] = 1
pairwiseDistanceEH2[(pairwiseDistanceEH2 .>= 250) | (pairwiseDistanceEH2 .<= 0)] = 0



### Observed patterns
migr_WH_NH = PresAbs_BR_NH_WH .- PresAbs_NB_NH_WH
migr_WH_SH = PresAbs_BR_SH_WH .- PresAbs_NB_SH_WH
migr_EH_NH = PresAbs_BR_NH_EH .- PresAbs_NB_NH_EH
migr_EH_SH = PresAbs_BR_SH_EH .- PresAbs_NB_SH_EH

migr_br_WH_NH = copy(migr_WH_NH)
migr_br_WH_NH[find(migr_br_WH_NH .== -1)] = 0
migr_nb_WH_NH = copy(migr_WH_NH)
migr_nb_WH_NH[find(migr_nb_WH_NH .== 1)] = 0
migr_nb_WH_NH = abs(migr_nb_WH_NH)
migr_br_WH_SH = copy(migr_WH_SH)
migr_br_WH_SH[find(migr_br_WH_SH .== -1)] = 0
migr_nb_WH_SH = copy(migr_WH_SH)
migr_nb_WH_SH[find(migr_nb_WH_SH .== 1)] = 0
migr_nb_WH_SH = abs(migr_nb_WH_SH)
migr_br_EH_NH = copy(migr_EH_NH)
migr_br_EH_NH[find(migr_br_EH_NH .== -1)] = 0
migr_nb_EH_NH = copy(migr_EH_NH)
migr_nb_EH_NH[find(migr_nb_EH_NH .== 1)] = 0
migr_nb_EH_NH = abs(migr_nb_EH_NH)
migr_br_EH_SH = copy(migr_EH_SH)
migr_br_EH_SH[find(migr_br_EH_SH .== -1)] = 0
migr_nb_EH_SH = copy(migr_EH_SH)
migr_nb_EH_SH[find(migr_nb_EH_SH .== 1)] = 0
migr_nb_EH_SH = abs(migr_nb_EH_SH)

perm_WH_NH = PresAbs_BR_NH_WH .* PresAbs_NB_NH_WH
perm_WH_SH = PresAbs_BR_SH_WH .* PresAbs_NB_SH_WH
perm_EH_NH = PresAbs_BR_NH_EH .* PresAbs_NB_NH_EH
perm_EH_SH = PresAbs_BR_SH_EH .* PresAbs_NB_SH_EH

breeding_WH = sum(migr_br_WH_NH, 2) + sum(migr_br_WH_SH, 2)
nonbreeding_WH = sum(migr_nb_WH_NH, 2) + sum(migr_nb_WH_SH, 2)
residentsObs_WH = sum(PresAbs_res_WH, 2) + sum(perm_WH_NH, 2) + sum(perm_WH_SH, 2)
prop_migr_obs_WH = (breeding_WH + nonbreeding_WH) ./ (breeding_WH + nonbreeding_WH + residentsObs_WH)
diff_obs_WH = breeding_WH - nonbreeding_WH

breeding_EH = sum(migr_br_EH_NH, 2) + sum(migr_br_EH_SH, 2)
nonbreeding_EH = sum(migr_nb_EH_NH, 2) + sum(migr_nb_EH_SH, 2)
residentsObs_EH = sum(PresAbs_res_EH, 2) + sum(perm_EH_NH, 2) + sum(perm_EH_SH, 2)
prop_migr_obs_EH = (breeding_EH + nonbreeding_EH) ./ (breeding_EH + nonbreeding_EH + residentsObs_EH)
diff_obs_EH = breeding_EH - nonbreeding_EH



### Simulate virtual ranges for a species

# Western Hemisphere
rangesize_sel = 131 #median(rangeSizes_WH)
rsimul = zeros(rangesize_sel)

j=1
while j < 400
	hex_1 = rand(1:1:2280)
	rangeSimulated = copy(hex_1)

	i = 1
	nicheneig_sel = 0
	while length(rangeSimulated) < rangesize_sel # & findin(true,isna(nicheneig_sel)) == false
    neig = sum(pairwiseDistanceWH2[:,rangeSimulated], 2)
    #neigs = [0]
    #for i in 1:length(neig)
    #  neigs = vcat(neigs, rep([i], Int(neig[i])))
    #end
    #neigs = neigs[2:end]
    neigs = findin(neig.>0, true)
    deleteat!(neigs, findin(neigs, rangeSimulated))

		# Climate driven selection of neighbour hexagons based on exponential decay
		if length(neigs) > 1
			nicheneig = euclidean([mean(temp_meanAnnual_std_WH[rangeSimulated]), mean(prec_meanAnnual_std_WH[rangeSimulated])], squeeze(reshape(hcat(temp_meanAnnual_std_WH[neigs], prec_meanAnnual_std_WH[neigs])[1,:], 2, 1), 2))
      for m in 2:length(neigs)
        nicheneig = hcat(nicheneig, euclidean([mean(temp_meanAnnual_std_WH[rangeSimulated]), mean(prec_meanAnnual_std_WH[rangeSimulated])], squeeze(reshape(hcat(temp_meanAnnual_std_WH[neigs], prec_meanAnnual_std_WH[neigs])[m,:], 2, 1), 2)))
      end

			if length(findin(isna(nicheneig), true)) > 0
				neigs = neigs[findin(isna(nicheneig), false)]
				nicheneig = nicheneig[findin(isna(nicheneig), false)]
      end

      beta = 3
			nicheneig_proba = exp(-beta .* nicheneig)
			nicheneig_proba = nicheneig_proba ./ sum(nicheneig_proba)
			nicheneig_sel = sample(neigs, WeightVec(squeeze(nicheneig_proba,1)), Int64(ceil(length(neigs)/4)), replace=false)
    end

    if length(neigs) == 1 nicheneig_sel = copy(neigs) end
		if length(neigs) == 0 nicheneig_sel = deleteat!(sortperm(pairwiseDistanceWH[:,hex_1]), findin(sortperm(pairwiseDistanceWH[:,hex_1]), rangeSimulated))[1] end

		rangeSimulated = vcat(rangeSimulated, nicheneig_sel)
		i = i+1
  end		# Get the neighbour hexagons and select the most suitable ones
	#rangeSimulated = rangeSimulated[which(is.na(range.simulated)==F)]
	rangeSimulated = rangeSimulated[1:rangesize_sel]

  rsimul = hcat(rsimul, rangeSimulated)
	j=j+1		# Generate potential ranges by expanding from a seed hexagon until 80% of the continent is occupied by at least one range
end
rsimul = rsimul[:,2:end]
rsimul = round(Int64, rsimul)

#writedlm("simulatedRangesWHnpp", rsimul, '\t')
rsimulWH = Array(readtable("simulatedRangesWHnpp.csv", header = false))



# Eastern Hemisphere
rangesize_sel = 180 #mean(rangeSizes_WH)
rsimul = zeros(rangesize_sel)

j=1
while j < 600
	hex_1 = rand(1:1:4969)
	rangeSimulated = copy(hex_1)

	i = 1
	nicheneig_sel = 0
	while length(rangeSimulated) < rangesize_sel # & findin(true,isna(nicheneig_sel)) == false
    neig = sum(pairwiseDistanceEH2[:,rangeSimulated], 2)
    #neigs = [0]
    #for i in 1:length(neig)
    #  neigs = vcat(neigs, rep([i], Int(neig[i])))
    #end
    #neigs = neigs[2:end]
    neigs = findin(neig.>0, true)
    deleteat!(neigs, findin(neigs, rangeSimulated))

		# Climate driven selection of neighbour hexagons based on exponential decay
		if length(neigs) > 1
			nicheneig = euclidean([mean(temp_meanAnnual_std_EH[rangeSimulated]), mean(prec_meanAnnual_std_EH[rangeSimulated])], squeeze(reshape(hcat(temp_meanAnnual_std_EH[neigs], prec_meanAnnual_std_EH[neigs])[1,:], 2, 1), 2))
      for m in 2:length(neigs)
        nicheneig = hcat(nicheneig, euclidean([mean(temp_meanAnnual_std_EH[rangeSimulated]), mean(prec_meanAnnual_std_EH[rangeSimulated])], squeeze(reshape(hcat(temp_meanAnnual_std_EH[neigs], prec_meanAnnual_std_EH[neigs])[m,:], 2, 1), 2)))
      end

			if length(findin(isna(nicheneig), true)) > 0
				neigs = neigs[findin(isna(nicheneig), false)]
				nicheneig = nicheneig[findin(isna(nicheneig), false)]
      end

      beta = 3
			nicheneig_proba = exp(-beta .* nicheneig)
			nicheneig_proba = nicheneig_proba ./ sum(nicheneig_proba)
			nicheneig_sel = sample(neigs, WeightVec(squeeze(nicheneig_proba,1)), Int64(ceil(length(neigs)/4)), replace=false)
    end

    if length(neigs) == 1 nicheneig_sel = copy(neigs) end
		if length(neigs) == 0 nicheneig_sel = deleteat!(sortperm(pairwiseDistanceEH[:,hex_1]), findin(sortperm(pairwiseDistanceEH[:,hex_1]), rangeSimulated))[1] end

		rangeSimulated = vcat(rangeSimulated, nicheneig_sel)
		i = i+1
  end		# Get the neighbour hexagons and select the most suitable ones
	#rangeSimulated = rangeSimulated[which(is.na(range.simulated)==F)]
	rangeSimulated = rangeSimulated[1:rangesize_sel]

  rsimul = hcat(rsimul, rangeSimulated)
	j=j+1		# Generate potential ranges by expanding from a seed hexagon until 80% of the continent is occupied by at least one range
end
rsimul = rsimul[:,2:end]
rsimul = round(Int64, rsimul)

#writedlm("simulatedRangesEHnpp", rsimul, '\t')
rsimulEH = Array(readtable("simulatedRangesEHnpp.csv", header = false))





### Parameters

mm = 0.00003 #0.0001016
rr = 0.6  # 0.3476501
tt = 25 #13.3983



## WESTERN HEMISPHERE

###  Costs associated with migration, thermoregulation and reproduction

migrationDistance_inKm = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
for i in 1:length(rsimulWH[1,:])
  for j in 1:length(rsimulWH[1,:])
      migrationDistance_inKm[i,j] = greatcircledistance(degtorad(mean(lonlatWH[rsimulWH[:,i],1])), degtorad(mean(lonlatWH[rsimulWH[:,i],2])), degtorad(mean(lonlatWH[rsimulWH[:,j],1])),degtorad(mean(lonlatWH[rsimulWH[:,j],2])))
  end
end
migrationCost = mm .* migrationDistance_inKm

thermoCost_NS = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
thermoCost_NW = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
for i in 1:length(rsimulWH[1,:])
  for j in 1:length(rsimulWH[1,:])
    thermoCost_NS[i,j] = (40 - tt - mean(tempNS_WH[rsimulWH[:,i]])) / tt
    thermoCost_NW[i,j] = (40 - tt - mean(tempNW_WH[rsimulWH[:,j]])) / tt
    if thermoCost_NS[i,j] < 0 thermoCost_NS[i,j] = 0 end
    if thermoCost_NW[i,j] < 0 thermoCost_NW[i,j] = 0 end
  end
end

migration_thermo_cost_brNS_NS = 1 .+ migrationCost .+ thermoCost_NS .+ rr
migration_thermo_cost_brNS_NW = 1 .+ migrationCost .+ thermoCost_NW
migration_thermo_cost_brNW_NS = 1 .+ migrationCost .+ thermoCost_NS
migration_thermo_cost_brNW_NW = 1 .+ migrationCost .+ thermoCost_NW .+ rr



S_July_WH = zeros(length(lonlatWH[:,1]))
S_January_WH = zeros(length(lonlatWH[:,1]))
ranges_selected_BR_WH = zeros(2280,1)
ranges_selected_NB_WH = zeros(2280,1)
sel_WH = [0]
selected_ranges = [0]
energyAvailable_NS = copy(energySupplyJulyWH)
energyAvailable_NW = copy(energySupplyJanuaryWH)

k=1
while sum(energyAvailable_NS) > (sum(energySupplyJulyWH) * 0.05) && sum(energyAvailable_NW) > (sum(energySupplyJanuaryWH) * 0.05) && sum(S_July_WH) < 500000  && sum(S_January_WH) < 500000

  supplyNS = zeros(length(rsimulWH[1,:]))
  for i in 1:length(rsimulWH[1,:])
    supplyNS[i] = mean(energyAvailable_NS[rsimulWH[:,i]])
    if supplyNS[i] == 0 supplyNS[i] = 0.1 end
  end
  supplyNW = zeros(length(rsimulWH[1,:]))
  for i in 1:length(rsimulWH[1,:])
    supplyNW[i] = mean(energyAvailable_NW[rsimulWH[:,i]])
    if supplyNW[i] == 0 supplyNW[i] = 0.1 end
  end

  state_proba_brNS = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
  state_proba_brNW = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
  for i in 1:length(rsimulWH[1,:])
    for j in 1:length(rsimulWH[1,:])
      #state_proba_brNS[i,j] = exp( - 10 * (migration_thermo_cost_brNS_NS[i,j] / supplyNS[i]) ) + exp( - 10 * (migration_thermo_cost_brNS_NW[i,j] / supplyNW[j]) )
      #state_proba_brNW[i,j] = exp( - 10 * (migration_thermo_cost_brNW_NS[i,j] / supplyNS[i]) ) + exp( - 10 * (migration_thermo_cost_brNW_NW[i,j] / supplyNW[j]) )
      state_proba_brNS[i,j] = (migration_thermo_cost_brNS_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNS_NW[i,j] / supplyNW[j])
      state_proba_brNW[i,j] = (migration_thermo_cost_brNW_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNW_NW[i,j] / supplyNW[j])
    end
  end

  # Find the most energy efficient strategy (i.e. the lowest value across both matrices)
  best_link1 = findmin(state_proba_brNS)
  best_link2 = findmin(state_proba_brNW)

  sele_WH = findmin([best_link1[1], best_link2[1]])[2]
  selected_link = [best_link1[2], best_link2[2]][sele_WH]
  sel_WH = vcat(sel_WH, sele_WH)
  selected_ranges = vcat(selected_ranges, selected_link)

  if sele_WH == 1
    best_link_indexes = findn(state_proba_brNS .== minimum(state_proba_brNS))
    rangeNB = best_link_indexes[2]
    rangeBR = best_link_indexes[1]
  else
    best_link_indexes = findn(state_proba_brNW .== minimum(state_proba_brNW))
    rangeNB = best_link_indexes[1]
    rangeBR = best_link_indexes[2]
  end

  #rangeNB = Int(ceil(selected_link/length(rsimul[1,:])))
  #rangeBR = selected_link%length(rsimul[1,:])
  #if rangeBR == 0 rangeBR = Int(selected_link/length(rsimul[1,:])) end
  ra = zeros(2280)
  ra[rsimulWH[:,rangeBR]] = 1
  ranges_selected_BR_WH = hcat(ranges_selected_BR_WH, ra)
  ra = zeros(2280)
  ra[rsimulWH[:,rangeNB]] = 1
  ranges_selected_NB_WH = hcat(ranges_selected_NB_WH, ra)

  if sele_WH == 1
    S_January_WH[ranges_selected_NB_WH[:,end] .== 1] = S_January_WH[ranges_selected_NB_WH[:,end] .== 1] .+ migration_thermo_cost_brNS_NW[rangeBR, rangeNB][1]
    S_July_WH[ranges_selected_BR_WH[:,end] .== 1] = S_July_WH[ranges_selected_BR_WH[:,end] .== 1] .+ migration_thermo_cost_brNS_NS[rangeBR, rangeNB][1]
  else
    S_January_WH[ranges_selected_BR_WH[:,end] .== 1] = S_January_WH[ranges_selected_BR_WH[:,end] .== 1] .+ migration_thermo_cost_brNW_NW[rangeNB, rangeBR][1]
    S_July_WH[ranges_selected_NB_WH[:,end] .== 1] = S_July_WH[ranges_selected_NB_WH[:,end] .== 1] .+ migration_thermo_cost_brNW_NS[rangeNB, rangeBR][1]
  end

  energyAvailable_NS = energySupplyJulyWH .- S_July_WH
  energyAvailable_NW = energySupplyJanuaryWH .- S_January_WH
  for i in 1:length(energyAvailable_NS)
    if energyAvailable_NS[i] < 0 energyAvailable_NS[i] = 0 end
    if energyAvailable_NW[i] < 0 energyAvailable_NW[i] = 0 end
  end

   k=k+1
 end

 ranges_selected_BR_WH = ranges_selected_BR_WH[:,2:end]
 ranges_selected_NB_WH = ranges_selected_NB_WH[:,2:end]

 migra = ranges_selected_BR_WH .- ranges_selected_NB_WH
 res = ranges_selected_BR_WH .+ ranges_selected_NB_WH
 migraBR = zeros(length(migra[:,1]), length(migra[1,:]))
 migraNB = zeros(length(migra[:,1]), length(migra[1,:]))
 residents = zeros(length(migra[:,1]), length(migra[1,:]))
 for i in 1:length(ranges_selected_BR_WH[:,1])
   for j in 1:length(ranges_selected_BR_WH[1,:])
     if migra[i,j] > 0 migraBR[i,j] = 1 end
     if migra[i,j] < 0 migraNB[i,j] = 1 end
     if res[i,j] == 2 residents[i,j] = 1 end
   end
end

richnessBRsimuWH = zeros(length(ranges_selected_BR_WH[:,1]))
richnessNBsimuWH = zeros(length(ranges_selected_NB_WH[:,1]))
richnessMigraBRsimuWH = zeros(length(ranges_selected_BR_WH[:,1]))
richnessMigraNBsimuWH = zeros(length(ranges_selected_NB_WH[:,1]))
richnessResidentsimuWH = zeros(length(ranges_selected_NB_WH[:,1]))
for i in 1:length(ranges_selected_BR_WH[:,1])
  richnessBRsimuWH[i] = sum(ranges_selected_BR_WH[i,:])
  richnessNBsimuWH[i] = sum(ranges_selected_NB_WH[i,:])
  richnessMigraBRsimuWH[i] = sum(migraBR[i,:])
  richnessMigraNBsimuWH[i] = sum(migraNB[i,:])
richnessResidentsimuWH[i] = sum(residents[i,:])
end
richnessDiffSimuWH = richnessBRsimuWH .- richnessNBsimuWH
PropMigrSimuWH = (richnessMigraBRsimuWH .+ richnessMigraNBsimuWH) ./ (richnessMigraBRsimuWH .+ richnessMigraNBsimuWH .+ richnessResidentsimuWH)


plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessBRsimuWH, style(default_point_size=1.8pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessNBsimuWH, style(default_point_size=1.8pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessMigraBRsimuWH, style(default_point_size=1.8pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessMigraNBsimuWH, style(default_point_size=1.8pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessResidentsimuWH, style(default_point_size=1.8pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessDiffSimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=Int(round((PropMigrSimu*100), 0)), style(default_point_size=1.3pt))

plot(x=lonlatWH[:,2], y=richnessMigraBRsimuWH, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,2], y=richnessMigraNBsimuWH, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,2], y=richnessDiffSimuWH, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,2], y=PropMigrSimuWH, style(default_point_size=1.3pt))

plot(x=lonlatWH[:,2], y=breeding_WH, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,2], y=nonbreeding_WH, style(default_point_size=1.3pt))




## EASTERN HEMISPHERE

  migrationDistance_inKm = zeros(length(rsimulEH[1,:]), length(rsimulEH[1,:]))
  for i in 1:length(rsimulEH[1,:])
    for j in 1:length(rsimulEH[1,:])
        migrationDistance_inKm[i,j] = greatcircledistance(degtorad(mean(lonlatEH[rsimulEH[:,i],1])), degtorad(mean(lonlatEH[rsimulEH[:,i],2])), degtorad(mean(lonlatEH[rsimulEH[:,j],1])),degtorad(mean(lonlatEH[rsimulEH[:,j],2])))
    end
  end
  migrationCost = mm .* migrationDistance_inKm

  thermoCost_NS = zeros(length(rsimulEH[1,:]), length(rsimulEH[1,:]))
  thermoCost_NW = zeros(length(rsimulEH[1,:]), length(rsimulEH[1,:]))
  for i in 1:length(rsimulEH[1,:])
    for j in 1:length(rsimulEH[1,:])
      thermoCost_NS[i,j] = (40 - tt - mean(tempNS_EH[rsimulEH[:,i]])) / tt
      thermoCost_NW[i,j] = (40 - tt - mean(tempNW_EH[rsimulEH[:,j]])) / tt
      if thermoCost_NS[i,j] < 0 thermoCost_NS[i,j] = 0 end
      if thermoCost_NW[i,j] < 0 thermoCost_NW[i,j] = 0 end
    end
  end

  migration_thermo_cost_brNS_NS = 1 .+ migrationCost .+ thermoCost_NS .+ rr
  migration_thermo_cost_brNS_NW = 1 .+ migrationCost .+ thermoCost_NW
  migration_thermo_cost_brNW_NS = 1 .+ migrationCost .+ thermoCost_NS
  migration_thermo_cost_brNW_NW = 1 .+ migrationCost .+ thermoCost_NW .+ rr


  S_July_EH = zeros(length(lonlatEH[:,1]))
  S_January_EH = zeros(length(lonlatEH[:,1]))
  ranges_selected_BR_EH = zeros(4969,1)
  ranges_selected_NB_EH = zeros(4969,1)
  sel_EH = [0]
  selected_ranges = [0]
  energyAvailable_NS = copy(energySupplyJulyEH)
  energyAvailable_NW = copy(energySupplyJanuaryEH)

  k=1
  while sum(energyAvailable_NS) > (sum(energySupplyJulyEH) * 0.05) && sum(energyAvailable_NW) > (sum(energySupplyJanuaryEH) * 0.05)  && sum(S_July_EH) < 750000  && sum(S_January_EH) < 750000

    supplyNS = zeros(length(rsimulEH[1,:]))
    for i in 1:length(rsimulEH[1,:])
      supplyNS[i] = mean(energyAvailable_NS[rsimulEH[:,i]])
      if supplyNS[i] == 0 supplyNS[i] = 0.1 end
    end
    supplyNW = zeros(length(rsimulEH[1,:]))
    for i in 1:length(rsimulEH[1,:])
      supplyNW[i] = mean(energyAvailable_NW[rsimulEH[:,i]])
      if supplyNW[i] == 0 supplyNW[i] = 0.1 end
    end

    state_proba_brNS = zeros(length(rsimulEH[1,:]), length(rsimulEH[1,:]))
    state_proba_brNW = zeros(length(rsimulEH[1,:]), length(rsimulEH[1,:]))
    for i in 1:length(rsimulEH[1,:])
      for j in 1:length(rsimulEH[1,:])
        #state_proba_brNS[i,j] = exp( - ((migration_thermo_cost_brNS_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNS_NW[i,j] / supplyNW[j])) )
        #state_proba_brNW[i,j] = exp( - ((migration_thermo_cost_brNW_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNW_NW[i,j] / supplyNW[j])) )
        state_proba_brNS[i,j] = (migration_thermo_cost_brNS_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNS_NW[i,j] / supplyNW[j])
        state_proba_brNW[i,j] = (migration_thermo_cost_brNW_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNW_NW[i,j] / supplyNW[j])
      end
    end

    # Find the most energy efficient strategy (i.e. the lowest value across both matrices)
    best_link1 = findmin(state_proba_brNS)
    best_link2 = findmin(state_proba_brNW)

    sele_EH = findmin([best_link1[1], best_link2[1]])[2]
    selected_link = [best_link1[2], best_link2[2]][sele_EH]
    sel_EH = vcat(sel_EH, sele_EH)
    selected_ranges = vcat(selected_ranges, selected_link)

    if sele_EH == 1
      best_link_indexes = findn(state_proba_brNS .== minimum(state_proba_brNS))
      rangeNB = best_link_indexes[2]
      rangeBR = best_link_indexes[1]
    else
      best_link_indexes = findn(state_proba_brNW .== minimum(state_proba_brNW))
      rangeNB = best_link_indexes[1]
      rangeBR = best_link_indexes[2]
    end

    ra = zeros(4969)
    ra[rsimulEH[:,rangeBR]] = 1
    ranges_selected_BR_EH = hcat(ranges_selected_BR_EH, ra)
    ra = zeros(4969)
    ra[rsimulEH[:,rangeNB]] = 1
    ranges_selected_NB_EH = hcat(ranges_selected_NB_EH, ra)

    if sele_EH == 1
      S_January_EH[ranges_selected_NB_EH[:,end] .== 1] = S_January_EH[ranges_selected_NB_EH[:,end] .== 1] .+ migration_thermo_cost_brNS_NW[rangeBR, rangeNB][1]
      S_July_EH[ranges_selected_BR_EH[:,end] .== 1] = S_July_EH[ranges_selected_BR_EH[:,end] .== 1] .+ migration_thermo_cost_brNS_NS[rangeBR, rangeNB][1]
    else
      S_January_EH[ranges_selected_BR_EH[:,end] .== 1] = S_January_EH[ranges_selected_BR_EH[:,end] .== 1] .+ migration_thermo_cost_brNW_NW[rangeNB, rangeBR][1]
      S_July_EH[ranges_selected_NB_EH[:,end] .== 1] = S_July_EH[ranges_selected_NB_EH[:,end] .== 1] .+ migration_thermo_cost_brNW_NS[rangeNB, rangeBR][1]
    end

    energyAvailable_NS = energySupplyJulyEH .- S_July_EH
    energyAvailable_NW = energySupplyJanuaryEH .- S_January_EH
    for i in 1:length(energyAvailable_NS)
      if energyAvailable_NS[i] < 0 energyAvailable_NS[i] = 0 end
      if energyAvailable_NW[i] < 0 energyAvailable_NW[i] = 0 end
    end

	  k=k+1
   end

   ranges_selected_BR_EH = ranges_selected_BR_EH[:,2:end]
   ranges_selected_NB_EH = ranges_selected_NB_EH[:,2:end]

   migra = ranges_selected_BR_EH .- ranges_selected_NB_EH
   res = ranges_selected_BR_EH .+ ranges_selected_NB_EH
   migraBR = zeros(length(migra[:,1]), length(migra[1,:]))
   migraNB = zeros(length(migra[:,1]), length(migra[1,:]))
   residents = zeros(length(migra[:,1]), length(migra[1,:]))
   for i in 1:length(ranges_selected_BR_EH[:,1])
     for j in 1:length(ranges_selected_BR_EH[1,:])
       if migra[i,j] > 0 migraBR[i,j] = 1 end
       if migra[i,j] < 0 migraNB[i,j] = 1 end
       if res[i,j] == 2 residents[i,j] = 1 end
     end
   end

  richnessBRsimuEH = zeros(length(ranges_selected_BR_EH[:,1]))
  richnessNBsimuEH = zeros(length(ranges_selected_NB_EH[:,1]))
  richnessMigraBRsimuEH = zeros(length(ranges_selected_BR_EH[:,1]))
  richnessMigraNBsimuEH = zeros(length(ranges_selected_NB_EH[:,1]))
  richnessResidentsimuEH = zeros(length(ranges_selected_NB_EH[:,1]))
  for i in 1:length(ranges_selected_BR_EH[:,1])
    richnessBRsimuEH[i] = sum(ranges_selected_BR_EH[i,:])
    richnessNBsimuEH[i] = sum(ranges_selected_NB_EH[i,:])
    richnessMigraBRsimuEH[i] = sum(migraBR[i,:])
    richnessMigraNBsimuEH[i] = sum(migraNB[i,:])
    richnessResidentsimuEH[i] = sum(residents[i,:])
  end
  richnessDiffSimuEH = richnessBRsimuEH .- richnessNBsimuEH
  PropMigrSimuEH = (richnessMigraBRsimuEH .+ richnessMigraNBsimuEH) ./ (richnessMigraBRsimuEH .+ richnessMigraNBsimuEH .+ richnessResidentsimuEH)


  plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=vcat(richnessBRsimuWH,richnessBRsimuEH), style(default_point_size=1.8pt))
  plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=vcat(richnessNBsimuWH,richnessNBsimuEH), style(default_point_size=1.8pt))
  plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=vcat(richnessMigraBRsimuWH,richnessMigraBRsimuEH), style(default_point_size=1.8pt))
  plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=vcat(richnessMigraNBsimuWH,richnessMigraNBsimuEH), style(default_point_size=1.8pt))
  plot(x=vcat(lonlatWH[:,1],lonlatEH[:,1]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), color=vcat(richnessResidentsimuWH,richnessResidentsimuEH), style(default_point_size=1.8pt))

  plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=vcat(lonlatWH[:,2],lonlatEH[:,2]), style(default_point_size=1.3pt))
  plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=vcat(richnessMigraNBsimuWH,richnessMigraNBsimuEH), style(default_point_size=1.3pt))
  plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=vcat(richnessDiffSimuWH,richnessDiffSimuEH), style(default_point_size=1.3pt))
  plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=vcat(PropMigrSimuWH,PropMigrSimuEH), style(default_point_size=1.3pt))

  plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=breeding_WH, style(default_point_size=1.3pt))
  plot(x=vcat(lonlatWH[:,2],lonlatEH[:,2]), y=nonbreeding_WH, style(default_point_size=1.3pt))










###  Costs associated with migration, thermoregulation and reproduction

migrationDistance_inKm = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
for i in 1:length(rsimulWH[1,:])
  for j in 1:length(rsimulWH[1,:])
    migrationDistance_inKm[i,j] = greatcircledistance(degtorad(mean(lonlatWH[rsimulWH[:,i],1])), degtorad(mean(lonlatWH[rsimulWH[:,i],2])), degtorad(mean(lonlatWH[rsimulWH[:,j],1])),degtorad(mean(lonlatWH[rsimulWH[:,j],2])))
  end
end
migrationCost = mm .* migrationDistance_inKm

thermoCost_NS = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
thermoCost_NW = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
for i in 1:length(rsimulWH[1,:])
  for j in 1:length(rsimulWH[1,:])
    thermoCost_NS[i,j] = (maximum([mean(tempNS_WH[rsimulWH[:,i]]), mean(tempNW_WH[rsimulWH[:,j]])]) - mean(tempNS_WH[rsimulWH[:,i]])) / tt
    thermoCost_NW[i,j] = (maximum([mean(tempNS_WH[rsimulWH[:,i]]), mean(tempNW_WH[rsimulWH[:,j]])]) - mean(tempNW_WH[rsimulWH[:,j]])) / tt
    if thermoCost_NS[i,j] < 0 thermoCost_NS[i,j] = 0 end
    if thermoCost_NW[i,j] < 0 thermoCost_NW[i,j] = 0 end
  end
end

migration_thermo_cost_brNS_NS = 1 .+ migrationCost .+ thermoCost_NS .+ rr
migration_thermo_cost_brNS_NW = 1 .+ migrationCost .+ thermoCost_NW
migration_thermo_cost_brNW_NS = 1 .+ migrationCost .+ thermoCost_NS
migration_thermo_cost_brNW_NW = 1 .+ migrationCost .+ thermoCost_NW .+ rr

migration_thermo_cost_brNS = migration_thermo_cost_brNS_NS .+ migration_thermo_cost_brNS_NW
migration_thermo_cost_brNW = migration_thermo_cost_brNW_NS .+ migration_thermo_cost_brNW_NW



### RUN THE MODEL

S_July_WH = zeros(length(lonlatWH[:,1]))
S_January_WH = zeros(length(lonlatWH[:,1]))
ranges_selected_BR_WH = zeros(2280,1)
ranges_selected_NB_WH = zeros(2280,1)
sel_WH = [0]
selected_ranges = [0]
energyAvailable_NS = copy(energySupplyJulyWH)
energyAvailable_NW = copy(energySupplyJanuaryWH)

k=1
while sum(energyAvailable_NS) > 5000 && sum(energyAvailable_NW) > 5000 #sum(S_January_WH) < sum(energySupplyJanuaryWH) && sum(S_July_WH) < sum(energySupplyJulyWH)

 supplyNS = zeros(length(rsimulWH[1,:]))
 for i in 1:length(rsimulWH[1,:])
   supplyNS[i] = mean(energyAvailable_NS[rsimulWH[:,i]])
   if supplyNS[i] == 0 supplyNS[i] = 0.1 end
 end
 supplyNW = zeros(length(rsimulWH[1,:]))
 for i in 1:length(rsimulWH[1,:])
   supplyNW[i] = mean(energyAvailable_NW[rsimulWH[:,i]])
   if supplyNW[i] == 0 supplyNW[i] = 0.1 end
 end

 state_proba_brNS = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
 state_proba_brNW = zeros(length(rsimulWH[1,:]), length(rsimulWH[1,:]))
 for i in 1:length(rsimulWH[1,:])
   for j in 1:length(rsimulWH[1,:])
     state_proba_brNS[i,j] = exp( - ((migration_thermo_cost_brNS_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNS_NW[i,j] / supplyNW[j])) )
     state_proba_brNW[i,j] = exp( - ((migration_thermo_cost_brNW_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNW_NW[i,j] / supplyNW[j])) )
   end
 end

 # Find the most energy efficient strategy (i.e. the lowest value across both matrices)
 best_link1 = findmax(state_proba_brNS)
 best_link2 = findmax(state_proba_brNW)

 sele_WH = findmax([best_link1[1], best_link2[1]])[2]
 selected_link = [best_link1[2], best_link2[2]][sele_WH]
 sel_WH = vcat(sel_WH, sele_WH)
 selected_ranges = vcat(selected_ranges, selected_link)

 if sele_WH == 1
   best_link_indexes = findn(state_proba_brNS .== maximum(state_proba_brNS))
   rangeNB = best_link_indexes[2]
   rangeBR = best_link_indexes[1]
 else
   best_link_indexes = findn(state_proba_brNW .== maximum(state_proba_brNW))
   rangeNB = best_link_indexes[1]
   rangeBR = best_link_indexes[2]
 end

 #rangeNB = Int(ceil(selected_link/length(rsimul[1,:])))
 #rangeBR = selected_link%length(rsimul[1,:])
 #if rangeBR == 0 rangeBR = Int(selected_link/length(rsimul[1,:])) end
 ra = zeros(2280)
 ra[rsimul[:,rangeBR]] = 1
 ranges_selected_BR_WH = hcat(ranges_selected_BR_WH, ra)
 ra = zeros(2280)
 ra[rsimul[:,rangeNB]] = 1
 ranges_selected_NB_WH = hcat(ranges_selected_NB_WH, ra)

 if sele_WH == 1
   S_January_WH[ranges_selected_NB_WH[:,end] .== 1] = S_January_WH[ranges_selected_NB_WH[:,end] .== 1] .+ migration_thermo_cost_brNS_NW[rangeBR, rangeNB][1]
   S_July_WH[ranges_selected_BR_WH[:,end] .== 1] = S_July_WH[ranges_selected_BR_WH[:,end] .== 1] .+ migration_thermo_cost_brNS_NS[rangeBR, rangeNB][1]
 else
   S_January_WH[ranges_selected_BR_WH[:,end] .== 1] = S_January_WH[ranges_selected_BR_WH[:,end] .== 1] .+ migration_thermo_cost_brNW_NW[rangeNB, rangeBR][1]
   S_July_WH[ranges_selected_NB_WH[:,end] .== 1] = S_July_WH[ranges_selected_NB_WH[:,end] .== 1] .+ migration_thermo_cost_brNW_NS[rangeNB, rangeBR][1]
 end

 energyAvailable_NS = energySupplyJulyWH .- S_July_WH
 energyAvailable_NW = energySupplyJanuaryWH .- S_January_WH
 for i in 1:length(energyAvailable_NS)
   if energyAvailable_NS[i] < 0 energyAvailable_NS[i] = 0 end
   if energyAvailable_NW[i] < 0 energyAvailable_NW[i] = 0 end
 end

	k=k+1
end


ranges_selected_BR_WH = ranges_selected_BR_WH[:,2:end]
ranges_selected_NB_WH = ranges_selected_NB_WH[:,2:end]
sel_WH = sel_WH[2:end]
#selected_ranges = selected_ranges[2:end]

migra = ranges_selected_BR_WH .- ranges_selected_NB_WH
res = ranges_selected_BR_WH .+ ranges_selected_NB_WH
migraBR = zeros(length(migra[:,1]), length(migra[1,:]))
migraNB = zeros(length(migra[:,1]), length(migra[1,:]))
residents = zeros(length(migra[:,1]), length(migra[1,:]))
for i in 1:length(ranges_selected_BR_WH[:,1])
 for j in 1:length(ranges_selected_BR_WH[1,:])
   if migra[i,j] > 0 migraBR[i,j] = 1 end
   if migra[i,j] < 0 migraNB[i,j] = 1 end
   if res[i,j] == 2 residents[i,j] = 1 end
 end
end

richnessBRsimu = zeros(length(ranges_selected_BR_WH[:,1]))
richnessNBsimu = zeros(length(ranges_selected_NB_WH[:,1]))
richnessMigraBRsimu = zeros(length(ranges_selected_BR_WH[:,1]))
richnessMigraNBsimu = zeros(length(ranges_selected_NB_WH[:,1]))
richnessResidentsimu = zeros(length(ranges_selected_NB_WH[:,1]))
for i in 1:length(ranges_selected_BR_WH[:,1])
 richnessBRsimu[i] = sum(ranges_selected_BR_WH[i,:])
 richnessNBsimu[i] = sum(ranges_selected_NB_WH[i,:])
 richnessMigraBRsimu[i] = sum(migraBR[i,:])
 richnessMigraNBsimu[i] = sum(migraNB[i,:])
 richnessResidentsimu[i] = sum(residents[i,:])
end
richnessDiffSimu = richnessBRsimu .- richnessNBsimu
PropMigrSimu = (richnessMigraBRsimu .+ richnessMigraNBsimu) ./ (richnessMigraBRsimu .+ richnessMigraNBsimu .+ richnessResidentsimu)

plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessBRsimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessNBsimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessMigraBRsimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessMigraNBsimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessResidentsimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessDiffSimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=Int(round((PropMigrSimu*100), 0)), style(default_point_size=1.3pt))
maximum(richnessDiffSimu)

plot(x=lonlatWH[:,2], y=energyAvailable_NS, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,2], y=energyAvailable_NW, style(default_point_size=1.3pt))

plot(x=lonlatWH[:,2], y=S_July_WH)
plot(x=lonlatWH[:,2], y=S_January_WH)


plot(x=lonlatWH[:,2], y=richnessMigraBRsimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,2], y=richnessMigraNBsimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,2], y=richnessDiffSimu, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,2], y=PropMigrSimu, style(default_point_size=1.3pt))

plot(x=lonlatWH[:,2], y=breeding_WH, style(default_point_size=1.3pt))
plot(x=lonlatWH[:,2], y=nonbreeding_WH, style(default_point_size=1.3pt))
## Plot observed patterns
# Breeding season
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=breeding_WH, style(default_point_size=1.3pt))
#plot(x=lonlatWH[:,2], y=breeding_WH)
# Non-breeding season
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=nonbreeding_WH, style(default_point_size=1.3pt))
#plot(x=lonlatWH[:,2], y=nonbreeding_WH)
#plot(x=lonlatWH[:,2], y=diff_obs_WH)


plot(x=richnessMigraBRsimu, y=breeding_WH, style(default_point_size=1.3pt))
cor(richnessMigraBRsimu, breeding_WH)
plot(x=richnessMigraNBsimu, y=nonbreeding_WH, style(default_point_size=1.3pt))
cor(richnessMigraNBsimu, nonbreeding_WH)


# Total proportion of migrants

numberMigraBRsimu = zeros(length(ranges_selected_BR_WH[1,:]))
numberMigraNBsimu = zeros(length(ranges_selected_NB_WH[1,:]))
for i in 1:length(ranges_selected_BR_WH[1,:])
 #richnessBRsimu[i] = sum(ranges_selected_BR_WH[:,1])
 #richnessNBsimu[i] = sum(ranges_selected_NB_WH[:,1])
 numberMigraBRsimu[i] = sum(migraBR[:,i])
 if numberMigraBRsimu[i] > 0 numberMigraBRsimu[i] = 1 end
 numberMigraNBsimu[i] = sum(migraNB[:,i])
 if numberMigraNBsimu[i] > 0 numberMigraNBsimu[i] = 1 end
end
sum(numberMigraBRsimu) / k
