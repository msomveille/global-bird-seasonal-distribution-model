#Pkg.add("DataFrames")
using DataFrames
#Pkg.add("Gadfly")
using Gadfly
#Pkg.add("Distances")
using Distances
#Pkg.add("StatsBase")
using StatsBase

Gadfly.push_theme(:dark)

cd("$(homedir())/Desktop/PhD/Chapter 4 – mechanistic model/Loading-data-NDVI")
pwd()

PresAbs_BR_NH = readtable("PresAbs_BR_NH.csv") ## dim = 7613, 1185
PresAbs_NB_NH = readtable("PresAbs_NB_NH.csv")
PresAbs_BR_SH = readtable("PresAbs_BR_SH.csv") ## dim = 7613, 256
PresAbs_NB_SH = readtable("PresAbs_NB_SH.csv")
PresAbs_res = readtable("PresAbs_res.csv")

# Get climate data
climateDataset = readtable("env_data_6months.csv")
climateData = climateDataset[(climateDataset[:Temp_NS] .!= 0) | (climateDataset[:Prec_NS] .!= 0) | (climateDataset[:Temp_NW] .!= 0) | (climateDataset[:Prec_NW] .!= 0), :]

PresAbs_BR_NH = PresAbs_BR_NH[(climateDataset[:Temp_NS] .!= 0) | (climateDataset[:Prec_NS] .!= 0) | (climateDataset[:Temp_NW] .!= 0) | (climateDataset[:Prec_NW] .!= 0), :]
PresAbs_NB_NH = PresAbs_NB_NH[(climateDataset[:Temp_NS] .!= 0) | (climateDataset[:Prec_NS] .!= 0) | (climateDataset[:Temp_NW] .!= 0) | (climateDataset[:Prec_NW] .!= 0), :]
PresAbs_BR_SH = PresAbs_BR_SH[(climateDataset[:Temp_NS] .!= 0) | (climateDataset[:Prec_NS] .!= 0) | (climateDataset[:Temp_NW] .!= 0) | (climateDataset[:Prec_NW] .!= 0), :]
PresAbs_NB_SH = PresAbs_NB_SH[(climateDataset[:Temp_NS] .!= 0) | (climateDataset[:Prec_NS] .!= 0) | (climateDataset[:Temp_NW] .!= 0) | (climateDataset[:Prec_NW] .!= 0), :]
PresAbs_res = PresAbs_res[(climateDataset[:Temp_NS] .!= 0) | (climateDataset[:Prec_NS] .!= 0) | (climateDataset[:Temp_NW] .!= 0) | (climateDataset[:Prec_NW] .!= 0), :]

## Subset the global hexagon grid (with species' presence/absence) into Western and Eastern Hemispheres
# Western-Northern Hemisphere
PresAbs_BR_NH_WH = PresAbs_BR_NH[climateData[:LONGITUDE] .<= -30, 2:1185]
PresAbs_NB_NH_WH = PresAbs_NB_NH[climateData[:LONGITUDE] .<= -30, 2:1185]
spprangeBR = zeros(1184)
spprangeNB = zeros(1184)
for i in 1:1184
  spprangeBR[i] = sum(array(PresAbs_BR_NH_WH[:,i]))
  spprangeNB[i] = sum(array(PresAbs_NB_NH_WH[:,i]))
end
speciestokeep = (spprangeBR .> 0) & (spprangeNB .> 0)
PresAbs_BR_NH_WH = PresAbs_BR_NH_WH[:,speciestokeep]
PresAbs_NB_NH_WH = PresAbs_NB_NH_WH[:,speciestokeep]
# Western-Southern Hemisphere
PresAbs_BR_SH_WH = PresAbs_BR_SH[climateData[:LONGITUDE] .<= -30, 2:256]
PresAbs_NB_SH_WH = PresAbs_NB_SH[climateData[:LONGITUDE] .<= -30, 2:256]
spprangeBR = zeros(255)
spprangeNB = zeros(255)
for i in 1:255
  spprangeBR[i] = sum(array(PresAbs_BR_SH_WH[:,i]))
  spprangeNB[i] = sum(array(PresAbs_NB_SH_WH[:,i]))
end
speciestokeep = (spprangeBR .> 0) & (spprangeNB .> 0)
PresAbs_BR_SH_WH = PresAbs_BR_SH_WH[:,speciestokeep]
PresAbs_NB_SH_WH = PresAbs_NB_SH_WH[:,speciestokeep]
# Eastern-Northern Hemisphere
PresAbs_BR_NH_EH = PresAbs_BR_NH[climateData[:LONGITUDE] .> -30, 2:1185]
PresAbs_NB_NH_EH = PresAbs_NB_NH[climateData[:LONGITUDE] .> -30, 2:1185]
spprangeBR = zeros(1184)
spprangeNB = zeros(1184)
for i in 1:1184
  spprangeBR[i] = sum(array(PresAbs_BR_NH_EH[:,i]))
  spprangeNB[i] = sum(array(PresAbs_NB_NH_EH[:,i]))
end
speciestokeep = (spprangeBR .> 0) & (spprangeNB .> 0)
PresAbs_BR_NH_EH = PresAbs_BR_NH_EH[:,speciestokeep]
PresAbs_NB_NH_EH = PresAbs_NB_NH_EH[:,speciestokeep]
# Eastern-Southern Hemisphere
PresAbs_BR_SH_EH = PresAbs_BR_SH[climateData[:LONGITUDE] .> -30, 2:256]
PresAbs_NB_SH_EH = PresAbs_NB_SH[climateData[:LONGITUDE] .> -30, 2:256]
spprangeBR = zeros(255)
spprangeNB = zeros(255)
for i in 1:255
  spprangeBR[i] = sum(array(PresAbs_BR_SH_EH[:,i]))
  spprangeNB[i] = sum(array(PresAbs_NB_SH_EH[:,i]))
end
speciestokeep = (spprangeBR .> 0) & (spprangeNB .> 0)
PresAbs_BR_SH_EH = PresAbs_BR_SH_EH[:,speciestokeep]
PresAbs_NB_SH_EH = PresAbs_NB_SH_EH[:,speciestokeep]

# Residents Western Hemisphere
PresAbs_res_WH = PresAbs_res[climateData[:LONGITUDE] .<= -30, 2:7820]
spprangeRes = zeros(7819)
for i in 1:7819
  spprangeRes[i] = sum(array(PresAbs_res_WH[:,i]))
end
speciestokeep = (spprangeRes .> 0)
PresAbs_res_WH = PresAbs_res_WH[:,speciestokeep]
# Residents Eastern Hemisphere
PresAbs_res_EH = PresAbs_res[climateData[:LONGITUDE] .> -30, 2:7820]
spprangeRes = zeros(7819)
for i in 1:7819
  spprangeRes[i] = sum(array(PresAbs_res_EH[:,i]))
end
speciestokeep = (spprangeRes .> 0)
PresAbs_res_EH = PresAbs_res_EH[:,speciestokeep]


# From climate dataset, get climate and NDVI data for Northern Summer (NS; April to September) and Northern Winter (NW; November to March)
precNS = log(climateData[:Prec_NS]+1)
precNW = log(climateData[:Prec_NW]+1)
tempNS = climateData[:Temp_NS]/10
tempNW = climateData[:Temp_NW]/10
ndviNS = log(climateData[:NDVI_NS]+1)
ndviNW = log(climateData[:NDVI_NW]+1)

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
tempNS_std_WH = tempNS_std[climateData[:LONGITUDE] .< -30, :]
tempNW_std_WH = tempNW_std[climateData[:LONGITUDE] .< -30, :]
precNS_std_WH = precNS_std[climateData[:LONGITUDE] .< -30, :]
precNW_std_WH = precNW_std[climateData[:LONGITUDE] .< -30, :]
tempNS_WH = tempNS[climateData[:LONGITUDE] .< -30, :]
tempNW_WH = tempNW[climateData[:LONGITUDE] .< -30, :]
precNS_WH = precNS[climateData[:LONGITUDE] .< -30, :]
precNW_WH = precNW[climateData[:LONGITUDE] .< -30, :]
NDVI_ns_WH = ndviNS[climateData[:LONGITUDE] .< -30, :]
NDVI_nw_WH = ndviNW[climateData[:LONGITUDE] .< -30, :]
tempNS_std_EH = tempNS_std[climateData[:LONGITUDE] .>= -30, :]
tempNW_std_EH = tempNW_std[climateData[:LONGITUDE] .>= -30, :]
precNS_std_EH = precNS_std[climateData[:LONGITUDE] .>= -30, :]
precNW_std_EH = precNW_std[climateData[:LONGITUDE] .>= -30, :]
tempNS_EH = tempNS[climateData[:LONGITUDE] .>= -30, :]
tempNW_EH = tempNW[climateData[:LONGITUDE] .>= -30, :]
precNS_EH = precNS[climateData[:LONGITUDE] .>= -30, :]
precNW_EH = precNW[climateData[:LONGITUDE] .>= -30, :]
NDVI_ns_EH = ndviNS[climateData[:LONGITUDE] .>= -30, :]
NDVI_nw_EH = ndviNW[climateData[:LONGITUDE] .>= -30, :]

temp_meanAnnual_std_WH = (tempNS_std_WH .+ tempNW_std_WH) ./ 2
prec_meanAnnual_std_WH = (precNS_std_WH .+ precNW_std_WH) ./ 2
temp_meanAnnual_std_EH = (tempNS_std_EH .+ tempNW_std_EH) ./ 2
prec_meanAnnual_std_EH = (precNS_std_EH .+ precNW_std_EH) ./ 2

energySupplyJulyWH = 82 * NDVI_ns_WH
energySupplyJanuaryWH = 82 * NDVI_nw_WH
energySupplyJulyEH = 82 * NDVI_ns_EH
energySupplyJanuaryEH = 82 * NDVI_nw_EH

lonlat = climateData[:,[16,15]]
lonlatWH = lonlat[climateData[:LONGITUDE] .< -30, :]
lonlatEH = lonlat[climateData[:LONGITUDE] .>= -30, :]

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

pairwiseDistanceWH = readtable("pairwise_dist_WH3.csv")
pairwiseDistanceWH = pairwiseDistanceWH[:,2:ncol(pairwiseDistanceWH)]
pairwiseDistanceWH = array(pairwiseDistanceWH)
pairwiseDistanceWH2 = copy(pairwiseDistanceWH)
pairwiseDistanceWH2[(pairwiseDistanceWH2 .< 250) & (pairwiseDistanceWH2 .> 0)] = 1
pairwiseDistanceWH2[(pairwiseDistanceWH2 .>= 250) | (pairwiseDistanceWH2 .<= 0)] = 0

pairwiseDistanceEH = readtable("pairwise_dist_EH3.csv")
pairwiseDistanceEH = pairwiseDistanceEH[:,2:ncol(pairwiseDistanceEH)]
pairwiseDistanceEH = array(pairwiseDistanceEH)
pairwiseDistanceEH2 = copy(pairwiseDistanceEH)
pairwiseDistanceEH2[(pairwiseDistanceEH2 .< 250) & (pairwiseDistanceEH2 .> 0)] = 1
pairwiseDistanceEH2[(pairwiseDistanceEH2 .>= 250) | (pairwiseDistanceEH2 .<= 0)] = 0


### Observed patterns
PresAbs_BR_NH_WH2 = PresAbs_BR_NH_WH[:,2:end][:, find((sum(array(PresAbs_NB_NH_WH[:,2:end]),1) .> 0) & (sum(array(PresAbs_BR_NH_WH[:,2:end]),1) .> 0))]
names(PresAbs_BR_NH_WH2)

PresAbs_NB_NH_WH2 = PresAbs_NB_NH_WH[:,2:end][:, find((sum(array(PresAbs_NB_NH_WH[:,2:end]),1) .> 0) & (sum(array(PresAbs_BR_NH_WH[:,2:end]),1) .> 0))]
PresAbs_BR_SH_WH2 = PresAbs_BR_SH_WH[:,2:end][:, find((sum(array(PresAbs_NB_SH_WH[:,2:end]),1) .> 0) & (sum(array(PresAbs_BR_SH_WH[:,2:end]),1) .> 0))]
PresAbs_NB_SH_WH2 = PresAbs_NB_SH_WH[:,2:end][:, find((sum(array(PresAbs_NB_SH_WH[:,2:end]),1) .> 0) & (sum(array(PresAbs_BR_SH_WH[:,2:end]),1) .> 0))]

PresAbs_BR_NH_EH2 = PresAbs_BR_NH_EH[:,2:end][:, find((sum(array(PresAbs_NB_NH_EH[:,2:end]),1) .> 0) & (sum(array(PresAbs_BR_NH_EH[:,2:end]),1) .> 0))]
PresAbs_NB_NH_EH2 = PresAbs_NB_NH_EH[:,2:end][:, find((sum(array(PresAbs_NB_NH_EH[:,2:end]),1) .> 0) & (sum(array(PresAbs_BR_NH_EH[:,2:end]),1) .> 0))]
PresAbs_BR_SH_EH2 = PresAbs_BR_SH_EH[:,2:end][:, find((sum(array(PresAbs_NB_SH_EH[:,2:end]),1) .> 0) & (sum(array(PresAbs_BR_SH_EH[:,2:end]),1) .> 0))]
PresAbs_NB_SH_EH2 = PresAbs_NB_SH_EH[:,2:end][:, find((sum(array(PresAbs_NB_SH_EH[:,2:end]),1) .> 0) & (sum(array(PresAbs_BR_SH_EH[:,2:end]),1) .> 0))]

migr_WH_NH = array(PresAbs_BR_NH_WH2) - array(PresAbs_NB_NH_WH2)
migr_WH_SH = array(PresAbs_BR_SH_WH2) - array(PresAbs_NB_SH_WH2)
migr_EH_NH = array(PresAbs_BR_NH_EH2) - array(PresAbs_NB_NH_EH2)
migr_EH_SH = array(PresAbs_BR_SH_EH2) - array(PresAbs_NB_SH_EH2)

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

perm_WH_NH = array(PresAbs_BR_NH_WH2) .* array(PresAbs_NB_NH_WH2)
perm_WH_SH = array(PresAbs_BR_SH_WH2) .* array(PresAbs_NB_SH_WH2)
perm_EH_NH = array(PresAbs_BR_NH_EH2) .* array(PresAbs_NB_NH_EH2)
perm_EH_SH = array(PresAbs_BR_SH_EH2) .* array(PresAbs_NB_SH_EH2)

PresAbs_res_WH2 = array(PresAbs_res_WH[:,2:end][:, find(sum(array(PresAbs_res_WH[:,2:end]),1) .> 0)])
PresAbs_res_EH2 = array(PresAbs_res_EH[:,2:end][:, find(sum(array(PresAbs_res_EH[:,2:end]),1) .> 0)])

breeding_WH = sum(migr_br_WH_NH, 2) + sum(migr_br_WH_SH, 2)
nonbreeding_WH = sum(migr_nb_WH_NH, 2) + sum(migr_nb_WH_SH, 2)
residentsObs_WH = sum(PresAbs_res_WH2, 2) + sum(perm_WH_NH, 2) + sum(perm_WH_SH, 2)
prop_migr_obs_WH = (breeding_WH + nonbreeding_WH) ./ (breeding_WH + nonbreeding_WH + residentsObs_WH)
diff_obs_WH = breeding_WH - nonbreeding_WH

breeding_EH = sum(migr_br_EH_NH, 2) + sum(migr_br_EH_SH, 2)
nonbreeding_EH = sum(migr_nb_EH_NH, 2) + sum(migr_nb_EH_SH, 2)
residentsObs_EH = sum(PresAbs_res_EH2, 2) + sum(perm_EH_NH, 2) + sum(perm_EH_SH, 2)
prop_migr_obs_EH = (breeding_EH + nonbreeding_EH) ./ (breeding_EH + nonbreeding_EH + residentsObs_EH)
diff_obs_EH = breeding_EH - nonbreeding_EH


#### RANGE SIZE
rangeSize_migra_BR_NH_WH = sum(array(PresAbs_BR_NH_WH2), 1)
rangeSize_migra_BR_SH_WH = sum(array(PresAbs_BR_SH_WH2), 1)
rangeSize_migra_NB_NH_WH = sum(array(PresAbs_NB_NH_WH2), 1)
rangeSize_migra_NB_SH_WH = sum(array(PresAbs_NB_SH_WH2), 1)
rangeSize_migra_BR_NH_EH = sum(array(PresAbs_BR_NH_EH2), 1)
rangeSize_migra_BR_SH_EH = sum(array(PresAbs_BR_SH_EH2), 1)
rangeSize_migra_NB_NH_EH = sum(array(PresAbs_NB_NH_EH2), 1)
rangeSize_migra_NB_SH_EH = sum(array(PresAbs_NB_SH_EH2), 1)
rangeSize_res_WH = sum(PresAbs_res_WH2, 1)
rangeSize_res_EH = sum(PresAbs_res_EH2, 1)
rangeSizes_WH = hcat(rangeSize_migra_BR_NH_WH, rangeSize_migra_BR_SH_WH, rangeSize_migra_NB_NH_WH, rangeSize_migra_NB_SH_WH, rangeSize_res_WH)
rangeSizes_EH = hcat(rangeSize_migra_BR_NH_EH, rangeSize_migra_BR_SH_EH, rangeSize_migra_NB_NH_EH, rangeSize_migra_NB_SH_EH, rangeSize_res_EH)




### Simulate virtual ranges for a species

rangesize_sel = 131 #mean(rangeSizes_WH)
rsimul = zeros(rangesize_sel)

j=1
while length(unique(rsimul)) / 2323 < 1
	hex_1 = rand(1:1:2323)
	while length(findin(hex_1, rsimul)) != 0
		hex_1 = rand(1:1:2323)
  end
	rangeSimulated = copy(hex_1)

	i = 1
	nicheneig_sel = 0
	while length(rangeSimulated) < rangesize_sel # & findin(true,isna(nicheneig_sel)) == false
    neig = sum(pairwiseDistanceWH2[:,rangeSimulated], 2)
    neigs = [0]
    for i in 1:length(neig)
      neigs = vcat(neigs, rep([i], Int(neig[i])))
    end
    neigs = neigs[2:end]
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

i = i+1
zerange = zeros(2323)
zerange[rsimul[:,i]] = 1
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=zerange, style(default_point_size=1.3pt))

#rsimul2 = copy(rsimul)


### Parameters

mm = 0.0001016
rr = 0.3476501
tt = 13.3983
#t1 = 26.63075
#t2 = 0.02 # 0.03
beta2 = 1


###  Costs associated with migration, thermoregulation and reproduction

migrationDistance_inKm = zeros(length(rsimul[1,:]), length(rsimul[1,:]))
for i in 1:length(rsimul[1,:])
  for j in 1:length(rsimul[1,:])
    migrationDistance_inKm[i,j] = greatcircledistance(degtorad(mean(lonlatWH[rsimul[:,i],1])), degtorad(mean(lonlatWH[rsimul[:,i],2])), degtorad(mean(lonlatWH[rsimul[:,j],1])),degtorad(mean(lonlatWH[rsimul[:,j],2])))
  end
end
migrationCost = mm .* migrationDistance_inKm

thermoCost_NS = zeros(length(rsimul[1,:]), length(rsimul[1,:]))
thermoCost_NW = zeros(length(rsimul[1,:]), length(rsimul[1,:]))
for i in 1:length(rsimul[1,:])
  for j in 1:length(rsimul[1,:])
    thermoCost_NS[i,j] = (maximum([mean(tempNS_WH[rsimul[:,i]]), mean(tempNW_WH[rsimul[:,j]])]) - mean(tempNS_WH[rsimul[:,i]])) / tt
    thermoCost_NW[i,j] = (maximum([mean(tempNS_WH[rsimul[:,i]]), mean(tempNW_WH[rsimul[:,j]])]) - mean(tempNW_WH[rsimul[:,j]])) / tt
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


body_mass = complete_cases!(readtable("Body_size_birds.csv"))
M = median(body_mass[1])



### RUN THE MODEL

S_July_WH = zeros(length(lonlatWH[:,1]))
S_January_WH = zeros(length(lonlatWH[:,1]))
ranges_selected_BR_WH = zeros(2323,1)
ranges_selected_NB_WH = zeros(2323,1)
sel_WH = [0]
selected_ranges = [0]
energyAvailable_NS = copy(energySupplyJulyWH)
energyAvailable_NW = copy(energySupplyJanuaryWH)

k=1
while sum(energyAvailable_NS) > 10000 && sum(energyAvailable_NW) > 10000 #sum(S_January_WH) < sum(energySupplyJanuaryWH) && sum(S_July_WH) < sum(energySupplyJulyWH)

 supplyNS = zeros(length(rsimul[1,:]))
 for i in 1:length(rsimul[1,:])
   supplyNS[i] = mean(energyAvailable_NS[rsimul[:,i]])
   if supplyNS[i] == 0 supplyNS[i] = 0.1 end
 end
 supplyNW = zeros(length(rsimul[1,:]))
 for i in 1:length(rsimul[1,:])
   supplyNW[i] = mean(energyAvailable_NW[rsimul[:,i]])
   if supplyNW[i] == 0 supplyNW[i] = 0.1 end
 end

 state_proba_brNS = zeros(length(rsimul[1,:]), length(rsimul[1,:]))
 state_proba_brNW = zeros(length(rsimul[1,:]), length(rsimul[1,:]))
 for i in 1:length(rsimul[1,:])
   for j in 1:length(rsimul[1,:])
     state_proba_brNS[i,j] = exp(-beta2 .* (migration_thermo_cost_brNS[i,j] ./ supplyNS[i]))
     state_proba_brNW[i,j] = exp(-beta2 .* (migration_thermo_cost_brNW[i,j] ./ supplyNW[j]))
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
 ra = zeros(2323)
 ra[rsimul[:,rangeBR]] = 1
 ranges_selected_BR_WH = hcat(ranges_selected_BR_WH, ra)
 ra = zeros(2323)
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
plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=PropMigrSimu, style(default_point_size=1.3pt))
maximum(richnessDiffSimu)

plot(x=lonlatWH[:,2], y=energyAvailable_NS)
plot(x=lonlatWH[:,2], y=energyAvailable_NW)

plot(x=lonlatWH[:,2], y=S_July_WH)
plot(x=lonlatWH[:,2], y=S_January_WH)

plot(x=lonlatWH[:,2], y=richnessMigraBRsimu)
plot(x=lonlatWH[:,2], y=richnessMigraNBsimu)
plot(x=lonlatWH[:,2], y=richnessDiffSimu)
plot(x=lonlatWH[:,2], y=PropMigrSimu)

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
sum(numberMigraNBsimu)






## Proba model

state_proba_brNS_pr = state_proba_brNS ./ (sum(state_proba_brNS) + sum(state_proba_brNW))
state_proba_brNW_pr = state_proba_brNW ./ (sum(state_proba_brNS) + sum(state_proba_brNW))

selected_link = wsample(1:(length(state_proba_brNS_pr[:]) + length(state_proba_brNW_pr[:])), vcat(state_proba_brNS_pr[:],state_proba_brNW_pr[:]), 1)[1]
if selected_link <= length(state_proba_brNS_pr[:])
  rangeNB = Int(ceil(selected_link/length(rsimul[1,:]))[1])
  rangeBR = selected_link%length(rsimul[1,:])
  if rangeBR == 0 rangeBR = Int(selected_link/length(rsimul[1,:])) end
  ra = zeros(2323)
  ra[rsimul[:,rangeBR]] = 1
  ranges_selected_BR_WH = hcat(ranges_selected_BR_WH, ra)
  ra = zeros(2323)
  ra[rsimul[:,rangeNB]] = 1
  ranges_selected_NB_WH = hcat(ranges_selected_NB_WH, ra)
  S_January_WH[ranges_selected_NB_WH[:,end] .== 1] = S_January_WH[ranges_selected_NB_WH[:,end] .== 1] .+ migration_thermo_cost_brNS_NW[rangeBR, rangeNB][1]
  S_July_WH[ranges_selected_BR_WH[:,end] .== 1] = S_July_WH[ranges_selected_BR_WH[:,end] .== 1] .+ migration_thermo_cost_brNS_NS[rangeBR, rangeNB][1]
else
  selected_link = selected_link - length(state_proba_brNS_pr[:])
  rangeBR = Int(ceil(selected_link/length(rsimul[1,:]))[1])
  rangeNB = selected_link%length(rsimul[1,:])
  if rangeNB == 0 rangeNB = Int(selected_link/length(rsimul[1,:])) end
  ra = zeros(2323)
  ra[rsimul[:,rangeBR]] = 1
  ranges_selected_BR_WH = hcat(ranges_selected_BR_WH, ra)
  ra = zeros(2323)
  ra[rsimul[:,rangeNB]] = 1
  ranges_selected_NB_WH = hcat(ranges_selected_NB_WH, ra)
  S_January_WH[ranges_selected_BR_WH[:,end] .== 1] = S_January_WH[ranges_selected_BR_WH[:,end] .== 1] .+ migration_thermo_cost_brNW_NW[rangeNB, rangeBR][1]
  S_July_WH[ranges_selected_NB_WH[:,end] .== 1] = S_July_WH[ranges_selected_NB_WH[:,end] .== 1] .+ migration_thermo_cost_brNW_NS[rangeNB, rangeBR][1]
end







LCT = zeros(length(rsimul[1,:]), length(rsimul[1,:]))
thermoCost_NS = zeros(length(rsimul[1,:]), length(rsimul[1,:]))
thermoCost_NW = zeros(length(rsimul[1,:]), length(rsimul[1,:]))
for i in 1:length(rsimul[1,:])
  for j in 1:length(rsimul[1,:])
    if maximum([mean(tempNS_WH[rsimul[:,i]]), mean(tempNW_WH[rsimul[:,j]])]) < t1 LCT[i,j] = maximum([mean(tempNS_WH[rsimul[:,i]]), mean(tempNW_WH[rsimul[:,j]])]) else LCT[i,j] = t1 end
    thermoCost_NS[i,j] = LCT[i,j] - mean(tempNS_WH[rsimul[:,i]])
    if thermoCost_NS[i,j] < 0 thermoCost_NS[i,j] = 0 end
    thermoCost_NS[i,j] = t2 * thermoCost_NS[i,j]
    thermoCost_NW[i,j] = LCT[i,j] - mean(tempNW_WH[rsimul[:,j]])
    if thermoCost_NW[i,j] < 0 thermoCost_NW[i,j] = 0 end
    thermoCost_NW[i,j] = t2 * thermoCost_NW[i,j]
  end
end
