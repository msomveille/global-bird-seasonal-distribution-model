using DataFrames
using HDF5
using JLD

function migraModel(n, mm, rr, tt, ee)

  #cd("/data/zool-avian-social-ecology/zool2068")

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

  rsimulWH = Array(readtable("simulatedRangesWH2.csv", header = false))
  rsimulEH = Array(readtable("simulatedRangesEH2.csv", header = false))
  dataloadedWH = Array(readtable("dataloadedWH.csv"))
  dataloadedEH = Array(readtable("dataloadedEH.csv"))

  lonlatWH = dataloadedWH[:,[1,2]]
  lonlatEH = dataloadedEH[:,[1,2]]
  tempNS_WH = dataloadedWH[:,3]
  tempNW_WH = dataloadedWH[:,4]
  tempNS_EH = dataloadedEH[:,3]
  tempNW_EH = dataloadedEH[:,4]
  NDVI_ns_WH = dataloadedWH[:,5]
  NDVI_nw_WH = dataloadedWH[:,6]
  NDVI_ns_EH = dataloadedEH[:,5]
  NDVI_nw_EH = dataloadedEH[:,6]


  resmec = (zeros(7352), zeros(4))

  for jj in 1:n

  ### Parameters

  #ee = 8.89 #sample(collect(0:0.000001:0.0002), 1)[1] # 8.89
  #mm = 1.5063e-7 * (259^-0.1082) * 1000 #sample(collect(0:0.000001:0.0002), 1)[1] # 1.5063e-7 * (259^-0.1082) * 1000
  #rr = 0.725 * (259^(-0.0872-(0.0727*log(10, 259)))) #sample(collect(0:0.01:1), 1)[1] # 0.725 * (259^(-0.0872-(0.0727*log(10, 259))))
  #tt = 8.1403 * (259^0.1364) #sample(collect(0:0.1:40), 1)[1] # 8.1403 * (259^0.1364)

  param = vcat(mm, rr, tt, ee)

  # Energy supply
  energySupplyJulyWH = (ee * NDVI_ns_WH).^2
  energySupplyJanuaryWH = (ee * NDVI_nw_WH).^2
  energySupplyJulyEH = (ee * NDVI_ns_EH).^2
  energySupplyJanuaryEH = (ee * NDVI_nw_EH).^2


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
  ranges_selected_BR_WH = zeros(2323,1)
  ranges_selected_NB_WH = zeros(2323,1)
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
    ra = zeros(2323)
    ra[rsimulWH[:,rangeBR]] = 1
    ranges_selected_BR_WH = hcat(ranges_selected_BR_WH, ra)
    ra = zeros(2323)
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

  resmecWH = hcat(richnessMigraBRsimuWH, richnessMigraNBsimuWH, richnessResidentsimuWH)

  #plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessBRsimuWH, style(default_point_size=1.3pt))
  #plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessNBsimuWH, style(default_point_size=1.3pt))
  #plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessMigraBRsimuWH, style(default_point_size=1.3pt))
  #plot(x=lonlatWH[:,1], y=lonlatWH[:,2], color=richnessMigraNBsimuWH, style(default_point_size=1.3pt))
  #plot(x=lonlatWH[:,2], y=richnessMigraBRsimuWH, style(default_point_size=1.3pt))
  #plot(x=lonlatWH[:,2], y=richnessMigraNBsimuWH, style(default_point_size=1.3pt))
  #plot(x=resmecha[:,2], y=resmecha[:,5], style(default_point_size=1.3pt))
  #plot(x=resmecha[:,2], y=resmecha[:,7], style(default_point_size=1.3pt))



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
  ranges_selected_BR_EH = zeros(5029,1)
  ranges_selected_NB_EH = zeros(5029,1)
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

    ra = zeros(5029)
    ra[rsimulEH[:,rangeBR]] = 1
    ranges_selected_BR_EH = hcat(ranges_selected_BR_EH, ra)
    ra = zeros(5029)
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

  resmecEH = hcat(richnessMigraBRsimuEH, richnessMigraNBsimuEH, richnessResidentsimuEH)

  resmec = (hcat(resmec[1], vcat(resmecWH, resmecEH)), hcat(resmec[2], param))

  mm = mm + ((1.5063e-7 * (259^-0.1082) * 1000 * 2)/9)

  end
  #resmec = resmec[2:end]
  return(resmec)
end


#resmech = migraModel(2, 0, 0, 1, 0.03)


#resmech1 = resmech[1][:,2:end]
#resmech2 = resmech[2][:,2:end]

#abc = sample(collect(1:1:10000000),1)[1]
#pathh = string("resmec", abc, ".jld")
#save(pathh, "data", resmech)

#resmec1 = load("resmec7285027.jld")
#resmec1["data"][1]

#resmecha = hcat(vcat(lonlatWH, lonlatEH), resmech1)
#plot(x=resmecha[:,1], y=resmecha[:,2], color=resmecha[:,3], style(default_point_size=1.3pt))
#plot(x=resmecha[:,1], y=resmecha[:,2], color=resmecha[:,4], style(default_point_size=1.3pt))
#plot(x=resmecha[:,1], y=resmecha[:,2], color=resmecha[:,5], style(default_point_size=1.3pt))
#plot(x=resmecha[:,2], y=resmecha[:,3], style(default_point_size=1.3pt))
#plot(x=resmecha[:,2], y=resmecha[:,4], style(default_point_size=1.3pt))
#plot(x=resmecha[:,2], y=resmecha[:,5], style(default_point_size=1.3pt))
#plot(x=resmecha[:,2], y=resmecha[:,7], style(default_point_size=1.3pt))
