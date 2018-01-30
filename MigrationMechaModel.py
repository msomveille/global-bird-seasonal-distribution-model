from sys import argv
import csv
import numpy as np
from pyemd import emd
from pyemd import emd_with_flow
from math import radians, cos, sin, asin, sqrt
import matplotlib.pyplot as plt
import statsmodels.api as sm
import copy
from scipy.stats.stats import pearsonr 
import os
import random

os.chdir('/Users/mariussomveille/Desktop/mecha_cluster/PythonModel')
#os.chdir('/data/zool-avian-social-ecology/zool2068/PythonModel')

#Observed presences-absences of species
PresAbs_BR_NH_WH = np.loadtxt('PresAbs_BR_NH_WH.csv', skiprows=1, usecols=range(0,468), delimiter=',')
PresAbs_NB_NH_WH = np.loadtxt('PresAbs_NB_NH_WH.csv', skiprows=1, usecols=range(0,468), delimiter=',')
PresAbs_BR_SH_WH = np.loadtxt('PresAbs_BR_SH_WH.csv', skiprows=1, usecols=range(0,158), delimiter=',')
PresAbs_NB_SH_WH = np.loadtxt('PresAbs_NB_SH_WH.csv', skiprows=1, usecols=range(0,158), delimiter=',')
PresAbs_res_WH = np.loadtxt('PresAbs_res_WH.csv', skiprows=1, usecols=range(0,3418), delimiter=',')
PresAbs_BR_NH_EH = np.loadtxt('PresAbs_BR_NH_EH.csv', skiprows=1, usecols=range(0,798), delimiter=',')
PresAbs_NB_NH_EH = np.loadtxt('PresAbs_NB_NH_EH.csv', skiprows=1, usecols=range(0,798), delimiter=',')
PresAbs_BR_SH_EH = np.loadtxt('PresAbs_BR_SH_EH.csv', skiprows=1, usecols=range(0,106), delimiter=',')
PresAbs_NB_SH_EH = np.loadtxt('PresAbs_NB_SH_EH.csv', skiprows=1, usecols=range(0,106), delimiter=',')
PresAbs_res_EH = np.loadtxt('PresAbs_res_EH.csv', skiprows=1, usecols=range(0,4408), delimiter=',')


#Environmental data  ---  NS:northern summer -- NW:northern winter -- WH:western Hemisphere -- EH:eastern Hemisphere
dataloadedWH = np.loadtxt('dataloadedWH.csv', skiprows=1, delimiter=',')
dataloadedEH = np.loadtxt('dataloadedEH.csv', skiprows=1, delimiter=',') 
tempNS_WH = dataloadedWH[:,2]
tempNW_WH = dataloadedWH[:,3]
tempNS_EH = dataloadedEH[:,2]
tempNW_EH = dataloadedEH[:,3]
NDVI_NS_WH = dataloadedWH[:,4]
NDVI_NW_WH = dataloadedWH[:,5]
NDVI_NS_EH = dataloadedEH[:,4]
NDVI_NW_EH = dataloadedEH[:,5]

lonWH = dataloadedWH[:,0]
latWH = dataloadedWH[:,1]
lonEH = dataloadedEH[:,0]
latEH = dataloadedEH[:,1]



#Compute distance matrices for the calculation of the Earth Mover's Distance
def haversine(lon1, lat1, lon2, lat2):
   """
   Calculate the great circle distance between two points 
   on the earth (specified in decimal degrees)
   """
   # convert decimal degrees to radians 
   lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
   # haversine formula 
   dlon = lon2 - lon1 
   dlat = lat2 - lat1 
   a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
   c = 2 * asin(sqrt(a)) 
   km = 6367 * c
   return km

T = 2000  ## threshold for ground distance
#distance_matrix_WH = np.zeros(shape=(len(lonWH),len(lonWH)))
#for i in range(0,len(lonWH)):
#   for j in range(0,len(lonWH)):
#       distance_matrix_WH[i,j] = haversine(lonWH[i], latWH[i], lonWH[j], latWH[j])
#       if distance_matrix_WH[i,j] > T:
#           distance_matrix_WH[i,j] = T

#distance_matrix_EH = np.zeros(shape=(len(lonEH),len(lonEH)))
#for i in range(0,len(lonEH)):
#   for j in range(0,len(lonEH)):
#       distance_matrix_EH[i,j] = haversine(lonEH[i], latEH[i], lonEH[j], latEH[j])
#       if distance_matrix_EH[i,j] > T:
#           distance_matrix_EH[i,j] = T
lon = np.concatenate((lonWH,lonEH))
lat = np.concatenate((latWH,latEH))
distance_matrix = np.zeros(shape=(len(lon),len(lon)))
for i in range(0,len(lon)):
   for j in range(0,len(lon)):
       distance_matrix[i,j] = haversine(lon[i], lat[i], lon[j], lat[j])
       if distance_matrix[i,j] > T:
           distance_matrix[i,j] = T



#Observed global spatial patterns
migr_WH_NH = PresAbs_BR_NH_WH - PresAbs_NB_NH_WH
migr_WH_SH = PresAbs_BR_SH_WH - PresAbs_NB_SH_WH
migr_EH_NH = PresAbs_BR_NH_EH - PresAbs_NB_NH_EH
migr_EH_SH = PresAbs_BR_SH_EH - PresAbs_NB_SH_EH
migr_br_WH_NH = copy.copy(migr_WH_NH)
migr_br_WH_NH[(migr_br_WH_NH == -1)] = 0
migr_nb_WH_NH = copy.copy(migr_WH_NH)
migr_nb_WH_NH[(migr_nb_WH_NH == 1)] = 0
migr_nb_WH_NH = np.abs(migr_nb_WH_NH)
migr_br_WH_SH = copy.copy(migr_WH_SH)
migr_br_WH_SH[(migr_br_WH_SH == -1)] = 0
migr_nb_WH_SH = copy.copy(migr_WH_SH)
migr_nb_WH_SH[(migr_nb_WH_SH == 1)] = 0
migr_nb_WH_SH = np.abs(migr_nb_WH_SH)
migr_br_EH_NH = copy.copy(migr_EH_NH)
migr_br_EH_NH[(migr_br_EH_NH == -1)] = 0
migr_nb_EH_NH = copy.copy(migr_EH_NH)
migr_nb_EH_NH[(migr_nb_EH_NH == 1)] = 0
migr_nb_EH_NH = np.abs(migr_nb_EH_NH)
migr_br_EH_SH = copy.copy(migr_EH_SH)
migr_br_EH_SH[(migr_br_EH_SH == -1)] = 0
migr_nb_EH_SH = copy.copy(migr_EH_SH)
migr_nb_EH_SH[(migr_nb_EH_SH == 1)] = 0
migr_nb_EH_SH = np.abs(migr_nb_EH_SH)
perm_WH_NH = PresAbs_BR_NH_WH * PresAbs_NB_NH_WH
perm_WH_SH = PresAbs_BR_SH_WH * PresAbs_NB_SH_WH
perm_EH_NH = PresAbs_BR_NH_EH * PresAbs_NB_NH_EH
perm_EH_SH = PresAbs_BR_SH_EH * PresAbs_NB_SH_EH

breeding_WH = migr_br_WH_NH.sum(axis=1) + migr_br_WH_SH.sum(axis=1)
nonbreeding_WH = migr_nb_WH_NH.sum(axis=1) + migr_nb_WH_SH.sum(axis=1)
residentsObs_WH = PresAbs_res_WH.sum(axis=1) + perm_WH_NH.sum(axis=1) + perm_WH_SH.sum(axis=1)
prop_migr_obs_WH = (breeding_WH + nonbreeding_WH) / (breeding_WH + nonbreeding_WH + residentsObs_WH)
diff_obs_WH = breeding_WH - nonbreeding_WH

breeding_EH = migr_br_EH_NH.sum(axis=1) + migr_br_EH_SH.sum(axis=1)
nonbreeding_EH = migr_nb_EH_NH.sum(axis=1) + migr_nb_EH_SH.sum(axis=1)
residentsObs_EH = PresAbs_res_EH.sum(axis=1) + perm_EH_NH.sum(axis=1) + perm_EH_SH.sum(axis=1)
prop_migr_obs_EH = (breeding_EH + nonbreeding_EH) / (breeding_EH + nonbreeding_EH + residentsObs_EH)
diff_obs_EH = breeding_EH - nonbreeding_EH

 
#Import simulated ranges (400 ranges randomly placed across the eastern Hemisphere and 600 ranges randomly placed across the western Hemisphere)
rsimulWH = np.loadtxt('simulatedRangesWH2.csv', delimiter=',')
rsimulEH = np.loadtxt('simulatedRangesEH2.csv', delimiter=',')



#Parameters
mm = float(argv[1]) #3.69E-05 #6.86E-05 #float(argv[1]) #0.0000645 #related to the cost of migration
rr = float(argv[2]) #0.810904615 #0.409030967 #float(argv[2]) #0.195 #related to the cost of reproduction
tt = float(argv[3]) #24.57998175 #23.56825437 #float(argv[3]) #13.4 #related to the cost of thermoregulation
ee = float(argv[4]) #76.01626261 #99.40440512 #float(argv[4]) #81.36 #to scale the energy supply from NDVI data

#output_data_file = open(os.path.join("MechaResults"+str(mm)+"_"+str(rr)+"_"+str(tt)+"_"+str(ee)+".csv"),'w')
#f = open(argv[5], 'w')

#Energy supply
energySupply_NS_WH = NDVI_NS_WH * ee
energySupply_NW_WH = NDVI_NW_WH * ee
energySupply_NS_EH = NDVI_NS_EH * ee
energySupply_NW_EH = NDVI_NW_EH * ee


### Run the model in Western Hemisphere ### --- To be repeated for the eastern Hemisphere

#Costs associated with migration, thermoregulation and reproduction
migrationDistance_inKm = np.zeros(shape=(len(rsimulWH[0,:]), len(rsimulWH[0,:])))
for i in range(0,len(rsimulWH[0,:])):
    for j in range(0,len(rsimulWH[0,:])):
        migrationDistance_inKm[i,j] = haversine(np.mean(lonWH[map(int,(rsimulWH[:,i]-1))]), np.mean(latWH[map(int,(rsimulWH[:,i]-1))]), np.mean(lonWH[map(int,(rsimulWH[:,j]-1))]), np.mean(latWH[map(int,(rsimulWH[:,j]-1))]))
migrationCost = mm * migrationDistance_inKm

thermoCost_NS = np.zeros(shape=(len(rsimulWH[0,:]), len(rsimulWH[0,:])))
thermoCost_NW = np.zeros(shape=(len(rsimulWH[0,:]), len(rsimulWH[0,:])))
for i in range(0,len(rsimulWH[0,:])):
    for j in range(0,len(rsimulWH[0,:])):
        thermoCost_NS[i,j] = (40 - tt - np.mean(tempNS_WH[map(int,(rsimulWH[:,i]-1))])) / tt
        thermoCost_NW[i,j] = (40 - tt - np.mean(tempNW_WH[map(int,(rsimulWH[:,j]-1))])) / tt
        if thermoCost_NS[i,j] < 0:
            thermoCost_NS[i,j] = 0
        if thermoCost_NW[i,j] < 0:
            thermoCost_NW[i,j] = 0

migration_thermo_cost_brNS_NS = 1 + migrationCost + thermoCost_NS + rr  #This is the total energetic requirements during the northern summer if the species breed during the northern summer  --- brNS:breed during the northern summer
migration_thermo_cost_brNS_NW = 1 + migrationCost + thermoCost_NW       
migration_thermo_cost_brNW_NS = 1 + migrationCost + thermoCost_NS       #This is the total energetic requirements during the northern summer if the species breed during the northern winter  --- brNW:breed during the northern winter  
migration_thermo_cost_brNW_NW = 1 + migrationCost + thermoCost_NW + rr  

#Declare objects
energyUsed_NS_WH = np.zeros(shape=len(lonWH))  #cumulative energy used by simulated species
energyUsed_NW_WH = np.zeros(shape=len(lonWH))
ranges_selected_BR_WH = np.zeros(shape=len(lonWH))
ranges_selected_NB_WH = np.zeros(shape=len(lonWH))
energyAvailable_NS_WH = copy.copy(energySupply_NS_WH)  #energy available as a function of energy supply and energy used by previously simulated species
energyAvailable_NW_WH = copy.copy(energySupply_NW_WH)
#rangesLatWH = np.zeros(shape=(10000, 2))  #to track the latitude of the selected ranges' centroid
#breedingSeasonWH = np.zeros(shape=10000)  #to track the breeding season selected

#Run the model
k=0
while (sum(energyAvailable_NS_WH) > (sum(energySupply_NS_WH) * 0.05)) & (sum(energyAvailable_NW_WH) > (sum(energySupply_NW_WH) * 0.05)):

  #Average energy available in each range
  supplyNS = np.zeros(len(rsimulWH[0,:]))
  for i in range(0,len(rsimulWH[0,:])):
      supplyNS[i] = np.mean(energyAvailable_NS_WH[map(int,(rsimulWH[:,i]-1))])
      if supplyNS[i] == 0:
          supplyNS[i] = 0.1
  supplyNW = np.zeros(len(rsimulWH[0,:]))
  for i in range(0,len(rsimulWH[0,:])):
      supplyNW[i] = np.mean(energyAvailable_NW_WH[map(int,(rsimulWH[:,i]-1))])
      if supplyNW[i] == 0:
          supplyNW[i] = 0.1

  #Compute an energy-efficiency score (E) associated with each distribution option (i.e. combination of breeding and non-breeding ranges for a given breeding season)      
  E_brNS = np.zeros(shape=(len(rsimulWH[0,:]), len(rsimulWH[0,:])))
  E_brNW = np.zeros(shape=(len(rsimulWH[0,:]), len(rsimulWH[0,:])))
  for i in range(0,len(rsimulWH[0,:])):
      for j in range(0,len(rsimulWH[0,:])):
          E_brNS[i,j] = (migration_thermo_cost_brNS_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNS_NW[i,j] / supplyNW[j])
          E_brNW[i,j] = (migration_thermo_cost_brNW_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNW_NW[i,j] / supplyNW[j])

  #Find the most energy efficient distributional strategy (i.e. the lowest value across both matrices) and update the energy available
  best_link1 = np.where(E_brNS == E_brNS.min()) 
  best_link2 = np.where(E_brNW == E_brNW.min()) 
  if E_brNS.min() < E_brNW.min():
      if len(best_link1[1]) > 1:
          sel = random.sample(range(0,len(best_link1[1])), 1)
          rangeBR = best_link1[0][sel]
          rangeNB = best_link1[1][sel]
      else:
          rangeBR = best_link1[0]
          rangeNB = best_link1[1]          
      energyUsed_NW_WH[map(int,rsimulWH[:,rangeNB]-1)] = energyUsed_NW_WH[map(int,rsimulWH[:,rangeNB]-1)] + migration_thermo_cost_brNS_NW[rangeBR, rangeNB]
      energyUsed_NS_WH[map(int,rsimulWH[:,rangeBR]-1)] = energyUsed_NS_WH[map(int,rsimulWH[:,rangeBR]-1)] + migration_thermo_cost_brNS_NS[rangeBR, rangeNB]      
      #breedingSeasonWH[k] = 0
  else:
      if len(best_link2[1]) > 1:
          sel = random.sample(range(0,len(best_link2[1])), 1)
          rangeBR = best_link2[1][sel]
          rangeNB = best_link2[0][sel]
      else:
         rangeBR = best_link2[1]
         rangeNB = best_link2[0]
      energyUsed_NW_WH[map(int,rsimulWH[:,rangeBR]-1)] = energyUsed_NW_WH[map(int,rsimulWH[:,rangeBR]-1)] + migration_thermo_cost_brNW_NW[rangeNB, rangeBR]
      energyUsed_NS_WH[map(int,rsimulWH[:,rangeNB]-1)] = energyUsed_NS_WH[map(int,rsimulWH[:,rangeNB]-1)] + migration_thermo_cost_brNW_NS[rangeNB, rangeBR]
      #breedingSeasonWH[k] = 1
  energyAvailable_NS_WH = energySupply_NS_WH - energyUsed_NS_WH
  energyAvailable_NW_WH = energySupply_NW_WH - energyUsed_NW_WH
  for i in range(0,len(energyAvailable_NS_WH)):
      if energyAvailable_NS_WH[i] < 0:
          energyAvailable_NS_WH[i] = 0
      if energyAvailable_NW_WH[i] < 0:
          energyAvailable_NW_WH[i] = 0

  #rangesLatWH[k,0] = np.mean(latWH[map(int, rsimulWH[:,rangeBR]-1)])
  #rangesLatWH[k,1] = np.mean(latWH[map(int, rsimulWH[:,rangeNB]-1)])

  #Update the newly occupied hexagons
  ra = np.zeros(shape=len(lonWH))
  ra[map(int,rsimulWH[:,rangeBR]-1)] = 1
  ranges_selected_BR_WH = np.column_stack((ranges_selected_BR_WH, ra))
  ra = np.zeros(shape=len(lonWH))
  ra[map(int,rsimulWH[:,rangeNB]-1)] = 1
  ranges_selected_NB_WH = np.column_stack((ranges_selected_NB_WH, ra))

  k=k+1

#rangesLatWH = rangesLatWH[0:k,:]
#breedingSeasonWH = breedingSeasonWH[0:k]


#Compute the simulated spatial patterns
ranges_selected_BR_WH = ranges_selected_BR_WH[:,1:len(ranges_selected_BR_WH[0,:])]
ranges_selected_NB_WH = ranges_selected_NB_WH[:,1:len(ranges_selected_NB_WH[0,:])]
migra = ranges_selected_BR_WH - ranges_selected_NB_WH
res = ranges_selected_BR_WH + ranges_selected_NB_WH
migraBR = np.zeros(shape=(len(migra[:,0]), len(migra[0,:])))
migraNB = np.zeros(shape=(len(migra[:,0]), len(migra[0,:])))
residents = np.zeros(shape=(len(migra[:,0]), len(migra[0,:])))
for i in range(0,len(ranges_selected_BR_WH[:,0])):
  for j in range(0,len(ranges_selected_BR_WH[0,:])):
      if migra[i,j] > 0:
          migraBR[i,j] = 1
      if migra[i,j] < 0:
          migraNB[i,j] = 1
      if res[i,j] == 2:
          residents[i,j] = 1
richnessBRsimuWH = np.zeros(shape=len(ranges_selected_BR_WH[:,0]))
richnessNBsimuWH = np.zeros(shape=len(ranges_selected_NB_WH[:,0]))
richnessMigraBRsimuWH = np.zeros(shape=len(ranges_selected_BR_WH[:,0]))
richnessMigraNBsimuWH = np.zeros(shape=len(ranges_selected_NB_WH[:,0]))
richnessResidentsimuWH = np.zeros(shape=len(ranges_selected_NB_WH[:,0]))
for i in range(0,len(ranges_selected_BR_WH[:,0])):
    richnessBRsimuWH[i] = sum(ranges_selected_BR_WH[i,:])
    richnessNBsimuWH[i] = sum(ranges_selected_NB_WH[i,:])
    richnessMigraBRsimuWH[i] = sum(migraBR[i,:])
    richnessMigraNBsimuWH[i] = sum(migraNB[i,:])
    richnessResidentsimuWH[i] = sum(residents[i,:])
richnessDiffSimuWH = richnessBRsimuWH - richnessNBsimuWH
PropMigrSimuWH = (richnessMigraBRsimuWH + richnessMigraNBsimuWH) / (richnessMigraBRsimuWH + richnessMigraNBsimuWH + richnessResidentsimuWH)

#BRpop_WH = np.zeros(shape=len(ranges_selected_NB_WH[0,:]))
#NBpop_WH = np.zeros(shape=len(ranges_selected_NB_WH[0,:]))
#respop_WH = np.zeros(shape=len(ranges_selected_NB_WH[0,:]))
#for i in range(0,len(ranges_selected_BR_WH[0,:])):
#    BRpop_WH[i] = sum(migraBR[:,i])
#    NBpop_WH[i] = sum(migraNB[:,i])
#    respop_WH[i] = 180 - BRpop_WH[i]
#BRpops_WH = sum(BRpop_WH)
#NBpops_WH = sum(NBpop_WH)
#respops_WH = sum(respop_WH)




### Run the model in Eastern Hemisphere ###

#Costs associated with migration, thermoregulation and reproduction
migrationDistance_inKm = np.zeros(shape=(len(rsimulEH[0,:]), len(rsimulEH[0,:])))
for i in range(0,len(rsimulEH[0,:])):
    for j in range(0,len(rsimulEH[0,:])):
        migrationDistance_inKm[i,j] = haversine(np.mean(lonEH[map(int,(rsimulEH[:,i]-1))]), np.mean(latEH[map(int,(rsimulEH[:,i]-1))]), np.mean(lonEH[map(int,(rsimulEH[:,j]-1))]), np.mean(latEH[map(int,(rsimulEH[:,j]-1))]))
migrationCost = mm * migrationDistance_inKm

thermoCost_NS = np.zeros(shape=(len(rsimulEH[0,:]), len(rsimulEH[0,:])))
thermoCost_NW = np.zeros(shape=(len(rsimulEH[0,:]), len(rsimulEH[0,:])))
for i in range(0,len(rsimulEH[0,:])):
    for j in range(0,len(rsimulEH[0,:])):
        thermoCost_NS[i,j] = (40 - tt - np.mean(tempNS_EH[map(int,(rsimulEH[:,i]-1))])) / tt
        thermoCost_NW[i,j] = (40 - tt - np.mean(tempNW_EH[map(int,(rsimulEH[:,j]-1))])) / tt
        if thermoCost_NS[i,j] < 0:
            thermoCost_NS[i,j] = 0
        if thermoCost_NW[i,j] < 0:
            thermoCost_NW[i,j] = 0

migration_thermo_cost_brNS_NS = 1 + migrationCost + thermoCost_NS + rr  #This is the total energetic requirements during the northern summer if the species breed during the northern summer  --- brNS:breed during the northern summer
migration_thermo_cost_brNS_NW = 1 + migrationCost + thermoCost_NW       
migration_thermo_cost_brNW_NS = 1 + migrationCost + thermoCost_NS       #This is the total energetic requirements during the northern summer if the species breed during the northern winter  --- brNW:breed during the northern winter  
migration_thermo_cost_brNW_NW = 1 + migrationCost + thermoCost_NW + rr  

#Declare objects
energyUsed_NS_EH = np.zeros(shape=len(lonEH))  #cumulative energy used by simulated species
energyUsed_NW_EH = np.zeros(shape=len(lonEH))
ranges_selected_BR_EH = np.zeros(shape=len(lonEH))
ranges_selected_NB_EH = np.zeros(shape=len(lonEH))
energyAvailable_NS_EH = copy.copy(energySupply_NS_EH)  #energy available as a function of energy supply and energy used by previously simulated species
energyAvailable_NW_EH = copy.copy(energySupply_NW_EH)
#rangesLatEH = np.zeros(shape=(10000, 2))  #to track the latitude of the selected ranges' centroid
#breedingSeasonEH = np.zeros(shape=10000)  #to track the breeding season selected

#Run the model
k=0
while (sum(energyAvailable_NS_EH) > (sum(energySupply_NS_EH) * 0.05)) & (sum(energyAvailable_NW_EH) > (sum(energySupply_NW_EH) * 0.05)):

  #Average energy available in each range
  supplyNS = np.zeros(len(rsimulEH[0,:]))
  for i in range(0,len(rsimulEH[0,:])):
      supplyNS[i] = np.mean(energyAvailable_NS_EH[map(int,(rsimulEH[:,i]-1))])
      if supplyNS[i] == 0:
          supplyNS[i] = 0.1
  supplyNW = np.zeros(len(rsimulEH[0,:]))
  for i in range(0,len(rsimulEH[0,:])):
      supplyNW[i] = np.mean(energyAvailable_NW_EH[map(int,(rsimulEH[:,i]-1))])
      if supplyNW[i] == 0:
          supplyNW[i] = 0.1

  #Compute an energy-efficiency score (E) associated with each distribution option (i.e. combination of breeding and non-breeding ranges for a given breeding season)      
  E_brNS = np.zeros(shape=(len(rsimulEH[0,:]), len(rsimulEH[0,:])))
  E_brNW = np.zeros(shape=(len(rsimulEH[0,:]), len(rsimulEH[0,:])))
  for i in range(0,len(rsimulEH[0,:])):
      for j in range(0,len(rsimulEH[0,:])):
          E_brNS[i,j] = (migration_thermo_cost_brNS_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNS_NW[i,j] / supplyNW[j])
          E_brNW[i,j] = (migration_thermo_cost_brNW_NS[i,j] / supplyNS[i]) + (migration_thermo_cost_brNW_NW[i,j] / supplyNW[j])

  #Find the most energy efficient distributional strategy (i.e. the lowest value across both matrices) and update the energy available
  best_link1 = np.where(E_brNS == E_brNS.min()) 
  best_link2 = np.where(E_brNW == E_brNW.min()) 
  if E_brNS.min() < E_brNW.min():
      if len(best_link1[1]) > 1:
          sel = random.sample(range(0,len(best_link1[1])), 1)
          rangeBR = best_link1[0][sel]
          rangeNB = best_link1[1][sel]
      else:
          rangeBR = best_link1[0]
          rangeNB = best_link1[1]          
      energyUsed_NW_EH[map(int,rsimulEH[:,rangeNB]-1)] = energyUsed_NW_EH[map(int,rsimulEH[:,rangeNB]-1)] + migration_thermo_cost_brNS_NW[rangeBR, rangeNB]
      energyUsed_NS_EH[map(int,rsimulEH[:,rangeBR]-1)] = energyUsed_NS_EH[map(int,rsimulEH[:,rangeBR]-1)] + migration_thermo_cost_brNS_NS[rangeBR, rangeNB]      
      #breedingSeasonEH[k] = 0
  else:
      if len(best_link2[1]) > 1:
          sel = random.sample(range(0,len(best_link2[1])), 1)
          rangeBR = best_link2[1][sel]
          rangeNB = best_link2[0][sel]
      else:
         rangeBR = best_link2[1]
         rangeNB = best_link2[0]
      energyUsed_NW_EH[map(int,rsimulEH[:,rangeBR]-1)] = energyUsed_NW_EH[map(int,rsimulEH[:,rangeBR]-1)] + migration_thermo_cost_brNW_NW[rangeNB, rangeBR]
      energyUsed_NS_EH[map(int,rsimulEH[:,rangeNB]-1)] = energyUsed_NS_EH[map(int,rsimulEH[:,rangeNB]-1)] + migration_thermo_cost_brNW_NS[rangeNB, rangeBR]
      #breedingSeasonWH[k] = 1
  energyAvailable_NS_EH = energySupply_NS_EH - energyUsed_NS_EH
  energyAvailable_NW_EH = energySupply_NW_EH - energyUsed_NW_EH
  for i in range(0,len(energyAvailable_NS_EH)):
      if energyAvailable_NS_EH[i] < 0:
          energyAvailable_NS_EH[i] = 0
      if energyAvailable_NW_EH[i] < 0:
          energyAvailable_NW_EH[i] = 0

  #rangesLatWH[k,0] = np.mean(latWH[map(int, rsimulWH[:,rangeBR]-1)])
  #rangesLatWH[k,1] = np.mean(latWH[map(int, rsimulWH[:,rangeNB]-1)])

  #Update the newly occupied hexagons
  ra = np.zeros(shape=len(lonEH))
  ra[map(int,rsimulEH[:,rangeBR]-1)] = 1
  ranges_selected_BR_EH = np.column_stack((ranges_selected_BR_EH, ra))
  ra = np.zeros(shape=len(lonEH))
  ra[map(int,rsimulEH[:,rangeNB]-1)] = 1
  ranges_selected_NB_EH = np.column_stack((ranges_selected_NB_EH, ra))

  k=k+1

#rangesLatWH = rangesLatWH[0:k,:]
#breedingSeasonWH = breedingSeasonWH[0:k]


#Compute the simulated spatial patterns
ranges_selected_BR_EH = ranges_selected_BR_EH[:,1:len(ranges_selected_BR_EH[0,:])]
ranges_selected_NB_EH = ranges_selected_NB_EH[:,1:len(ranges_selected_NB_EH[0,:])]
migra = ranges_selected_BR_EH - ranges_selected_NB_EH
res = ranges_selected_BR_EH + ranges_selected_NB_EH
migraBR = np.zeros(shape=(len(migra[:,0]), len(migra[0,:])))
migraNB = np.zeros(shape=(len(migra[:,0]), len(migra[0,:])))
residents = np.zeros(shape=(len(migra[:,0]), len(migra[0,:])))
for i in range(0,len(ranges_selected_BR_EH[:,0])):
  for j in range(0,len(ranges_selected_BR_EH[0,:])):
      if migra[i,j] > 0:
          migraBR[i,j] = 1
      if migra[i,j] < 0:
          migraNB[i,j] = 1
      if res[i,j] == 2:
          residents[i,j] = 1
richnessBRsimuEH = np.zeros(shape=len(ranges_selected_BR_EH[:,0]))
richnessNBsimuEH = np.zeros(shape=len(ranges_selected_NB_EH[:,0]))
richnessMigraBRsimuEH = np.zeros(shape=len(ranges_selected_BR_EH[:,0]))
richnessMigraNBsimuEH = np.zeros(shape=len(ranges_selected_NB_EH[:,0]))
richnessResidentsimuEH = np.zeros(shape=len(ranges_selected_NB_EH[:,0]))
for i in range(0,len(ranges_selected_BR_EH[:,0])):
    richnessBRsimuEH[i] = sum(ranges_selected_BR_EH[i,:])
    richnessNBsimuEH[i] = sum(ranges_selected_NB_EH[i,:])
    richnessMigraBRsimuEH[i] = sum(migraBR[i,:])
    richnessMigraNBsimuEH[i] = sum(migraNB[i,:])
    richnessResidentsimuEH[i] = sum(residents[i,:])
richnessDiffSimuEH = richnessBRsimuEH - richnessNBsimuEH
PropMigrSimuEH = (richnessMigraBRsimuEH + richnessMigraNBsimuEH) / (richnessMigraBRsimuEH + richnessMigraNBsimuEH + richnessResidentsimuEH)

#BRpop_EH = np.zeros(shape=len(ranges_selected_NB_EH[0,:]))
#NBpop_EH = np.zeros(shape=len(ranges_selected_NB_EH[0,:]))
#respop_EH = np.zeros(shape=len(ranges_selected_NB_EH[0,:]))
#for i in range(0,len(ranges_selected_BR_EH[0,:])):
#    BRpop_EH[i] = sum(migraBR[:,i])
#    NBpop_EH[i] = sum(migraNB[:,i])
#    respop_EH[i] = 180 - BRpop_EH[i]
#BRpops_EH = sum(BRpop_EH)
#NBpops_EH = sum(NBpop_EH)
#respops_EH = sum(respop_EH)


richnessMigraBRsimu = np.concatenate((richnessMigraBRsimuWH, richnessMigraBRsimuEH))
richnessMigraNBsimu = np.concatenate((richnessMigraNBsimuWH, richnessMigraNBsimuEH))
richnessResidentsimu = np.concatenate((richnessResidentsimuWH, richnessResidentsimuEH))
breedingrichness = np.concatenate((breeding_WH, breeding_EH)) 
nonbreedingrichness = np.concatenate((nonbreeding_WH, nonbreeding_EH)) 
residentsrichness = np.concatenate((residentsObs_WH, residentsObs_EH))

#bestfit_results = np.vstack((richnessMigraBRsimu,richnessMigraNBsimu,richnessResidentsimu,breedingrichness,nonbreedingrichness,residentsrichness)).T

#Compute Earth Mover's Distances

EMD_BR = emd_with_flow(richnessMigraBRsimu, breedingrichness, distance_matrix)
flow=0
for i in range(0,len(distance_matrix)):
    flow = flow + sum(EMD_BR[1][i])
EMD_BR2 = EMD_BR[0] / flow

EMD_NB = emd_with_flow(richnessMigraNBsimu, nonbreedingrichness, distance_matrix)
flow=0
for i in range(0,len(distance_matrix)):
    flow = flow + sum(EMD_NB[1][i])
EMD_NB2 = EMD_NB[0] / flow

EMD_res = emd_with_flow(richnessResidentsimu, residentsrichness, distance_matrix)
flow=0
for i in range(0,len(distance_matrix)):
    flow = flow + sum(EMD_res[1][i])
EMD_res2 = EMD_res[0] / flow

EMD_tot = (EMD_BR2 + EMD_NB2 + EMD_res2) / 3

#EMD_BR_WH = emd(richnessMigraBRsimuWH, breeding_WH, distance_matrix_WH) / BRpops_WH
#EMD_NB_WH = emd(richnessMigraNBsimuWH, nonbreeding_WH, distance_matrix_WH) / NBpops_WH
#EMD_res_WH = emd(richnessResidentsimuWH, residentsObs_WH, distance_matrix_WH) / respops_WH
#EMD_BR_EH = emd(richnessMigraBRsimuEH, breeding_EH, distance_matrix_EH) / BRpops_EH
#EMD_NB_EH = emd(richnessMigraNBsimuEH, nonbreeding_EH, distance_matrix_EH) / NBpops_EH
#EMD_res_EH = emd(richnessResidentsimuEH, residentsObs_EH, distance_matrix_EH) / respops_EH
#EMD_tot = pow((EMD_BR_WH * EMD_NB_WH * EMD_res_WH * EMD_BR_EH * EMD_NB_EH * EMD_res_EH), 1.0/6)
#EMD_tot = EMD_BR_WH + EMD_NB_WH + EMD_res_WH + EMD_BR_EH + EMD_NB_EH + EMD_res_EH

#f.write(str(EMD_tot))

#exit(0)

print EMD_tot


#output_data_file.write(str(mm)+","+str(rr)+","+str(tt)+","+str(ee)+","+str(EMD_BR_WH)+","+str(EMD_NB_WH)+","+str(EMD_res_WH)+"\n")
#output_data_file.close()

#np.savetxt("model_outputs/nomigra_rangesBR_WH.csv", ranges_selected_BR_WH, delimiter=',')
#np.savetxt("model_outputs/nomigra_rangesNB_WH.csv", ranges_selected_NB_WH, delimiter=',')
#np.savetxt("model_outputs/nomigra_rangesBR_EH.csv", ranges_selected_BR_EH, delimiter=',')
#np.savetxt("model_outputs/nomigra_rangesNB_EH.csv", ranges_selected_NB_EH, delimiter=',')