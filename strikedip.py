## Python Code to Calculate Strike and Dip From a Cluster of Points
# Ryan R. Maurer - 19 July 2024
import csv
import numpy as np
import math
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt2
from statsmodels.stats.weightstats import DescrStatsW

#Source data, comment out any that you dont want to process.
src = [
        'PointClouds/Cm11bSD04.csv'
       ]
dec = -9.3 #Magnetic declination
exportToCSV = True
csvName = 'Exports/csv/LC_Cm11b-1.csv'
pltitle = "Laurel Caverns - Cm11b-1"
#Empty arrays
strike = []
dip = []
area=[]
    
for clouds in src:
    rawdata = []
    
    #Get Data from CSV
    with open(clouds) as csvfile:
        reader = csv.reader(csvfile,dialect='excel')
        for row in reader:
            nums = [float(x) for x in row]
            rawdata.append(nums)
    
    #Get origin and tie in shots    
    origin = rawdata[0]
    tieInShot1 = rawdata[1]
    tieInShot2 = rawdata[2]

    #Find average tie in shot data
    tieInD = (tieInShot1[0]+tieInShot2[0])*.5
    tieInA = (tieInShot1[1]+tieInShot2[1])*.5+dec
    tieInI = (tieInShot1[2]+tieInShot2[2])*.5

    #Find displacement from tie in shots.
    tieInDelta = [tieInD*math.cos(math.radians(tieInI))*math.sin(math.radians(tieInA)),
                tieInD*math.cos(math.radians(tieInI))*math.cos(math.radians(tieInA)),
                tieInD*math.sin(math.radians(tieInI))]
    
    #Calculate the origin of the cloud.
    cloudOrigin = np.add(origin,tieInDelta)
    
    
    cloudCoords = []
    cloudData = []
    combIndex =  []
    
    #get cloud data
    for x in range(3,len(rawdata)):
        cloudData.append(rawdata[x])
    
    #for each shot, 
    for shot in cloudData:
        
        #calculate the X, Y, and Z coordinates relative to the shot origin.
        delta = [shot[0]*math.cos(math.radians(shot[2]))*math.sin(math.radians((shot[1]+dec))),
                shot[0]*math.cos(math.radians(shot[2]))*math.cos(math.radians((shot[1]+dec))),
                shot[0]*math.sin(math.radians(shot[2]))]
        
        #Add the point and origin cordinates to fix point to the global datum
        point = np.add(delta,cloudOrigin)
        
        #Add point to the cloud cluster.
        cloudCoords.append(point)
    
    #generate an array from [0, 1, 2, ... number of points]
    i = 0
    for x in cloudCoords:
        combIndex.append(i)
        i=i+1

    #Find every 3 point combination
    combs = combinations(combIndex, 3)
    
    #For each combination
    for comb in combs:
        
        #Get plane coordinates
        A = cloudCoords[comb[0]]
        B = cloudCoords[comb[1]]
        C = cloudCoords[comb[2]]
        
        #calculate plane vectors AB and AC
        AB = [B[0]-A[0],B[1]-A[1],B[2]-A[2]]
        AC = [C[0]-A[0],C[1]-A[1],C[2]-A[2]]
        
        #Find area of the plane for weighted standard deviation
        areaVec = .5*np.abs(np.add(np.cross(A,B),np.cross(C,C),np.cross(C,A)))**.5
        a = math.sqrt(sum(pow(el,2) for el in areaVec))
        area.append(a)
        
        #calculate N, E, V
        N = AB[1]*AC[2]-AC[1]*AB[2]
        E = -(AB[0]*AC[2]-AB[2]*AC[0])
        V = AB[1]*AC[0]-AB[0]*AC[1]
        #calculate raw strike and dip
        delta = math.degrees(math.acos(V/(N**2+E**2+V**2)**.5))
        if E>0 and N>0: 
            upsilon = math.degrees(math.atan(E/N))
        if E>0 and N<0: 
            upsilon = math.degrees(math.atan(-N/E))+90
        if E<0 and N<0: 
            upsilon = math.degrees(math.atan(N/E))+180
        if E<0 and N>0: 
            upsilon = math.degrees(math.atan(N/-E))+270
            
       
        #correct the strike and dip based on quadrant and append to global array
        if delta > 90 :
            dip.append(np.abs(delta-180))
            
            upsilon = upsilon+90
            if upsilon >360:
                upsilon = upsilon-360
            upsilon = np.abs(upsilon-180)
            strike.append(upsilon)
        else:
            dip.append(delta)
            upsilon = upsilon+90
            if upsilon >360:
                upsilon = upsilon-360
            strike.append(upsilon)
            
       


#Generate histogram data
# nstrike, binsstrike, strikepatches = plt.hist(strike, bins=100)
nstrike, binsstrike, strikepatches = plt.hist(strike, bins=100,range=[0,90])
# strikeModeIndex = nstrike.argmax()
strikeModeIndex = nstrike.argmax()
#Generate histogram data
# nstrike, binsstrike, strikepatches = plt.hist(strike, bins=100)
ndip, binsdip, dippatches = plt.hist(dip, bins=100, range=[0,30])
# strikeModeIndex = nstrike.argmax()
dipModeIndex = ndip.argmax()

#Plotting
print("OUTPUT")
print("------------------")
print("Number of Shots:  ", len(cloudCoords))
print("Number of planes: ", len(dip))
print("Mode Dip:         ", str((binsdip[dipModeIndex] + binsdip[dipModeIndex+1])/2))
print("Std Dip:          ", DescrStatsW(dip, weights=area, ddof=1).std)
print("Mode Strike:      ", str((binsstrike[strikeModeIndex] + binsstrike[strikeModeIndex+1])/2))
print("Std Strike:       ",  DescrStatsW(strike, weights=area, ddof=1).std)
print('Average Area:     ', np.mean(area))
plt.xlabel("degrees")
plt.title("Total Clusters: "+str(len(src))+"  Total Planes: "+str(len(dip)))
plt.suptitle(pltitle)
strikeLabel = "Strike - Mode: "+str(np.round((binsstrike[strikeModeIndex] + binsstrike[strikeModeIndex+1])/2,2))
dipLabel = "Dip - Mode: "+str(np.round((binsdip[dipModeIndex] + binsdip[dipModeIndex+1])/2))
plt.legend([strikeLabel,dipLabel])
plt.show()

if exportToCSV:
    with open(csvName, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        
        spamwriter.writerow(['SHOTS'])
        spamwriter.writerow([src])
        spamwriter.writerow(['STRIKE','DIP','AREA'])
        i = 0
        for strk in strike:
            spamwriter.writerow([strk,dip[i],area[i]])
            i= i+1
     
        

