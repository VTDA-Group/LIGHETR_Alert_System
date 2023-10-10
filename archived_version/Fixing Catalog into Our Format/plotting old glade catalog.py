import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/Users/sky5265/Documents/GitHub/Astro_Code')
from Astro_useful_funcs import *
from Analysis_useful_funcs import *
import time
import matplotlib.pyplot as plt


#I'll be making 2 of these files
	#one for galaxies that are only visible to HET
	#one for all galaxies 
 

#read in the file. I only need columns:
'''
indexing from 0
col1 = glade no
col9 = RA (degrees)
col10 = Dec (degrees)
col11 = apparent B mag
col19 = apparent K mag
col33 = distance in Mpc
'''

cat1 = pd.read_csv("./GLADE2.3HETd.csv", sep=',',usecols = [1,2,3,4,5],names=['RAJ2000','DEJ2000','d','B_Abs','K_Abs'],header=0,dtype=np.float64)

plt.scatter(cat1['RAJ2000'], cat1['DEJ2000'])
plt.xlabel("RA")
plt.ylabel("Dec")
plt.savefig("3HETd catalog.pdf")

print("Max declination: "+str(max(cat1['DEJ2000'])))
print("Min declination: "+str(min(cat1['DEJ2000'])))
print("Max RA: "+str(max(cat1['RAJ2000'])))
print("Min RA: "+str(min(cat1['RAJ2000'])))
