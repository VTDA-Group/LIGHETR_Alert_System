import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/Users/sky5265/Documents/GitHub/Astro_Code')
from Astro_useful_funcs import *
from Analysis_useful_funcs import *
import time
import pdb

'''
I want to have these parameters in my catalog of HET-visible galaxies:
- Galaxy Name
- RA (degrees)
- Dec (degrees)
- Distance in Mpc
- Apparent B magnitude
- Absolute B magnitude
- Stellar Mass (units of 10^10 solar mass)

I'll read in these columns:
indexing from 0
col1 = glade no
col9 = RA (degrees)
col10 = Dec (degrees)
col11 = apparent B mag
col14 = absolute B mag
col19 = apparent K mag
col29 = CMB frame redshift
col33 = distance in Mpc
col35 = dist flag
col36 = stellar mass

'''

chunks=pd.read_table("GLADE_full.txt", chunksize=1000000,sep=' ', usecols = [1,8,9,10, 13,18, 28,32,34, 35], header=None)
df=pd.DataFrame()
df=pd.concat(chunks)

df.columns = ['Glade_Name','RAJ2000','DEJ2000','B_app', 'B_abs','K_app', 'z_CMB','dist_Mpc','dist_flag','M*']

print("There are "+str(len(df))+" total galaxies in this catalog")
#print("There are "+str(len(df.loc[(df['DEJ2000'] > -12) & (df['DEJ2000'] < 74)]))+" total galaxies in this catalog within the HET field of view")

HET_visible_galaxies = df.loc[(df['DEJ2000'] > -12) & (df['DEJ2000'] < 74) & (df['dist_flag']>0) & (df['dist_Mpc'].notnull()) & (df['M*'].notnull())]
All_visible_galaxies = df.loc[(df['dist_flag']>0) & (df['dist_Mpc'].notnull()) & (df['M*'].notnull())]
print("There are "+str(len(HET_visible_galaxies))+" HET visible galaxies")

print(HET_visible_galaxies)

RA = HET_visible_galaxies['RAJ2000']
Dec = HET_visible_galaxies['DEJ2000']

#confirm that RA ranges from 0 to 360
print("max RA = "+str(max(RA)))
print("min RA = "+str(min(RA)))

#confirm dec ranges from -90 to 90
print("max dec = "+str(max(Dec)))
print("min dec = "+str(min(Dec)))

HET_visible_galaxies.to_csv("Glade_HET_Visible_Galaxies.csv", sep=',')
All_visible_galaxies.to_csv("Glade_Visible_Galaxies.csv", sep=',')
