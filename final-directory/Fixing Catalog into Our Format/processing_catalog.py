import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/Users/sky5265/Documents/GitHub/Astro_Code')
from Astro_useful_funcs import *
from Analysis_useful_funcs import *
import time
import pdb


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

def get_dist_mod(dist_Mpc):
    return 5*np.log10(dist_Mpc*1.0E6)-5

def get_abs_from_app_mag(app_mag, dist_Mpc):
    return app_mag - get_dist_mod(dist_Mpc)
    

chunks=pd.read_table("Full_Glade_Catalog.txt", chunksize=1000000,sep=' ', usecols = [9,10,11,19,33], header=None)
df=pd.DataFrame()
df=pd.concat(chunks)

df.columns = ['RAJ2000','DEJ2000','B','K','dist_Mpc']

print(df)


HET_visible_galaxies = df.loc[(df['DEJ2000'] > -12) & (df['DEJ2000'] < 74)]
print(HET_visible_galaxies)
RA = HET_visible_galaxies['RAJ2000']
Dec = HET_visible_galaxies['DEJ2000']
B_app = HET_visible_galaxies['B']
K_app = HET_visible_galaxies['K']
dist_Mpc = HET_visible_galaxies['dist_Mpc']
pdb.set_trace()

B_abs = np.asarray(get_abs_from_app_mag(B_app[i], dist_Mpc[i]) for i in range(len(B_app)))
K_abs = np.asarray(get_abs_from_app_mag(K_app[i], dist_Mpc[i]) for i in range(len(K_app)))



for i in range(len(B_app)):
    HET_visible_galaxies.loc[i] = [RA[i], Dec[i], B_abs[i], K_abs[i], dist_Mpc[i]]



HET_visible_galaxies.to_csv("Glade_HET_Visible_Galaxies.csv", sep=',')


print(len(df))

#df = dd.read_csv("Full_Glade_Catalog.txt", sep = ' ')

#cat1 = pd.read_csv("Full_Glade_Catalog.txt", sep=' ')
print(df)

#cat1 = pd.read_csv("Full_Glade_Catalog.txt", sep=' ',usecols = [1,9,10,11,19,33],names=['GladeName','RAJ2000','DEJ2000','d','B_app','K_Aapp','dist_Mpc'],header=0,dtype=np.float64)



#convert apparent magnitudes into apparent magnitudes





print(B_abs)
