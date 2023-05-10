import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from scipy.stats import norm
from astropy.utils.data import download_file
from astropy.io import ascii
from astropy.table import Table
from astropy.table import Column
from astroquery.vizier import Vizier
from scipy.special import gammaincinv
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
import astropy.constants as c
import pandas as pd
from ligo.skymap.distance import conditional_pdf
import pdb
import matplotlib.pyplot as plt
from get_LST import *
from make_phaseii import *



def cdf(pdf):
    #Calculate contour in probability
    sortedpix = np.argsort(pdf)[::-1]
    cumsum = np.cumsum(pdf[sortedpix])
    cls = np.empty_like(pdf)
    cls[sortedpix] = cumsum*100
    return cls

#def get_probability_index(cat, probb, distmu, distsigma, distnorm, pixarea, nside, probability):
def get_probability_index(cat, probb, distmu, distsigma, distnorm, pixarea, nside):
    
    '''
    This will take a pandas-read in csv file, and will return a ordered list of galaxies within that catalog that are ordered by probability map
    '''
    
    print("cat1: "+str(cat))
    theta = 0.5*np.pi - cat['DEJ2000']*np.pi/180
    theta = np.asarray([float(i) for i in theta])
    print("theta: "+str(theta))
    
    phi = cat['RAJ2000']*np.pi/180
    phi = np.asarray([float(i) for i in phi])
    cls = cdf(probb)

    print("cls: "+str(cls))

    ipix = hp.ang2pix(nside, theta, phi)
    print("ipix: "+str(ipix))
    cls = cls[ipix]

    dist = cat['d']
    
    logdp_dV = np.log(conditional_pdf(dist,distmu[ipix],distsigma[ipix],distnorm[ipix]).tolist()) - np.log(pixarea)
    
    #logdp_dV = np.log(probability[ipix]) + np.log(conditional_pdf(dist,distmu[ipix],distsigma[ipix],distnorm[ipix]).tolist()) - np.log(pixarea)

    #cutting to select only 90 % confidence in position
    cattop = cat[cls<90]
    logdp_dV= logdp_dV[cls<90]
    s_lumK = 10**(-0.4*cat['B'][cls<90])
    s_lumK = s_lumK/s_lumK.sum()
    #s_lumB = 10**(-0.4*cat1['B_Abs'][cls>90])
    #s_lumB = s_lumB/s_lumB.sum()
    cls = cls[cls<90]
    #only using K for now
    logdp_dV = np.log(s_lumK) + logdp_dV

    #Now working only with event with overall probability 99% lower than the most probable
    top99i = logdp_dV-np.max(logdp_dV) > np.log(1/100)

    cattop = cattop[top99i]
    logdp_dV = logdp_dV[top99i]
    cls = cls[top99i]

    #sorting by probability
    isort = np.argsort(logdp_dV)[::-1]
    
    cattop = Table.from_pandas(cattop.iloc[isort])
    logptop = logdp_dV.iloc[isort]
    cls = cls[isort]
    
    
    return cattop, logptop, cls



def get_galaxy_list(Skymap_fits_file = "60029.41050925926/flattened_multiorder_fits_MS230326j.fits", HET_visible_galaxy_file = "Glade_HET_Visible_Galaxies.csv", All_galaxy_file = "Glade_HET_Visible_Galaxies.csv"):

    if Skymap_fits_file is not None:
        locinfo, header = hp.read_map(Skymap_fits_file, field=range(4), h=True)
        probb, distmu, distsigma, distnorm = locinfo
        # Getting healpix resolution and pixel area in deg^2
        npix = len(probb)
        nside = hp.npix2nside(npix)
        # Area per pixel in steradians
        pixarea = hp.nside2pixarea(nside)


    #read in galaxy catalog using csv file

    cat1 = pd.read_csv("Glade_HET_Visible_Galaxies.csv", sep=',',usecols = [1,2,3,4,5],names=['RAJ2000','DEJ2000','B','K','d'],header=0,dtype=np.float64)

    #cattop, logptop, cls = get_probability_index(cat1, probb, distmu, distsigma, distnorm, pixarea, nside, probability)
    cattop, logptop, cls = get_probability_index(cat1, probb, distmu, distsigma, distnorm, pixarea, nside)
    #edit catalog
    
    index = Column(name='index',data=np.arange(len(cattop)))
    logprob = Column(name='LogProb',data=logptop)
    exptime = Column(name='exptime',data=60*20*np.ones(len(cattop)))
    contour = Column(name='contour',data = cls)
    Nvis = Column(name='Nvis',data=np.ones(len(cattop)))
    cattop.add_columns([index,logprob,exptime,Nvis,contour])
    ascii.write(cattop['index','RAJ2000','DEJ2000','exptime','Nvis','LogProb','contour'], 'Test_Found_Galaxies_list.dat', overwrite=True)
    
    print("read in galaxies: "+str(cat1))
    print("ordered galaxies by probability: "+str(cattop))
    plt.scatter(cattop['RAJ2000'][:50], cattop['DEJ2000'][:50])
    plt.savefig("top 50 galaxies in list.pdf")
    
    return cattop,logptop


cattop, logptop  = get_galaxy_list()
mincontour = get_LST(savedir = '',targf = 'Test_Found_Galaxies_list.dat')
make_phaseii(lstfile = 'LSTs_Found.out', savedir = '')
