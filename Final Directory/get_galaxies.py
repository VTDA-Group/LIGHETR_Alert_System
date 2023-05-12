#Mainly a simplified copy of HET_obs.py
#https://github.com/sjanowiecki/HET_observability
#and the ligo skymaps tutorials
#https://github.com/gw-odw/odw-2018/tree/master/skymaps

import healpy as hp # for working with HEALPix files
import numpy as np # needed for vector operations
from scipy.stats import norm # probability functions
from astropy.utils.data import download_file
from astropy.io import ascii
import argparse
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
import pdb

def parseargs():

    class GetLoc(argparse.Action):
        def __init__(self, option_strings, dest, nargs=None, **kwargs):
            if nargs is not None:
                raise ValueError("nargs not allowed")
            super(GetLoc, self).__init__(option_strings, dest, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            url = (values)
            filename = download_file(url, cache=True)
            setattr(namespace, self.dest, filename)

    parser = argparse.ArgumentParser(description='FIND GALAXIES TO OBSERVE IN TWO CATALOGS')
    parser.add_argument('--http', dest='fits', default='https://dcc.ligo.org/public/0146/G1701985/001/LALInference_v2.fits.gz', action=GetLoc, help='HTTPS link to LIGO event localization. It will download the file if not cached.')
    parser.add_argument('-cat', dest='cat', default='MANGROVE', help='Specify which catalog to use: MANGROVE or GLADE')
    args = parser.parse_args()

    return args

def cdf(pdf):
    #Calculate contour in probability
    sortedpix = np.argsort(pdf)[::-1]
    cumsum = np.cumsum(pdf[sortedpix])
    cls = np.empty_like(pdf)
    cls[sortedpix] = cumsum*100
    return cls

def get_probability_index(cat, probb, distmu, distsigma, distnorm, pixarea, nside, probability, visible_mask, never_visible_mask):
    
    '''
    This will take a pandas-read in csv file, and will return a ordered list of galaxies within that catalog that are ordered by probability map
    '''
    theta = 0.5*np.pi - cat['DEJ2000']*np.pi/180
    theta = np.asarray([float(i) for i in theta])
    phi = cat['RAJ2000']*np.pi/180
    phi = np.asarray([float(i) for i in phi])
    cls = cdf(probb)


    ipix = hp.ang2pix(nside, theta, phi)
    pix = visible_mask
    

    
    
    #accounting for probability distribution along the sky
    pixel_log_prob = np.log(probability[ipix])
    
    #accounting for distance of galaxies and distance localization from LIGO
    dist = cat['dist_Mpc'] #catalog galaxy distances in megaparsec
    c_pdf = conditional_pdf(dist,distmu[ipix],distsigma[ipix],distnorm[ipix]).tolist()
    distance_log_prob = np.log(c_pdf)
    
    
    #accounting for mass of galaxies
    stellar_mass = np.asarray(cat['M*'])

    #getting probability density by accounting for distance, galaxy stellar mass, and sky location
    logdp_dV = pixel_log_prob + distance_log_prob + np.log(stellar_mass) - np.log(pixarea)
    #logdp_dV = pixel_log_prob + distance_log_prob - np.log(pixarea)
    
    
    
    cls = cls[ipix]
    #cutting to select only 90 % confidence in position
    cattop = cat[cls<90]
    logdp_dV= logdp_dV[cls<90]
    cls = cls[cls<90]
    #Now working only with event with overall probability 99% lower than the most probable
    
    gg = np.where((logdp_dV is not np.nan) & (logdp_dV-np.max(logdp_dV) > np.log(1/100)))
    top99i = logdp_dV-np.max(logdp_dV) > np.log(1/1000)

    cattop = cattop[top99i]
    logdp_dV = logdp_dV[top99i]
    cls = cls[top99i]
    #print("cls: "+str(cls))
    #sorting by probability
    isort = np.argsort(logdp_dV)[::-1]
    
    cattop = Table.from_pandas(cattop.iloc[isort])
    logptop = logdp_dV[isort]
    cls = cls[isort]
    
    
    return cattop, logptop, cls

def write_catalog(params, savedir=''):
    fits = params['skymap_fits']
    event = params['superevent_id']
    probability = params['skymap_array']
    visible_mask = params['visible_mask']
    never_visible_mask = params['never_visible_mask']
    
    
    # Reading in the skymap prob and header
    locinfo, header = hp.read_map(fits, field=range(4), h=True)
    probb, distmu, distsigma, distnorm = locinfo
    # Getting healpix resolution and pixel area in deg^2
    npix = len(probb)
    nside = hp.npix2nside(npix)
    # Area per pixel in steradians
    pixarea = hp.nside2pixarea(nside)
    # Get the catalog
    
    
    #working with list of galaxies visble to HET
    cat1 = pd.read_csv("Glade_HET_Visible_Galaxies.csv", sep=',',header=0,dtype=np.float64)
    #plt.show()
    #print("cat1: "+str(cat1))
    cattop, logptop, cls = get_probability_index(cat1, probb, distmu, distsigma, distnorm, pixarea, nside, probability, visible_mask, never_visible_mask)
    #print("cattop: "+str(cattop))

    index = Column(name='index',data=np.arange(len(cattop)))
    logprob = Column(name='LogProb',data=logptop)
    exptime = Column(name='exptime',data=60*20*np.ones(len(cattop)))
    contour = Column(name='contour',data = cls)
    Nvis = Column(name='Nvis',data=np.ones(len(cattop)))
    cattop.add_columns([index,logprob,exptime,Nvis,contour])
    ascii.write(cattop['index','RAJ2000','DEJ2000','dist_Mpc','exptime','Nvis','LogProb','contour'], savedir+'HET_Visible_Galaxies_prob_list.dat', overwrite=True)
    
    
    return cattop,logptop

def main():




    args = parseargs()
    prob, header = hp.read_map(args.fits, h=True)
    header = dict(header)
    params = {'skymap_fits':args.fits,'skymap_array':prob,'GraceID':header['OBJECT']}
    if args.cat == 'MANGROVE' or args.cat == 'GLADE':
        write_catalog(params,args.cat)
    else:
        print('Must specify GLADE or MANGROVE as catalog.')

if __name__== "__main__":
    main()

