#Mainly a simplified copy of HET_obs.py
#https://github.com/sjanowiecki/HET_observability

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
from utils import *

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

def get_probability_index(cat, probb, distmu, distsigma, distnorm, pixarea, nside, probability):
    
    '''
    This will take a pandas-read in csv file, and will return a ordered list of galaxies within that catalog that are ordered by probability map
    '''
    theta = 0.5*np.pi - cat['DEJ2000']*np.pi/180
    theta = np.asarray([float(i) for i in theta])
    phi = cat['RAJ2000']*np.pi/180
    phi = np.asarray([float(i) for i in phi])
    cls = cdf(probb)

    distinfo = np.array([distmu, distsigma, distnorm]).T
    
    #accounting for probability distribution along the sky
    ipix = hp.ang2pix(nside, theta, phi)
    
    ipix_reduced = ipix[probability[ipix] > 0] # speedup
    
    if len(ipix_reduced) == 0:
        return np.array([]), np.array([]), np.array([])
        
    cat = cat[probability[ipix] > 0]
   
    pixel_binary_matrix, unique_ipix = get_gal_binary_matrix(ipix_reduced)
    pixel_prob = probability[unique_ipix]
    distinfo = distinfo[ipix_reduced]
    
    gal_probs = distribute_pixel_prob(pixel_prob, pixel_binary_matrix, cat, distinfo)
    logdp_dV = np.log(gal_probs)
    

    cls = cls[ipix_reduced]
    #cutting to select only 90 % confidence in position
    cattop = cat.iloc[cls<90]
    logdp_dV= logdp_dV[cls<90]
    cls = cls[cls<90]
    #Now working only with event with overall probability 99% lower than the most probable
    
    #top99i = (~np.isnan(logdp_dV)) & (logdp_dV-np.max(logdp_dV) > np.log(1/100))
    #top99i = ~np.isnan(logdp_dV)
    top99i = (logdp_dV-np.nanmax(logdp_dV)) > np.log(1/100)
    
    cattop = cattop.iloc[top99i]
    
    logdp_dV = logdp_dV[top99i]
    cls = cls[top99i]
    
    # sorting happens after aggregation anyways
    """
    #sorting by probability
    isort = np.argsort(logdp_dV)[::-1]
    
    logptop = logdp_dV[isort]
    cls = cls[isort]
    cattop = cattop.iloc[isort]
    """
    
    return cattop, logptop, cls


def aggregate_galaxies(chunksize, probb, distmu, distsigma, distnorm, pixarea, nside, probability, HET_specific_constraints):
    """
    Helper function to attempt galaxy file imports with different chunksizes.
    """
    if HET_specific_constraints:
        reader = pd.read_csv("Glade_HET_Visible_Galaxies.csv", chunksize=chunksize, sep=',',header=0,dtype=np.float64)
    else:
        reader = pd.read_csv("Glade_Visible_Galaxies.csv", chunksize=chunksize, sep=',',header=0,dtype=np.float64)
    #plt.show()
    
    ras = np.array([])
    decs = np.array([])
    dists = np.array([])
    logptop = np.array([])
    cls = np.array([])
    
    for chunk in reader:
        cattop, l, c = get_probability_index(chunk, probb, distmu, distsigma, distnorm, pixarea, nside, probability)
        if len(cattop) == 0:
            continue
        
        ras = np.append(ras, np.array(cattop['RAJ2000']))
        decs = np.append(decs, np.array(cattop['DEJ2000']))
        dists = np.append(dists, np.array(cattop['dist_Mpc']))
        logptop = np.append(logptop, l)
        cls = np.append(cls, c)
        
    return ras, decs, dists, logptop, cls
 
 
def write_catalog(params, savedir='', HET_specific_constraints = True):
    fits = params['skymap_fits']
    
    if HET_specific_constraints:
        probability = params['skymap_array_HET']
    else:
        probability = params['skymap_array']
    
    
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
    chunksize = 50000 # initial chunksize
    while True:
        try:
            ras, decs, dists, logptop, cls = aggregate_galaxies(
                chunksize,
                probb,
                distmu,
                distsigma,
                distnorm,
                pixarea,
                nside,
                probability,
                HET_specific_constraints
            )
            break
        except:
            chunksize /= 5 # if crashed from chunksize

    # 90% confidence cut redundant, removed
    
    top99i = (logptop is not np.nan) & (logptop-np.max(logptop) > np.log(1/100))

    logptop = logptop[top99i]
    cls = cls[top99i]
    ras = ras[top99i]
    decs = decs[top99i]
    dists = dists[top99i]
    
    #sorting by probability
    isort = np.argsort(logptop)[::-1]
    
    num_keep = 1000 # number of galaxies in output
    logptop = logptop[isort][:num_keep]
    cls = cls[isort][:num_keep]
    ras = ras[isort][:num_keep]
    decs = decs[isort][:num_keep]
    dists = dists[isort][:num_keep]
    
    
    index = Column(name='index',data=np.arange(len(ras)))
    ra_col = Column(name='RAJ2000',data=ras)
    dec_col = Column(name='DEJ2000',data=decs)
    dist_col = Column(name='dist_Mpc',data=dists)
    logprob = Column(name='LogProb',data=logptop)
    exptime = Column(name='exptime',data=60*20*np.ones(len(ras)))
    contour = Column(name='contour',data = cls)
    Nvis = Column(name='Nvis',data=np.ones(len(ras)))
    cattop = Table()
    cattop.add_columns([index,ra_col, dec_col, dist_col,logprob,exptime,Nvis,contour])
    
    if HET_specific_constraints:
        ascii.write(cattop['index','RAJ2000','DEJ2000','dist_Mpc','exptime','Nvis','LogProb','contour'], savedir+'HET_Visible_Galaxies_prob_list.dat', overwrite=True)
    else:
        ascii.write(cattop['index','RAJ2000','DEJ2000','dist_Mpc','exptime','Nvis','LogProb','contour'], savedir+'Visible_Galaxies_prob_list.dat', overwrite=True)

    return cattop, logptop

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

