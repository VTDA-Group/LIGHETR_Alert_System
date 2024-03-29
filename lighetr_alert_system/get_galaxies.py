#Mainly a simplified copy of HET_obs.py
#https://github.com/sjanowiecki/HET_observability
import os
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

from lighetr_alert_system.utils import *

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
    clz = np.empty_like(pdf)
    clz[sortedpix] = cumsum*100
    return clz

def get_probability_index(cat, cls_all, distmu, distsigma, distnorm, pixarea, nside, probability):
    
    '''
    This will take a pandas-read in csv file, and will return a ordered list of galaxies within that catalog that are ordered by probability map
    '''
    theta = 0.5*np.pi - cat['DEJ2000']*np.pi/180
    phi = cat['RAJ2000']*np.pi/180

    distinfo = np.array([distmu, distsigma, distnorm]).T
    ipix_reduced = hp.ang2pix(nside, theta, phi)
    
    # get all pixels of 90% confidence region
    #all_banana_pixs = np.where((probability > 0.) & (cls_all < 90.))[0]
    
    unique_ipix, pixel_mapping = np.unique(ipix_reduced, return_inverse=True)
    pixel_prob = probability[unique_ipix]
    distinfo = distinfo[ipix_reduced]
    
    # missing prob is probability occupied by < 100% completeness, galaxies we cannot see
    gal_probs, missing_prob = distribute_pixel_prob(
                                                pixel_prob,
                                                pixel_mapping,
                                                cat,
                                                distinfo
                                            )
    gal_sum = np.sum(gal_probs)
    print("norm", gal_sum, missing_prob)
    # we can normalize the non-missing prob to sum up to (1 - missing_prob)
    gal_probs /= (gal_sum + missing_prob)
    missing_prob /= (gal_sum + missing_prob)
    print("norm", np.sum(gal_probs), missing_prob)
    logdp_dV = np.log(gal_probs)
    #logdp_dV= logdp_dV[cls<90]
    clz = cls_all[ipix_reduced]
    
    top99i = np.where((logdp_dV-np.nanmax(logdp_dV)) > np.log(1/100))[0]
    
    cattop = cat.iloc[top99i]
    logptop = logdp_dV[top99i]
    clz = clz[top99i]
    
    #sorting by probability
    isort = np.argsort(logptop)[::-1]
    
    num_keep = 1000 # number of galaxies in output
    
    if len(logptop) > num_keep:
        isort = isort[:num_keep]

    logptop = logptop[isort]
    clz = clz[isort]
    cattop = cattop.iloc[isort]
    
    return cattop, logptop, clz


def aggregate_galaxies(chunksize, probb, distmu, distsigma, distnorm, pixarea, nside, probability):
    """
    Helper function to attempt galaxy file imports with different chunksizes.
    """
    
    dir_path = os.path.dirname(os.path.realpath(__file__))
    reader = pd.read_csv(
        os.path.join(
            dir_path,
            "../data/Glade_Visible_Galaxies.csv"
        ), chunksize=chunksize, sep=',',header=0
    )

    clz = cdf(probb)
    cattop = None
    for chunk in reader:
        dists = chunk['dist_Mpc']
        theta = 0.5*np.pi - chunk['DEJ2000']*np.pi/180
        phi = chunk['RAJ2000']*np.pi/180
        ipix = hp.ang2pix(nside, theta, phi)
        clz_red = clz[ipix]

        # 3 sigma cut
        keep_idxs = np.where(
                        ((dists - distmu[ipix])**2 <= (3 * distsigma[ipix])**2) \
                        & (probability[ipix] > 0.) \
                        & (clz_red < 90.) \
                        #& (dists < 400.) # Mpc
                    )[0] # dont include masked out pixels
        #keep_idxs = np.where(probability[ipix] > 0.)[0]
        
        if cattop is None:
            cattop = chunk.iloc[keep_idxs]
        else:
            cattop = pd.concat((cattop, chunk.iloc[keep_idxs]), axis=0)
            
        
    #accounting for probability distribution along the sky
    cattop, l, c = get_probability_index(cattop, clz, distmu, distsigma, distnorm, pixarea, nside, probability)
        
    return cattop, l, c
 

def write_catalogs(observatories, savedir=''):
    
    alert = observatories[0].alert
    
    fits = alert.singleorder_file # unmodified skymap
    # Reading in the skymap prob and header
    locinfo, header = hp.read_map(fits, field=range(4), h=True)
    probb, distmu, distsigma, distnorm = locinfo
    # Getting healpix resolution and pixel area in deg^2
    npix = len(probb)
    nside = hp.npix2nside(npix)
    # Area per pixel in steradians
    pixarea = hp.nside2pixarea(nside)
        
    cattop_all = []
    logptop_all = []
    
    for i, obs in enumerate(observatories):
        if obs.frac_visible == 0.0:
            cattop_all.append([])
            logptop_all.append([])
            continue
            
        probability = obs.skymap
        #working with list of galaxies visble to HET
        chunksize = 50000 # initial chunksize
        while True:
            #try:
            cattop, logptop, clz = aggregate_galaxies(
                chunksize,
                probb,
                distmu,
                distsigma,
                distnorm,
                pixarea,
                nside,
                probability,
            )
            break
            #except:
            #    chunksize /= 5 # if crashed from chunksize
            #    if chunksize < 1000:
            #        raise ValueError("Failure to generate galaxy catalog")

        index = Column(name='index',data=np.arange(len(cattop)))
        logprob = Column(name='LogProb',data=logptop)
        exptime = Column(name='exptime',data=60*20*np.ones(len(cattop)))
        contour = Column(name='contour',data = clz)
        Nvis = Column(name='Nvis',data=np.ones(len(cattop)))

        cattop = Table.from_pandas(cattop)
        cattop.add_columns([index,logprob,exptime,Nvis,contour])

        ascii.write(
            cattop['index','RAJ2000','DEJ2000','dist_Mpc', 'M*', 'exptime','Nvis','LogProb','contour'], 
            os.path.join(
                savedir, f'{observatories[i].name}_Visible_Galaxies_prob_list.dat'
            ),
            overwrite=True
        )
        
        cattop_all.append(cattop)
        logptop_all.append(logptop)

    return cattop_all, logptop_all


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

