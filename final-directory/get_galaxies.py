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
import astropy_healpix as ah
from astropy.table import QTable
from astropy.io import fits

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

def get_90_prob_region(m):
    """
    Uses multi-order sky map to get 90% confidence region indices.
    """
    level, ipix = ah.uniq_to_level_ipix(m['UNIQ'])
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    prob = pixel_area * m['PROBDENSITY']
    cumprob = np.cumsum(prob)
    p90i = cumprob.searchsorted(0.9)
    return p90i, cumprob

def get_idx_from_ang(ra, dec, lvls, ipix):
    """
    Get uniq IDs for set of theta and phi values.
    """
    ra = ra * u.deg
    dec = dec * u.deg

    nside = ah.level_to_nside(lvls)
    match_ipix = ah.lonlat_to_healpix(ra, dec, nside, order='nested')
    i = np.flatnonzero(ipix == match_ipix)[0]
    
    return i
    
def get_probability_index(
    cat,
    p90i,
    cumprob,
    distmu,
    distsigma,
    distnorm,
    pixarea,
    lvls,
    ipix,
    probability
):
    
    '''
    This will take a pandas-read in csv file, and will return a ordered list of galaxies within that catalog that are ordered by probability map
    '''
    
    dec = np.array(cat['DEJ2000']).astype(float)
    ra = np.array(cat['RAJ2000']).astype(float)
    
    gal_i = np.array([get_idx_from_ang(ra[i], dec[i], lvls, ipix) for i in range(len(ra))])

    dist = cat['d']
    logdp_dV = np.log(probability[gal_i]) + np.log(conditional_pdf(dist,distmu[gal_i],distsigma[gal_i],distnorm[gal_i]).tolist()) - np.log(pixarea)

    #cutting to select only 90 % confidence in position
    cattop = cat[gal_i < p90i]
    logdp_dV = logdp_dV[gal_i < p90i]
    gal_cls = cumprob[gal_i]
    cls = gal_cls[gal_i < p90i]
    s_lumK = 10**(-0.4*cattop['B'])
    s_lumK = s_lumK/s_lumK.sum()
    #s_lumB = 10**(-0.4*cat1['B_Abs'][cls>90])
    #s_lumB = s_lumB/s_lumB.sum()
    
    #only using K for now
    logdp_dV = np.log(s_lumK) + logdp_dV
    
    return cattop, logdp_dV, cls

def write_catalog(params, savedir=''):
    fits_f = params['skymap_fits']
    event = params['superevent_id']
    probability = params['skymap_array']
    
    """
    # Reading in the skymap prob and header
    locinfo, header = hp.read_map(fits, field=range(4), h=True)
    probb, distmu, distsigma, distnorm = locinfo
    """
    
    m = QTable.read(fits_f)
    
    m.sort('PROBDENSITY', reverse=True) # MUST STAY AT TOP OF FUNCTION
    distmu = m['DISTMU']
    distsigma = m['DISTSIGMA']
    distnorm = m['DISTNORM']
    
    header = fits.open(fits_f)[1].header
    
    level, ipix = ah.uniq_to_level_ipix(m['UNIQ'])
    nside = ah.level_to_nside(level)
    pixarea = ah.nside_to_pixel_area(ah.level_to_nside(level))
    UNIQ = m['UNIQ']
    
    p90i, cumprob = get_90_prob_region(m)
    
    print("Beginning CSV read to pandas")
    #working with list of galaxies visble to HET
    reader = pd.read_csv("Glade_HET_Visible_Galaxies.csv", chunksize=100000, sep=',',usecols = [1,2,3,4,5],names=['RAJ2000','DEJ2000','B','K','d'],header=0,dtype=np.float64)
    #plt.show()
    cattop = np.array([])
    logdp_dV = np.array([])
    cls = np.array([])
    
    for chunk in reader:
        ct, logptop, con = get_probability_index(chunk, p90i, cumprob, distmu, distsigma, distnorm, pixarea, level, ipix, probability)
        cattop = np.append(cattop, ct)
        logdp_dV = np.append(logdp_dV, logptop)
        cls = np.append(cls, con)
        
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
    
    index = Column(name='index',data=np.arange(len(cattop)))
    logprob = Column(name='LogProb',data=logptop)
    exptime = Column(name='exptime',data=60*20*np.ones(len(cattop)))
    contour = Column(name='contour',data = cls)
    Nvis = Column(name='Nvis',data=np.ones(len(cattop)))
    cattop.add_columns([index,logprob,exptime,Nvis,contour])
    ascii.write(cattop['index','RAJ2000','DEJ2000','exptime','Nvis','LogProb','contour'], savedir+'HET_Visible_Galaxies_prob_list.dat', overwrite=True)
    
    #should find the number of galaxies that will be visible to HET, compared to the number of total galaxies within the region
    num_galaxies_visible_HET = len(index)
    
    #divide this number by the number of galaxies within the region corresponding to the skymap
    
    
    return cattop,logptop,num_galaxies_visible_HET

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

