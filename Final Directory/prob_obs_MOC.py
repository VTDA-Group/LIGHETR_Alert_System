import healpy as hp
import astropy.coordinates
import astropy.time
import astropy.units as u
import astropy_healpix as ah
import numpy as np
#from astropy.coordinates import Angle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mhealpy import HealpixMap
import sys


maxhetdec_deg = 74
minhetdec_deg = -12
HET_loc = (-104.01472,30.6814,2025)


def get_uniq_from_ang(ra, dec, lvls, ipix):
    """
    Get uniq IDs for set of theta and phi values.
    """
    ra = ra * u.deg
    dec = dec * u.deg

    nside = ah.level_to_nside(lvls)
    match_ipix = ah.lonlat_to_healpix(ra, dec, nside, order='nested')
    i = np.flatnonzero(ipix == match_ipix)[0]
    
    return ah.level_ipix_to_uniq(lvls[i], ipix[i])

def convert_uniq_to_ra_dec(uniq):
    """
    Convert uniq IDs to ra dec values in degrees.
    """
    lvl, ipix = ah.uniq_to_level_ipix(uniq)
    nside = ah.level_to_nside(lvl)
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')
    return ra.degree, dec.degree
    
def query_strip_uniq(min_dec, max_dec, lvls, ipix):
    """
    Convert uniq IDs to ra dec values in degrees.
    """
    nside = ah.level_to_nside(lvls)
    
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')
    uniq = ah.level_ipix_to_uniq(lvls, ipix)
    
    uniq_strip = uniq[(dec.degree > min_dec) & (dec.degree < max_dec)]
            
    return np.unique(uniq_strip)
    
def query_disc_uniq(ra, dec, radius, lvls, ipix):
    """
    Convert uniq IDs to ra dec values in degrees.
    """
    nside = ah.level_to_nside(lvls)
    
    ra_all, dec_all = ah.healpix_to_lonlat(ipix, nside, order='nested')
    uniq = ah.level_ipix_to_uniq(lvls, ipix)
    
    uniq_strip = uniq[(ra_all.degree - ra)**2 + (dec_all.degree - dec)**2 <= radius**2]
    
    if len(uniq_strip) == 0: # pixels too large around disc
        return np.array([get_uniq_from_ang(ra, dec, lvls, ipix),])
        
    return np.unique(uniq_strip)

def get_night_times(t, observatory):
    """
    Get time values where Sun is down.
    """
    delta_time = np.linspace(0, 24, 1000)*u.hour
    times24 = t + delta_time
    frames24 = astropy.coordinates.AltAz(obstime = times24, location=observatory)
    sunaltazs24 = astropy.coordinates.get_sun(times24).transform_to(frames24)
    
    is_nighttime = sunaltazs24.alt<-18*u.deg
    night_times = times24[is_nighttime] - t
    night_times.format = 'sec'
    return night_times.value
    
def get_90_prob_region(m):
    """
    Uses multi-order sky map to get 90% confidence region indices.
    """
    level, ipix = ah.uniq_to_level_ipix(m['UNIQ'])
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    prob = pixel_area * m['PROBDENSITY']
    cumprob = np.cumsum(prob)
    p90i = cumprob.searchsorted(0.9)
    return np.arange(0, p90i)
    
    
def prob_observable(m, header, time, savedir, plot = True):
    """
    Determine the integrated probability contained in a gravitational-wave
    sky map that is observable with HET at a particular time. Needs hetpix.dat,
    pixels and lonlat of HET pupil, in the directory.
    """
    
    # all imports
    hetedge = np.loadtxt('hetedge.dat')
    hetpupil = np.loadtxt('hetpix.dat')

    # Obtain properties of skymap
    m.sort('PROBDENSITY', reverse=True) # MUST STAY AT TOP OF FUNCTION
    level, ipix = ah.uniq_to_level_ipix(m['UNIQ'])
    nside = ah.level_to_nside(level)
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    UNIQ = m['UNIQ']
    mplot = np.zeros(len(UNIQ))
   

    # get information about observatory
    observatory = astropy.coordinates.EarthLocation(
        lat=HET_loc[1]*u.deg,
        lon=HET_loc[0]*u.deg,
        height=HET_loc[2]*u.m
    )
    
    # information about current time
    t = astropy.time.Time(time,scale='utc',location=HET_loc)
    LST = t.sidereal_time('mean').deg
    
    # Get full pixel coverage of HET over 24 hours
    hetfullpix_uniq = query_strip_uniq(minhetdec_deg, maxhetdec_deg, level, ipix)
    hetfullpix = np.array([np.where(u == UNIQ)[0][0] for u in hetfullpix_uniq]) # TODO: how to vectorize this?


    # get CURRENT HET pupil location pixels
    HETphi = (hetpupil[:,1]+LST)%360
    HETtheta = hetpupil[:,2]
    
    newuniq = np.array([get_uniq_from_ang(HETphi[i], HETtheta[i], level, ipix) for i in range(len(HETtheta))])

    newuniq = np.unique(newuniq)
    newpix = np.array([np.where(u == UNIQ)[0][0] for u in newuniq]) # TODO: how to vectorize this?

    # Alt/az reference frame at the observatory, in this time
    frame = astropy.coordinates.AltAz(obstime=t, location=observatory)

    # Look up (celestial) spherical polar coordinates of HEALPix grid.
    ra, dec = convert_uniq_to_ra_dec(UNIQ)
    radecs = astropy.coordinates.SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        
    # Transform grid to alt/az coordinates at observatory, in this time
    altaz = radecs.transform_to(frame)

    #Get RA,DEC of the sun in this time
    sun = astropy.coordinates.get_sun(time)
    sun_altaz = sun.transform_to(frame)
    night_times = get_night_times(t, observatory)
    
    """
    # HET pixel location when the night begins
    LST_nightstart = nightstart.sidereal_time('mean').deg
    HETphi_nightstart = (hetpupil[:,1]+LST_nightstart )%360 * u.deg
    newuniq_nightstart = get_uniq_from_ang(HETtheta, HETphi_nightstart)
    newpix_nightstart = np.array([np.where(u == UNIQ)[0] for u in newuniq_nightstart])
    """
    
    # How likely is it that the (true, unknown) location of the source
    # is within the area that is visible, in this time and within 24 hours? 
    # Demand that it falls in the HETDEX pupil, that the sun is at least 18 
    # degrees below the horizon and that the airmass (secant of zenith angle 
    # approximation) is at most 2.5.

    # determine 90% probability region
    p90i = get_90_prob_region(m)

    # plotting function
    if plot:

        #SUN CIRCLE OF 18 DEGREES
        radius = 18
        uniq_sun = query_disc_uniq(sun.ra.degree, sun.dec.degree, radius, level, ipix)
        ipix_sun = np.array([np.where(u == UNIQ)[0][0] for u in uniq_sun])

        #Coloring the plot, order important here!
        mplot[:] = 1
        mplot[altaz.secz > 2.5] = 0.1
        mplot[altaz.alt < 0] = 0.99

        mplot[newpix] = 0.2
        mplot[p90i] = 0.4
        mplot[ipix_sun] = 0.6
        
        #hp.mollview(mplot, coord='C',cmap= 'nipy_spectral', cbar=False, max=1, title='HET NOW')
        #hp.graticule(local=True)
        mhealpy_map = HealpixMap(mplot, UNIQ, density = True, scheme='NUNIQ')
        mhealpy_map.plot(coord = 'C', cmap = 'nipy_spectral', cbar=False, vmin=0., vmax=1)
        ax1 = [2,4,6,8,10,12,14,16,18,20,22,24]
        for ii in ax1:
            hp.projtext(ii/24.*360-1,-5,str(ii)+'h',lonlat=True)
        ax2 = [60,30,0,-30,-60]
        for ii in ax2:
            hp.projtext(360./2,ii,'   '+str(ii)+'°',lonlat=True)
        plt.savefig(savedir+'HET_visibility_figure.pdf')

    # intersection of the 90% confidence region and what
    # HET can see NOW and over the next 24 hours
    
    mask_arraynow = np.intersect1d(p90i, newpix)
    mask_arrayfull = np.intersect1d(p90i, hetfullpix)
    
    # get exact RA/DEC intersection of 90% and HET
    ra90, dec90 = convert_uniq_to_ra_dec(UNIQ[p90i]) # IN DEGREES
    dec90HETi = (dec90 > minhetdec_deg) * (dec90 < maxhetdec_deg)
    
    dec90HET = dec90[dec90HETi]
    ra90HET = ra90[dec90HETi]
    
    
    # get exact contour of HET region edge
    hetedgef = lambda x: np.interp(x,hetedge[:,0],hetedge[:,1]) # input theta, get corresponding phi of edge
    x = (ra90HET - LST)%360 #hetedgef uses HA, so we re-adjust
    intersect_time = (x - hetedgef(dec90HET))*3600*12/180 # distance until reaching edge, converted to time in seconds
    s_start, s_end = np.min(intersect_time), np.max(intersect_time)
    
    #if the region doesn't intersect HET now
    if len(np.intersect1d(p90i,newpix)) == 0:
    
        #if the region doesn't intersect HET at all
        if len(np.intersect1d(p90i,hetfullpix)) == 0:
            return 0 , 0 , -99, 0, np.array([]), np.array([]), hetedgef(dec90HET)
 
        assert s_end > s_start
        assert s_start > 0.
        
        # This code block should never proc...
        """
        if s_start < 0:
            print("POTENTIAL BUG: wsecs < 0")
            #it's inside the pupil...
            hetedgef2 = lambda x: np.interp(x,hetedge[:,0],hetedge[:,2])
            y = (hetedgef2(theta90HET)+180)%360 -180
            x = (x+180)%360-180
            wsecs = np.min(x - y)*3600*12/180
        """
        good_observing_times = night_times[(night_times > s_start) & (night_times < s_end)]
        
    else:
        good_observing_times = night_times[night_times < s_end]
        
    if len(good_observing_times) == 0:
        return 0 , 0 , -99, 0, np.array([]), np.array([]), hetedgef(dec90HET)
        
    timetill90 = np.min(good_observing_times)/3600
    
    indiv_probs = pixel_area * m['PROBDENSITY']
    prob = np.sum(indiv_probs[mask_arraynow]).value
    probfull = np.sum(indiv_probs[mask_arrayfull]).value
    
    # NOTE: this is a modification of the skymap. Only do at end
    m['PROBDENSITY'][np.setdiff1d(np.arange(len(m)),mask_arrayfull,assume_unique=True)] = np.min(m['PROBDENSITY'])

    # Done!
    return prob, probfull, timetill90, m, mask_arraynow, mask_arrayfull, hetedgef(dec90HET) # import the intersections too


def main():
        skymap, header = hp.read_map(sys.argv[1],
                                     h=True, verbose=False)
        header = {'GraceID': 'TEST'}
        time = astropy.time.Time.now()
        prob, probfull, timetill90, m = prob_observable(skymap, header, time, plot=True)
        return timetill90
        
if __name__=="__main__":
    main()
