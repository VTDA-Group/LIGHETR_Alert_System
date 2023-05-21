import healpy as hp
import astropy.coordinates
import astropy.time
import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
from astropy.time import TimeDelta


maxdec = 90
mindec = -90
mindec_rad = (90-maxdec)*np.pi/180
maxdec_rad = (90-mindec)*np.pi/180


loc = (-104.01472,30.6814,2025)

def make_visibility_figure(savedir, time, m, plot = True):
    """
    This function makes the visibility figure at a given time.
    """
    
    # Determine resolution of sky map
    mplot = np.copy(m)
    npix = len(m)
    nside = hp.npix2nside(npix)
    
    fullpix = hp.query_strip(nside, mindec_rad, \
                            maxdec_rad)
    observatory = astropy.coordinates.EarthLocation(lat=loc[1]*u.deg, lon=loc[0]*u.deg, height=loc[2]*u.m)

    # Find pixels of HET pupil in this time
    t = astropy.time.Time(time,scale='utc',location=loc)
    LST = t.sidereal_time('mean').deg
    
    
    # Alt/az reference frame at the observatory, in this time
    frame = astropy.coordinates.AltAz(obstime=t, location=observatory)

    # Look up (celestial) spherical polar coordinates of HEALPix grid.
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    # Convert to RA, Dec.
    radecs = astropy.coordinates.SkyCoord(
        ra=phi*u.rad, dec=(0.5*np.pi - theta)*u.rad)
        
    
    # Transform grid to alt/az coordinates at observatory, in this time
    altaz = radecs.transform_to(frame)

    #Get RA,DEC of the sun in this time
    sun = astropy.coordinates.get_sun(time)
    #Get RA,DEC of the moon in this time
    moon = astropy.coordinates.get_moon(time)
    
    
    # Where is the sun in the Texas sky, in this time?
    sun_altaz = sun.transform_to(frame)
    # Where is the moon in the Texas sky, in this time?
    moon_altaz = moon.transform_to(frame)
    
    
    msortedpix = np.flipud(np.argsort(m))
    cumsum = np.cumsum(m[msortedpix])
    cls = np.empty_like(m)
    cls[msortedpix] = cumsum*100
    p90i = np.where(cls <= 90)[0]
    
    
    p90i_visible = None
    
    if sun_altaz.alt < -18*u.deg:
        p90i_visible = np.intersect1d(p90i,fullpix)

    
    if plot:
        if not os.path.exists(savedir+"visibilities/"):
            os.mkdir(savedir+"visibilities/")
        
        #SUN CIRCLE OF 18 DEGREES
        radius_sun = 18
        phis = Angle(sun.ra).radian
        thetas = 0.5*np.pi-Angle(sun.dec).radian
        radius = np.deg2rad(radius_sun)
        xyz = hp.ang2vec(thetas, phis)
        ipix_sun = hp.query_disc(nside, xyz, radius)

        #Moon CIRCLE OF 5 DEGREES
        #Can try: https://stackoverflow.com/questions/70534054/plot-illuminated-percentage-of-the-moon, to get what fraction of the moon is covered
        radius_moon = 5
        phis = Angle(moon.ra).radian
        thetas = 0.5*np.pi-Angle(moon.dec).radian
        radius = np.deg2rad(radius_moon)
        xyz = hp.ang2vec(thetas, phis)
        ipix_moon = hp.query_disc(nside, xyz, radius)
        

        #Coloring the plot, order important here!
        mplot[:] = 1 #background color, everything is the grey color
        mplot[altaz.secz > 2.5] = 0.99 # salmon-region where airmass is too high
        mplot[altaz.alt < 0] = 0.99 #salmon-region below the horizon
        
        
        mplot[p90i] = 0.4 #blue, location of 90% region contours
        if sun_altaz.alt < -18*u.deg:
            mplot[p90i_visible] = 0.5 #location of 90% region contours that are visible while the sun is down
        mplot[ipix_sun] = 0.8 #yellow, location of the sun
        mplot[ipix_moon] = 0.3 #location of the moon. Note this doesn't take into account the actual brightness of the moon. It just gives the location
        hp.mollview(mplot, coord='C',cmap= 'nipy_spectral', cbar=False, max=1, title='HET NOW')
        hp.graticule(local=True)
        ax1 = [2,4,6,8,10,12,14,16,18,20,22,24]
        for ii in ax1:
            hp.projtext(ii/24.*360-1,-5,str(ii)+'h',lonlat=True)
        ax2 = [60,30,0,-30,-60]
        for ii in ax2:
            hp.projtext(360./2,ii,'   '+str(ii)+'°',lonlat=True)
            
        time_hr = time.datetime.hour
        plt.savefig(savedir+"visibilities/"+'visibility_figure_'+str(round(time_hr, 2))+'_.pdf', bbox_inches = 'tight')
    
    
    return p90i_visible

    


def prob_observable(m, header, time, savedir, plot = True):

    
    # Determine resolution of sky map
    mplot = np.copy(m)
    npix = len(m)
    nside = hp.npix2nside(npix)

    # Geodetic coordinates of MacDonald Obs
    fullpix = hp.query_strip(nside, mindec_rad, \
                            maxdec_rad)

    observatory = astropy.coordinates.EarthLocation(
        lat=loc[1]*u.deg, lon=loc[0]*u.deg, height=loc[2]*u.m)

    # Find pixels of HET pupil in this time
    t = astropy.time.Time(time,scale='utc',location=loc)
    LST = t.sidereal_time('mean').deg

    
    
    # Alt/az reference frame at the observatory, in this time
    frame = astropy.coordinates.AltAz(obstime=t, location=observatory)

    # Look up (celestial) spherical polar coordinates of HEALPix grid.
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    # Convert to RA, Dec.
    radecs = astropy.coordinates.SkyCoord(
        ra=phi*u.rad, dec=(0.5*np.pi - theta)*u.rad)
        
    
    # Transform grid to alt/az coordinates at observatory, in this time
    altaz = radecs.transform_to(frame)

    #Get RA,DEC of the sun in this time
    sun = astropy.coordinates.get_sun(time)
    #Get RA,DEC of the moon in this time
    moon = astropy.coordinates.get_moon(time)
    
    
    # Where is the sun in the sky, in this time?
    sun_altaz = sun.transform_to(frame)
    # Where is the moon in the sky, in this time?
    moon_altaz = moon.transform_to(frame)
    
    msortedpix = np.flipud(np.argsort(m))
    cumsum = np.cumsum(m[msortedpix])
    cls = np.empty_like(m)
    cls[msortedpix] = cumsum*100
    p90i = np.where(cls <= 90)
    
    
    
    #SUN CIRCLE OF 18 DEGREES
    radius_sun = 18
    phis = Angle(sun.ra).radian
    thetas = 0.5*np.pi-Angle(sun.dec).radian
    radius = np.deg2rad(radius_sun)
    xyz = hp.ang2vec(thetas, phis)
    ipix_sun = hp.query_disc(nside, xyz, radius)

    #Moon CIRCLE OF 5 DEGREES
    radius_moon = 5
    phis = Angle(moon.ra).radian
    thetas = 0.5*np.pi-Angle(moon.dec).radian
    radius = np.deg2rad(radius_moon)
    xyz = hp.ang2vec(thetas, phis)
    ipix_moon = hp.query_disc(nside, xyz, radius)

    #Coloring the plot, order important here!
    mplot[:] = 1 #background color, everything is the grey color
    mplot[altaz.secz > 2.5] = 0.99 # salmon-region where airmass is too high
    mplot[altaz.alt < 0] = 0.99 #salmon-region below the horizon
    
    
    mplot[p90i] = 0.4 #blue, location of 90% region contours
    #for i in range(24*6):
    delta_time = np.linspace(0, 24, 24)*u.hour
    for dt in delta_time:
        future_time = time + dt
    
        p90i_visible = make_visibility_figure(savedir=savedir, m = m, time = future_time, plot = False)
        
        if p90i_visible is not None:
            #print("m of the visible region is to be: "+str(m[p90i_HET_visible])+", mplot is being set to 0.5")
            mplot[p90i_visible] = 0.5
            pass
    
    #find mask of where mplot still is 0.4, this is the part of the 90% region that is never visible to HET at nighttime, because it was never set to 0.5
    
    
    never_visible_mask = np.where(mplot == 0.4)[0]
    mplot[never_visible_mask] = 0.7
    visible_mask = np.where(mplot == 0.5)[0]
    
    frac_visible = len(visible_mask)/(len(never_visible_mask) + len(visible_mask))
    
    
    #m[never_visible_mask] = 0.0
    m[mplot!=0.5] = 0.0
    
    
    
        

    
    mplot[ipix_sun] = 0.8 #yellow, location of the sun
    mplot[ipix_moon] = 0.3 #location of the moon. Note this doesn't take into account the actual brightness of the moon. It just gives the location
    hp.mollview(mplot, coord='C',cmap= 'nipy_spectral', cbar=False, max=1, title='Visible NOW')
    hp.graticule(local=True)
    ax1 = [2,4,6,8,10,12,14,16,18,20,22,24]
    for ii in ax1:
        hp.projtext(ii/24.*360-1,-5,str(ii)+'h',lonlat=True)
    ax2 = [60,30,0,-30,-60]
    for ii in ax2:
        hp.projtext(360./2,ii,'   '+str(ii)+'°',lonlat=True)
        
    plt.savefig(savedir+"Full_Visibility.pdf")
    

    

    

    # Done!
    
    
    
    return m, frac_visible

def main():
        skymap, header = hp.read_map(sys.argv[1],
                                     h=True, verbose=False)
        header = {'GraceID': 'TEST'}
        time = astropy.time.Time.now()
        prob, probfull, timetill90, m = prob_observable(skymap, header, time, plot=True)
        print(timetill90)
        return timetill90
if __name__=="__main__":
    main()
