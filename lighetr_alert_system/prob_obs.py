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

def helper_function(observatory, time, m):
    mplot = np.copy(m)
    npix = len(m)
    nside = hp.npix2nside(npix)
 
    # Find pixels of HET pupil in this time
    t = astropy.time.Time(time,scale='utc',location=observatory.loc)
    LST = t.sidereal_time('mean').deg
    
    if observatory.pupil is not None:
        phi_obs = ((observatory.pupil[:,1]+LST)%360)*np.pi/180
        theta_obs = (90-observatory.pupil[:,2])*np.pi/180
        newpix = hp.ang2pix(nside, theta_obs, phi_obs)
    else:
        newpix = hp.query_strip(
            nside,
            observatory.min_dec_rad,
            observatory.max_dec_rad
        )
    
    # Alt/az reference frame at the observatory, in this time
    frame = astropy.coordinates.AltAz(obstime=t, location=observatory.loc)

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
    moon = astropy.coordinates.get_body("moon", time)
    
    
    # Where is the sun in the Texas sky, in this time?
    sun_altaz = sun.transform_to(frame)
    # Where is the moon in the Texas sky, in this time?
    moon_altaz = moon.transform_to(frame)
    
    msortedpix = np.flipud(np.argsort(m))
    cumsum = np.cumsum(m[msortedpix])
    cls = np.empty_like(m)
    cls[msortedpix] = cumsum*100
    if cls[msortedpix[0]] > 90:
        p90i = [0,]
    else:
        p90i = np.where(cls <= 90)[0]

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
    mplot[newpix] = 0.2 #black, color of HET pixels
    mplot[p90i] = 0.4 #blue, location of 90% region contours

    if sun_altaz.alt < -18*u.deg:
        p90i_obs_visible = np.intersect1d(p90i, newpix)
        mplot[p90i_obs_visible] = 0.5
        
    mplot[altaz.secz > 2.5] = 0.99 # salmon-region where airmass is too high
    mplot[altaz.alt < 0] = 0.99 #salmon-region below the horizon
    mplot[ipix_sun] = 0.8 #yellow, location of the sun
    mplot[ipix_moon] = 0.3 #location of the moon. Note this doesn't take into account the actual brightness of the moon. It just gives the location
        
    return mplot, p90i
    
    
def make_visibility_figure(savedir, observatory, mplot, time, combined=False):
    """
    This function makes the visibility figure at a given time.
    """
    if not os.path.exists(os.path.join(savedir,f"{observatory.name}_visibilities/")):
        os.mkdir(savedir+f"{observatory.name}_visibilities/")

    hp.mollview(mplot, coord='C',cmap= 'nipy_spectral', cbar=False, max=1, title=f'{observatory.name}')
    hp.graticule(local=True)
    ax1 = [2,4,6,8,10,12,14,16,18,20,22,24]
    for ii in ax1:
        hp.projtext(ii/24.*360-1,-5,str(ii)+'h',lonlat=True)
    ax2 = [60,30,0,-30,-60]
    for ii in ax2:
        hp.projtext(360./2,ii,'   '+str(ii)+'°',lonlat=True)

    time_hr = time.datetime.hour
    
    if combined:
        plt.savefig(
            os.path.join(
                savedir,
                f"{observatory.name}_Full_Visibility.pdf"
            )
        )
    else:
        plt.savefig(
            os.path.join(
                savedir,
                f"{observatory.name}_visibilities",
                f"{observatory.name}_visibilities_visibility_figure_{round(time_hr, 2)}_.pdf",
            ),
            bbox_inches = 'tight'
        )
    plt.close()

    

def prob_observable(m, observatory, time, savedir, plot = True, plot_timestamps = False):
    t = astropy.time.Time(time,scale='utc',location=observatory.loc)
    delta_time = np.linspace(0, 24, 100)*u.hour
    times24 = t + delta_time
    # not used currently
    """
    frames24 = astropy.coordinates.AltAz(obstime = times24, location=observatory)
    sunaltazs24 = astropy.coordinates.get_sun(times24).transform_to(frames24)
    timetilldark = 0*u.hour
    timetillbright = 0*u.hour

    nightstart = times24[sunaltazs24.alt<-18*u.deg][0]
    nighttimemask = np.array((sunaltazs24.alt<-18*u.deg))*1
    if (sun_altaz.alt > -18*u.deg):
        nightend = times24[(np.roll(nighttimemask, 1) - nighttimemask) != 0][1]
        timetilldark = (nightstart-t)
        timetilldark.format = 'sec'
        LST = nightstart.sidereal_time('mean').deg
        HETphi = ((hetpupil[:,1]+LST)%360)*np.pi/180
        newpix = hp.ang2pix(nside, HETtheta, HETphi)
    else:
        nightend = times24[(np.roll(nighttimemask, 1) - nighttimemask) != 0][0]
        timetillbright = (nightend-t)
        timetillbright.format = 'sec'
        
    """
    timetill90 = 0
    mplot_combined, p90i = helper_function(observatory, t, m)
    
    for future_time in times24:
        mplot, _ = helper_function(observatory, future_time, m)
        if plot_timestamps:
            make_visibility_figure(savedir, observatory, mplot, future_time)
        mplot_combined[mplot == 0.5] = 0.5

    if plot:
        make_visibility_figure(savedir, observatory, mplot_combined, t, combined=True)
        
    visible_mask = np.where(mplot_combined == 0.5)[0]
    frac_visible = len(visible_mask)/len(mplot[p90i])
    m[mplot_combined != 0.5] = 0.0 # for speedup, also mask out every pixel not in 90% region at all    

    return timetill90, frac_visible


def main():
        skymap, header = hp.read_map(sys.argv[1],
                                     h=True, verbose=False)
        header = {'GraceID': 'TEST'}
        time = astropy.time.Time.now()
        timetill90, frac_visible = prob_observable(skymap, header, time, plot=True)
        print(frac_visible)
        return frac_visible
    
if __name__=="__main__":
    main()
