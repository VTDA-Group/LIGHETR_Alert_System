import healpy as hp
import astropy.coordinates
import astropy.time
import numpy as np
from astropy.coordinates import Angle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
from astropy.time import TimeDelta

maxhetdec = 74
minhetdec = -12


#maxhetdec = 90
#minhetdec = -25
minhetdec_rad = (90-maxhetdec)*np.pi/180
maxhetdec_rad = (90-minhetdec)*np.pi/180


HET_loc = (-104.01472,30.6814,2025)

def make_visibility_figure(savedir, time, m, plot = True):
    """
    This function makes the visibility figure at a given time.
    """
    
    # Determine resolution of sky map
    mplot = np.copy(m)
    npix = len(m)
    nside = hp.npix2nside(npix)
    
    nside_HET = nside
    
    hetpupil = np.loadtxt('hetpix.dat')
    hetfullpix = hp.query_strip(nside_HET, minhetdec_rad, \
                            maxhetdec_rad)
    observatory = astropy.coordinates.EarthLocation(lat=HET_loc[1]*u.deg, lon=HET_loc[0]*u.deg, height=HET_loc[2]*u.m)

    # Find pixels of HET pupil in this time
    t = astropy.time.Time(time,scale='utc',location=HET_loc)
    LST = t.sidereal_time('mean').deg
    HETphi = ((hetpupil[:,1]+LST)%360)*np.pi/180
    HETtheta = (90-hetpupil[:,2])*np.pi/180
    newpix = hp.ang2pix(nside_HET, HETtheta, HETphi)
    
    newpixp = newpix
    
    
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
    p90i = np.where(cls <= 90)
    
    
    p90i_HET_visible = None
    
    if sun_altaz.alt < -18*u.deg:
        p90i_HET_visible = np.intersect1d(p90i,newpix)
        #p90i_HET_visible = np.intersect1d(p90i,hetfullpix)

    
    if plot:
        if not os.path.exists(savedir+"HET_visibilities/"):
            os.mkdir(savedir+"HET_visibilities/")
        
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
        
        
        mplot[newpixp] = 0.2 #purple, color of HET pixels
        mplot[p90i] = 0.4 #blue, location of 90% region contours
        if sun_altaz.alt < -18*u.deg:
            mplot[p90i_HET_visible] = 0.5 #location of 90% region contours that are visible to HET while the sun is down
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
        plt.savefig(savedir+"HET_visibilities/"+'HET_visibility_figure_'+str(round(time_hr, 2))+'_.pdf', bbox_inches = 'tight')
    
    
    return p90i_HET_visible

    


def prob_observable(m, header, time, savedir, plot = True):

    
    # Determine resolution of sky map
    mplot = np.copy(m)
    npix = len(m)
    nside = hp.npix2nside(npix)

    # Geodetic coordinates of MacDonald Obs
    hetpupil = np.loadtxt('hetpix.dat')
    hetfullpix = hp.query_strip(nside, minhetdec_rad, \
                            maxhetdec_rad)

    observatory = astropy.coordinates.EarthLocation(
        lat=HET_loc[1]*u.deg, lon=HET_loc[0]*u.deg, height=HET_loc[2]*u.m)

    # Find pixels of HET pupil in this time
    t = astropy.time.Time(time,scale='utc',location=HET_loc)
    LST = t.sidereal_time('mean').deg
    HETphi = ((hetpupil[:,1]+LST)%360)*np.pi/180
    HETtheta = (90-hetpupil[:,2])*np.pi/180
    newpix = hp.ang2pix(nside, HETtheta, HETphi)
    
    newpixp = newpix

    
    
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
    p90i = np.where(cls <= 90)

    # How likely is it that the (true, unknown) location of the source
    # is within the area that is visible, in this time and within 24 hours? 
    # Demand that it falls in the HETDEX pupil, that the sun is at least 18 
    # degrees below the horizon and that the airmass (secant of zenith angle 
    # approximation) is at most 2.5.
    
    
    delta_time = np.linspace(0, 24, 1000)*u.hour
    times24 = t + delta_time
    frames24 = astropy.coordinates.AltAz(obstime = times24, location=observatory)
    sunaltazs24 = astropy.coordinates.get_sun(times24).transform_to(frames24)
    timetilldark = 0*u.hour
    timetillbright = 0*u.hour

    nightstart = times24[sunaltazs24.alt<-18*u.deg][0]
    nighttimemask = np.array((sunaltazs24.alt<-18*u.deg))*1
    if (sun_altaz.alt > -18*u.deg):
        nightend = times24[(np.roll(nighttimemask, 1) - nighttimemask) != 0][1]
    else:
        nightend = times24[(np.roll(nighttimemask, 1) - nighttimemask) != 0][0]
    nightime = nightend - nightstart
    nightime.format = 'sec'
    #Moving to start of the night if in daytime
    if (sun_altaz.alt > -18*u.deg):
        timetilldark = (nightstart-t)
        timetilldark.format = 'sec'
        LST = nightstart.sidereal_time('mean').deg
        HETphi = ((hetpupil[:,1]+LST)%360)*np.pi/180
        newpix = hp.ang2pix(nside, HETtheta, HETphi)
    else:
        timetillbright = (nightend-t)
        timetillbright.format = 'sec'


    theta90, phi90 = hp.pix2ang(nside, p90i)
    #mask skymap pixels by hetdex accesible region
    theta90HETi = (theta90 > minhetdec_rad)*(theta90 < maxhetdec_rad)
    print(theta90HETi.sum(),theta90.min(),theta90.max())
    theta90HET = theta90[theta90HETi]
    phi90HET = phi90[theta90HETi]
    timetill90 = 0

    '''
    if len(np.intersect1d(p90i,newpix)) == 0: #if the region doesn't intersect HET now

        #if the region doesn't intersect HET at all
        if len(np.intersect1d(p90i,hetfullpix)) == 0:
            return  -99, 0, 0
        hetedge = np.loadtxt('hetedge.dat')
        hetedgef = lambda x: np.interp(x,hetedge[:,0],hetedge[:,1])
        y = theta90HET*180/np.pi #DEC SKYMAP
        x = (phi90HET*180/np.pi-LST)%360 #RA SKYMAP ZEROING HET PUPIL
        wsecs = np.min(x - hetedgef(y))*3600*12/180
        if wsecs < 0:
            #it's inside the pupil...
            hetedgef2 = lambda x: np.interp(x,hetedge[:,0],hetedge[:,2])
            y = (hetedgef2(y)+180)%360 -180
            x = (x+180)%360-180
            wsecs = np.min(x - y)*3600*12/180

        if timetilldark == 0:
            if wsecs > timetillbright.value:
                return  -99, 0, 0
        else:
            if wsecs > nightime.value:
                return  -99, 0, 0
        timetill90 = (wsecs+timetilldark.value)/3600
    elif timetilldark.value > 0:
        timetill90 = timetilldark.value/3600
    '''

    
    
    
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
    
    
    mplot[newpixp] = 0.2 #black, color of HET pixels
    mplot[p90i] = 0.4 #blue, location of 90% region contours
    #for i in range(24*6):
    delta_time = np.linspace(0, 24, 24)*u.hour
    for dt in delta_time:
        future_time = time + dt
    
        p90i_HET_visible = make_visibility_figure(savedir=savedir, m = m, time = future_time, plot = False)
        
        if p90i_HET_visible is not None:
            #print("m of the visible region is to be: "+str(m[p90i_HET_visible])+", mplot is being set to 0.5")
            mplot[p90i_HET_visible] = 0.5
    
    #find mask of where mplot still is 0.4, this is the part of the 90% region that is never visible to HET at nighttime, because it was never set to 0.5
    
    
    never_visible_mask = np.where(mplot == 0.4)[0]
    mplot[never_visible_mask] = 0.7
    visible_mask = np.where(mplot == 0.5)[0]
    
    frac_visible = len(visible_mask)/(len(never_visible_mask) + len(visible_mask))
    
    
    #m[never_visible_mask] = 0.0
    m[mplot != 0.5] = 0.0 # for speedup, also mask out every pixel not in 90% region at all

    
    
        

    
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
        
    plt.savefig(savedir+"HET_Full_Visibility.pdf")
    

    

    

    # Done!
    
    
    
    return timetill90, m, frac_visible

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
