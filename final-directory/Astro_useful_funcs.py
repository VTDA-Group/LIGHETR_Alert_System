import math
import random
import warnings

G_SI = 6.674E-11 #Universal gravitational constant G in SI units N*m^2/kg^2
G_cgs = 6.674E-8 #Universal gravitational constant G in cgs units dyne*cm^2/g^2
c_ms = 3.0E8 #speed of light in meters/sec

H_0_km_s_Mpc = 72.0 #hubble constant in units of km/s/Mpc
H_0_hz = 2.269E-18 #hubble constant in units of 1/sec, this is useful for SI unit calculations
H_0 = H_0_km_s_Mpc
H_0_SI = H_0_hz

h = H_0_km_s_Mpc/100 #dimensionless h constant for Hubble constant
rho_crit_SI = H_0_SI**2*3/(8*np.pi*G_SI) #critical density of the universe in SI units
rho_crit_cgs = 1.9E-29*h #in cgs units of g/cm3

Omega_M = 0.3147 #mass density of the universe
Omega_Lambda = 0.6847 #dark energy density of universe
Omega_k = 0.001
t_0_Gyr = 13.797 #age of universe in Gyr
t_0 = t_0_Gyr

def get_plotting_marker(band, instrument = None, telescope = None):
    '''inputs: band-the name of the band that we want to assign a color and marker to, instrument-the instrument on which the observation was done (defaults to none-instrument does not need to be provided), telescope-the telescope on which the observation was done (defaults to none-instrument does not need to be provided)
       outputs: color-color of marker to be plotted, marker-marker in matplotlib that needs to be plotted
       This function takes a band in which a photmetric observation was done and returns the color and marker in which we want to plot this observation on a light curve. This will allow for standardization for LC plotting over all of my figures. Returns as a tuple: (color, marker)'''
       
    
    key = band
    
    colors_to_attempt = {'g':'limegreen', 'gp':'limegreen', 'up':'darkslategrey', 'ip':'indianred', 'rp':'crimson', 'G': 'green', 'r':'crimson', 'R': 'red', 'b':'royalblue', 'B':'blue', 'cyan':'cyan', 'orange':'orange', 'u':'darkslategrey', 'U':'darkslategrey', 'UVM2':'violet', 'UVW2':'slateblue', 'UVW1': 'midnightblue', 'i':'indianred', 'z':'brown', 'v':'orangered', 'V':'palevioletred', 'generic':'black', 'w':'navy'}
    colors_assigned = {}
    used_up_colors = []

    detector_dictionary = {"UVM2":{"telescope":"Swift", "instrument":"UVOT"}, "UVW1":{"telescope":"Swift", "instrument":"UVOT"}, "UVW2":{"telescope":"Swift", "instrument":"UVOT"}}
    
    instrument = None
    telescope = None

    for detector_key in detector_dictionary.keys():
        for color_key in colors_to_attempt.keys():
            if color_key+"-" in detector_key or color_key == detector_key:
                instrument = detector_dictionary[detector_key]["instrument"]
                telescope = detector_dictionary[detector_key]["telescope"]
    
    telescope_groups = []
    symbols_to_attempt = {'Swift': '^', 'ZTF': 'o', 'P48':'o', 'ATLAS': 'D', 'Siding_Spring_1m':'s', 'Thacher':'P', 'generic':'.'}
    symbols_assigned = {}
    used_up_symbols = []
    telescope_names = {'Siding_Spring_1m': 'Las Cumbres', 'P48': 'ZTF'}
    
    curr_col = ''
    curr_symbol = ''
    for color_key in colors_to_attempt.keys():
        if color_key+"-" in key or color_key == key:
            colors_assigned[key] = colors_to_attempt[color_key]
            used_up_colors.append(colors_to_attempt[color_key])
            curr_col = color_key
            
    if key not in colors_assigned.keys():
        print("Key not in colors assgined: "+str(key))
        colors_assigned[key] = colors_to_attempt['generic']
        
    
            
    for symbol_key in symbols_to_attempt.keys():
        if symbol_key in key or symbol_key == key:
            symbols_assigned[key] = symbols_to_attempt[symbol_key]
            used_up_symbols.append(symbols_to_attempt[symbol_key])
            curr_symbol = symbol_key
            
    if key not in symbols_assigned.keys() and telescope is not None and telescope in symbols_to_attempt.keys():
        symbols_assigned[key] = symbols_to_attempt[telescope]
        used_up_symbols.append(symbols_to_attempt[telescope])
        curr_symbol = telescope
        
    elif key not in symbols_assigned.keys() and instrument is not None and instrument in symbols_to_attempt.keys():
        symbols_assigned[key] = symbols_to_attempt[instrument]
        used_up_symbols.append(symbols_to_attempt[instrument])
        curr_symbol = instrument

    if key not in symbols_assigned.keys():
        print("Key not in symbols assgined: "+str(key))
        symbols_assigned[band] = symbols_assigned['generic']
        
    print("Colors assgined: "+str(colors_assigned))
    print("symbols assgined: "+str(symbols_assigned))
    
    return colors_assigned[band], symbols_assigned[band]


def find_nearest(wv,value):
    '''inputs: wv-list of data to look in, value-the value we are looking for in wv
       outputs: idx-index in wv that is closest to value
       this finds the index of the closest element to value in list wv'''
    idx= (np.abs(np.asarray(wv) - value)).argmin()
    return idx
    
def gaussian(x, mu, sig, A = -np.inf):
    '''inputs: x-list or array of independent variables over which to build a gaussian, mu-center of gaussian, sig-standard deviation of gaussian, A-vertical scaling constant of gaussian (default to 1.0)
       outputs: gaussian-array of dependent variable defining gaussian with input parameters
       This function passes out a gaussian with the parameters that are passed in'''
    
    
    if A == -np.inf:
        GG = (1/(sig*np.sqrt(2*np.pi))) * np.asarray(np.exp(-np.power(np.asarray(x) - mu, 2.) / (2 * np.power(sig, 2.))))
    else:
        GG = np.asarray(A*np.exp(-np.power(np.asarray(x) - mu, 2.) / (2 * np.power(sig, 2.))))
    return GG
    
def quadratic_std(x, a, b, c):
    '''inputs: x-list or array of independent variables over which to build a quadratic in standard form (of the form ax^2+bx+c), a-coefficient in front of x^2, b-coefficient in front of x, c-extra constant added to quadratic
       outputs: quad-array of dependent variable defining gaussian with input parameters
       This function passes out a quadratic constructed in standard form ax^2+bx+c with the parameters that are passed in'''
    return np.asarray(a*np.asarray(x)**2+b*np.asarray(x)+c)
    
def quadratic_vertex(x, a, h, k):
    '''inputs: x-list or array of independent variables over which to build a quadratic in vertex form (of the form a(x-h)^2+k), a-coefficient in front of x^2, h-x coordinate of vertex, k-y coordinate of vertex
       outputs: quad-array of dependent variable defining gaussian with input parameters
       This function passes out a quadratic constructed in standard form a*(x-h)^2+k with the parameters that are passed in'''
    return np.asarray(a*(np.asarray(x)-h)**2+k)
    
def quadratic_zeros(x, a, x1, x2):
    '''inputs: x-list or array of independent variables over which to build a quadratic in zero point form (a*(x-x1)(x-x2)), a-coefficient in front of x^2, x1-x coordinate first zero, x2-x coordinate of second zero
       outputs: quad-array of dependent variable defining gaussian with input parameters
       This function passes out a quadratic constructed in standard form a*(x-x1)(x-x2) with the parameters that are passed in'''
    return np.asarray(a*(np.asarray(x)-x1)*(np.asarray(x)-x2))
    
def double_gaussian(x, mu1, sig1, mu2, sig2, A1 = 1.0, A2 = 1.0):
    '''inputs: x-list or array of independent variables over which to build a gaussian, mu1-center of gaussian1, sig1-standard deviation of gaussian1, A1-vertical scaling constant of gaussian1 (default to 1.0), mu2-center of gaussian2, sig2-standard deviation of gaussian2, A2-vertical scaling constant of gaussian2 (default to 1.0),
       outputs: gaussian-array of dependent variable defining gaussian with input parameters
       This function passes out a gaussian with the parameters that are passed in'''
    return np.asarray(A1*np.exp(-np.power(np.asarray(x) - mu1, 2.) / (2 * np.power(sig1, 2.)))) + np.asarray(A2*np.exp(-np.power(np.asarray(x) - mu2, 2.) / (2 * np.power(sig2, 2.))))
    
def noisy_gaussian(x, mu, sig, noise_strength, A = 1.0):
    '''inputs: x-list or array of independent variables over which to build a gaussian, mu-center of gaussian, sig-standard deviation of gaussian, A-vertical scaling constant of gaussian (default to 1.0), noise_strength-strength of noise parameter in the gaussian that is created
       outputs: gaussian-array of dependent variable defining gaussian with input parameters
       This function passes out a randomized gaussian with the parameters that are passed in, with an extra noise term controlled by noise_strength'''
    
    '''
    gg = noise_strength*np.random.random_sample(size = len(x))
    hh = (1+gg)*np.asarray(A*np.exp(-np.power(np.asarray(x) - mu, 2.) / (2 * np.power(sig, 2.))))
    return gg*hh'''
    
    g = np.asarray(A*np.exp(-np.power(np.asarray(x) - mu, 2.) / (2 * np.power(sig, 2.))))
    tor = []
    for i in g:
        tor.append(i*(1+random.random()*noise_strength))
    return np.asarray(tor)
    
def noisy_double_gaussian(x, mu1, sig1, mu2, sig2, noise_strength, A1 = 1.0, A2 = 1.0):
    '''inputs: x-list or array of independent variables over which to build a gaussian, mu1-center of gaussian1, sig1-standard deviation of gaussian1, A1-vertical scaling constant of gaussian1 (default to 1.0), mu2-center of gaussian2, sig2-standard deviation of gaussian2, A2-vertical scaling constant of gaussian2 (default to 1.0), noise_strength-strength of noise parameter in the gaussian that is created
       outputs: gaussian-array of dependent variable defining gaussian with input parameters
       This function passes out a randomized gaussian with the parameters that are passed in, with an extra noise term controlled by noise_strength'''
    
    '''
    gg = noise_strength*np.random.random_sample(size = len(x))
    hh = (1+gg)*np.asarray(A*np.exp(-np.power(np.asarray(x) - mu, 2.) / (2 * np.power(sig, 2.))))
    return gg*hh'''
    
    g = np.asarray(A1*np.exp(-np.power(np.asarray(x) - mu1, 2.) / (2 * np.power(sig1, 2.)))) + np.asarray(A2*np.exp(-np.power(np.asarray(x) - mu2, 2.) / (2 * np.power(sig2, 2.))))
    tor = []
    for i in g:
        tor.append(i*(1+random.random()*noise_strength))
    return np.asarray(tor)
    
    
def convolve_analytical(seed_func, kernel_func, x_half_width, x_prec):
    '''inputs: seed_func-a python function that gives the mathematical seed function that we want to convolve, kernel_func-a python function that gives the mathematical kernel function that we want to convolve our seed function with, x_half_width-the half width over the independent variable that we want to convolve over, x_prec-precision of the integral
       outputs: Xs-an array of the independent variable that is twice as long as the region over which the kernel function and seed function are defined. This is the region over which the convolved function will be defined, xs-array of the independent variable over which the kernel function and seed function are defined, convolved-array of final convolved function
       This function takes the convolution of seed_func and kernel func. Note that seed_fun and kernel func must be python functions that take a single variable as a parameter.'''
    
    
    Xs = np.asarray(np.arange(-2.*x_half_width, 2.*x_half_width, x_prec))
    xs = np.asarray(np.arange(-1.*x_half_width, x_half_width,  x_prec))
    
    convolved = []
    for i in range(len(Xs)):
        X = Xs[i]
        local_sum = 0.0
        
        for j in range(len(xs)):
            x = xs[j]
            if j < len(xs)-1:
                width = xs[j+1]-xs[j]
            else:
                width = xs[-1]-xs[-2]
                
            local_sum += kernel_func(x)*seed_func(X-x)*width
        convolved.append(local_sum)
    
    return {"Xs": Xs, "xs": xs, "convolved": np.asarray(convolved)}
    
    
    
    
def convolve_empirical(x, y, kernel_func, width = -np.inf, x_prec = -np.inf):

    '''inputs: x-an array of x values (independent variables) over which the seed function is empirically defined, y-an array of y values (dependent values) that hold the actual seed function values, kernel_func-a python function that gives the mathematical kernel function that we want to convolve our seed function with, x_prec-precision of the integral (defaults to -infinity, which the function interpret as estimating the accuracy of the precision from the precision of the x array)
       outputs: Xs-an array of the independent variable that is twice as long as the region over which the kernel function and seed function are defined. This is the region over which the convolved function will be defined, xs-array of the independent variable over which the kernel function and seed function are defined, kernel_arr-array of kernel values defined over xs, convolved-array of final convolved function, norm_convolved-array of final convolved function divided by max value of convolved,
       This function takes the convolution of seed_func and kernel func. Note that seed_fun is described by an array for x and y and kernel func must be a python function that take a single variable as a parameter.'''
    
    midpoint = np.median(x)
    x-=np.median(x)
    
    if width == -np.inf:
        width = 2.
    
    xs = x
    ys = y
    
    
    med = np.median(x)
    lw = med-np.min(x)
    uw = np.max(x)-med
    
    if x_prec == -np.inf:
        x_prec = (np.max(x)-np.min(x))/len(x)
        
    l = list(np.arange(np.min(x)-(width-1)*lw, np.min(x), x_prec))
    Xs = l
    Ys = [0. for i in l]
    
    for i in range(len(x)):
        Xs.append(x[i])
        Ys.append(y[i])
        
    u = list(np.arange(np.max(x), np.max(x)+(width-1)*uw, x_prec))
    for i in range(len(u)):
    
        Xs.append(u[i])
        Ys.append(0.)
    
    Xs = np.asarray(Xs)
    Ys = np.asarray(Ys)
    
    
    convolved = []
    kernel_arr = []
    for i in range(len(Xs)):
        X = Xs[i]
        local_sum = 0.0
        add_to_kernel = False
        if len(kernel_arr) == 0:
            add_to_kernel = True
        for j in range(len(xs)):

            x = xs[j]
            if j < len(xs)-1:
                width = xs[j+1]-xs[j]
            else:
                width = xs[-1]-xs[-2]
                
            idx = find_nearest(X-x, Xs)
            #print("")
            if add_to_kernel:
                kernel_arr.append(kernel_func(x))
            
            local_sum += kernel_func(x)*Ys[idx]*width
        convolved.append(local_sum)
    
    Xs += midpoint
    xs += midpoint
    return {"Xs": Xs, "xs": xs, "kernel_arr": np.asarray(kernel_arr), "convolved": np.asarray(convolved), "norm_convolved": np.asarray(convolved)/np.max(convolved)}
    
    
def correct_cosmological_redshift(wavelengths, fluxes, z):
    '''inputs: wavelengths-list of wavelengths for the spectrum that is passed in, fluxes-corresponding list of fluxes for the spectrum, z-cosmological redshift of the obeject emitting this spectrum
       outputs: wavelengths_prog_frame-list of wavelengths of the spectrum in progenitor frame, fluxes_prog_frame-list of corresponding fluxes of the spectrum in progenitor frame
       This function takes an observed spectrum and converts it into a progenitor frame spectrum, given a redshift z'''
    fluxes_prog_frame = [i*(1+z) for i in fluxes]
    wavelengths_prog_frame = [i/(1+z) for i in wavelengths]
    
    return wavelengths_prog_frame, fluxes_prog_frame

def power_from_obs_flux(obs_flux, distance):
    '''inputs: obs_flux-the observed flux observation of an object, distance-distance to the object
       outputs: power-total power output of object
       This function gives the total power output of an object in units of [obs_flux]*[distance]^2--as in, whatever units are passed in are the units that come out. Note that this just gives the total power output, not the luminosity at some distance or some measure of apparent magnitude or absolute magnitude. Also note that this function inherently assumes this object is an isotropic emitter.'''
    
    power = obs_flux*4*np.pi*distance**2.0
    return power
    

def distance_modulus_from_dist(d_cm = None, d_m = None, d_Km = None, d_AU = None, d_ly = None, d_pc = None,  d_Kpc = None, d_Mpc = None, d_Gpc = None):
    '''inputs: every possible length
       outputs: distance modulus-the corresponding distance modulus to the object, given its distance. This is taken to be its lumninosity distance
       This function gives distance modulus to an object, from its distance'''
    dist_pc = length_conversions(d_cm = d_cm, d_m = d_m, d_Km = d_Km, d_AU = d_AU, d_ly = d_ly, d_pc = d_pc, d_Kpc = d_Kpc, d_Mpc = d_Mpc, d_Gpc = d_Gpc)["pc"]
    dist_mod = 5*np.log10(dist_pc/10)
    
    return dist_mod


def get_absolute_mag(app_mag, dist_parsec):
	'''inputs: app_mag-apparent magnitude of an object, dist_parsec-distance to the object in parsecs
	   outputs: absolute_mag-absolute magnitude of the object, dist_mod-distance modulus of the object
	   this function uses the equation m-M = 5log(d/10 parsec) to get the distance modulus and the absolute magnitude of an object, given its apparent magnitude and distance'''
	
	dist_mod = 5*np.log10(dist_parsec/10)
	absolute_mag = -1.0*dist_mod + app_mag
	
	return {"absolute_mag": absolute_mag, "dist_mod":dist_mod}
 
def get_apparent_mag(abs_mag, dist_parsec):
    '''inputs: abs_mag-absolute magnitude of an object, dist_parsec-distance to the object in parsecs
       outputs: app_mag-apparent magnitude of the object, dist_mod-distance modulus of the object
       this function uses the equation m-M = 5log(d/10 parsec) to get the distance modulus and the apparent magnitude of an object, given its absolute magnitude and distance in parsec'''
    
    dist_mod = 5*np.log10(dist_parsec/10)
    app_mag = abs_mag + dist_mod
    
    return {"apparent_mag": app_mag, "dist_mod":dist_mod}


def fnu_from_flambda(flambda, lam, c = 3.0E18):
	'''inputs: flambda-f_lambda that we want to get the f_nu of, lam-wavelength at which this is happening (angstroms is default), c-optional speed of light (this is useful to specify units of length, if something other than angstroms is needed, specify speed of light in these length units)
	   outputs: fnu-f_nu from converting flambda above, units are returned in same units of flambda as passed in, multiplied by length units passed in divided by seconds
	   this converts flambda into fnu for a given measurement. The default length units are angstroms and can be changed by changing the units of c, default times units are seconds, and cannot be changed. The function returns fnu in the units of flambda multiplied by the length units specified and divided 
by seconds'''

	fnu = flambda*lam**2/c
	return fnu


def flambda_from_fnu(fnu, lam, c = 3.0E18):
    '''inputs: fnu-f_nu that we want to get the f_lambda of, lam-wavelength at which this is happening (angstroms is default), c-optional speed of light (this is useful to specify units of length, if something other than angstroms is needed, specify speed of light in these length units)
       outputs: flambda-f_lambda from converting nu above, units are returned in same units of flambda as passed in, multiplied by seconds divided by length units passed in
       this converts flambda into fnu for a given measurement. The default length units are angstroms and can be changed by changing the units of c, default times units are seconds, and cannot be changed. The function returns flambda in the units of fnu multiplied by seconds divided by the length units passed in'''

    flambda = fnu*c/lam**2
    return flambda


def flux_nu_jansky_from_flux_nu_cgs(flux_nu_cgs):
    '''inputs: flux_nu_cgs-flux in units of erg/cm^2/s/Hz
       outputs: flux_nu_jansky-flux in units of janskies
    '''
    return flux_nu_cgs * 10.**23
    
def flux_nu_cgs_from_flux_nu_jansky(flux_nu_jansky):
    '''input: flux_nu_jansky-flux in units of janskies
       output: flux_nu_cgs-flux in units of erg/cm^2/s/Hz
    '''
    return flux_nu_jansky * 10.**-23


def get_AB_mag_from_flux_nu(flux_nu_cgs = -np.pi, flux_nu_jansky = -np.pi):
    '''inputs: flux_nu_cgs-flux in units of erg/cm^2/s/Hz, flux_nu_jansky-flux in units of janskies
       outputs: AB_mag-AB magnitude measure of magnitude from flux that is passed in
       This function gets the AB magnitude of a f_nu measurement. Either flux_nu_cgs or flux_nu_jansky needs to be passed in (or both). If neither is passed in, an error is raised
    '''
    if flux_nu_cgs == -np.pi and flux_nu_jansky == -np.pi:
        print("Error: flux not passed in")
        return np.pi
    elif flux_nu_cgs != -np.pi and flux_nu_jansky != -np.pi and flux_nu_jansky_from_flux_nu_cgs(flux_nu_cgs) != flux_nu_jansky:
        print("Error: fluxes not consistent")
        return np.pi
    elif flux_nu_cgs != -np.pi:
        flux_nu_jansky = flux_nu_jansky_from_flux_nu_cgs(flux_nu_cgs)
    
    AB_mag = -2.5*np.log10(flux_nu_jansky)+8.9
    
    return AB_mag
    
def get_flux_nu_from_AB_mag(AB_mag):
    '''inputs: AB_mag-AB magnitude measure of magnitude
       outputs: flux_nu_cgs-flux in units of erg/cm^2/s/Hz, flux_nu_jansky-flux in units of janskies
       This function gets the AB magnitude of a f_nu measurement. Either flux_nu_cgs or flux_nu_jansky needs to be passed in (or both). If neither is passed in, an error is raised
    '''
    
    flux_nu_jansky = 10.0**((AB_mag-8.9)/-2.5)
    flux_nu_cgs = flux_nu_cgs_from_flux_nu_jansky(flux_nu_jansky)
    
    return {"flux_nu_jansky": flux_nu_jansky, "flux_nu_cgs": flux_nu_cgs}
    
def get_AB_mag_from_flux_lam(lam, flux_lam_cgs = -np.pi):
    '''inputs: lam-wavelength of this f_lam measurement in Angstroms, flux_lam_cgs-flux in units of erg/cm^2/s/angstrom
       outputs: AB_mag-AB magnitude measure of magnitude from flux that is passed in
       This function gets the AB magnitude of a f_lam measurement using the lam that is passed in. Note that AB_mag is defined natively in f_nu space, so a wavelength must be given to properly convert to f_nu space. This turns f_lam into f_nu_jansky, then gets AB-mag from that
    '''
    if flux_lam_cgs == -np.pi:
        print("Error: flux not passed in")
        return np.pi

    flux_nu_cgs = fnu_from_flambda(flambda = flux_lam_cgs, lam = lam)
    flux_nu_jansky = flux_nu_jansky_from_flux_nu_cgs(flux_nu_cgs = flux_nu_cgs)
    AB_mag = get_AB_mag_from_flux_nu(flux_nu_jansky = flux_nu_jansky)
    
    return AB_mag
        
def get_tangential_vel(proper_motion, dist_parsec):
    '''inputs: proper_motion-proper motion of object in arsec/yr, dist_parsec-distance to object in parsec
       outputs: tang_vel-tangential velocity of object in km/sec
       this gets the tangential velocity of an object in km/sec, given its proper motion in arcsec/yr and distance in parsecs
    '''
    tang_vel = 4.74*proper_motion*dist_parsec
    return tang_vel
   
def get_dist_tang_vel_method(radial_vel, lam_rad, proper_motion):
    '''inputs: radial_vel-radial velocity of object in km/sec, lam_rad between direction of radial velocity and direction of overal motion in radians, proper_motion-proper motion of object in arcsec/yr
       outputs: dist-distance to object in parsec
       this gets the distance to an object using the moving cluster method, given the radial velocity, angle between radial velocity and overall convergence direction, and proper motion
    '''
    
    dist = radial_vel*np.tan(lam_rad)/(4.74*proper_motion)
    
    return dist

def radian_from_degree(angle_degree):
    '''inputs: angle_degree-an angle measured in degrees
       outputs: angle_rad-angle measured in radians
       gets radians from degrees of an angle
    '''
    return np.pi*angle_degree/180
    
def deg_min_sec_from_degree(degree):
    '''inputs: degree-angle measured in degrees
       outputs: deg_min_sec-angle measured in units of degrees, minutes, seconds
       This function converts an angle in units of degrees to degrees minutes and seconds'''
       
    deg = int(degree)
    
    deg_frac = degree-deg
    minutes = deg_frac*60.
    
    minutes_frac = minutes-int(minutes)
    secs = minutes_frac*60.
    
    if deg < 0:
        minutes *= -1.0
        secs *= -1.0
    
    return str(deg)+":"+str(int(minutes))+":"+str(round(secs, 4))
    
    
    
    
def rad_from_degree_min_sec(angle_degree_min_sec):
    '''inputs: angle_degree_min_sec-angle measured in units of degrees, arcminutes, and arcseconds
       outputs: angle_rad-corresponding angle in radians
       This function gets radians from an angle given in units of degreesArcminArcsec'''
       
    GG = angle_degree_min_sec.split(":")
    deg_floor = float(GG[0])
    arc_min = float(GG[1])
    arc_sec = float(GG[2])
    
    if deg_floor < 0:
        arc_min *= -1.0
        arc_sec *= -1.0

    angle_deg = deg_floor + arc_min/60.+arc_sec/3600.
    
    return radian_from_degree(angle_deg)

def energy_from_wavelength(lam_m):
    '''inputs: lam_m-wavelength in METERS that we want the energy of
       outputs: lam_nm-wavelength of light in nanometers, lam_m-wavelength of light in meters, lam_a-wavelength of light in angstroms, energy_erg-energy of this light in ergs, energy_ev-energy of this light in ev, energy_joule-energy of this light in joules, freq-frequency of this light in hz, correct_color-string representing what band this would be (like 'radio', 'visible', 'gamma' or whatever), t_eff-effective temperature of this light, assume h*frequency = boltzmann constant*t_eff
       This function takes a wavelength in meters and gives back all the info of this light-including the energy in various units, frequency, effective temperature, and where on the EM spectrum it lies
    '''
    c_m = 3.0E8 #speed of light in meters/sec
    freq = c_m/lam_m #frequency of light in hz
    
    h = 6.626E-34 #planck constant in SI units
    energy_joule = h*freq #energy of light in Joules
    
    energy_ev = energy_joule*6.242E18 #energy of light in ev
    energy_erg = energy_joule*1.0E7 #energy of light in ergs
    
    kB_SI = 1.381E-23 #boltzmann constant in SI units
    t_eff = energy_joule/kB_SI #effective temperature of this light in Kelvin
    
    spectrum_wavelength_bounds = {"Gamma": [0, 1.0E-12], "Hard X-Ray": [1.0E-12, 1.0E-10], "Soft X-Ray": [1.0E-10, 3E-9], "Extreme Ultraviolet": [3E-9, 30E-9], "Near Ultraviolet": [30E-9, 380E-9], "Violet": [380E-9, 450E-9], "Blue":[450E-9, 485E-9], "Cyan":[485E-9, 500E-9], "Green":[500E-9, 565E-9], "Yellow":[565E-9, 590E-9], "Orange":[590E-9, 625E-9], "Red":[625E-9, 750E-9], "Near Infrared":[750E-9, 3E-6], "Mid Infrared":[3E-6,15E-6], "Far Infrared":[15E-6, 1E-3], "Microwave (EHF)":[1E-3, 10E-3], "Microwave (SHF)":[1E-2, 10E-2], "Microwave (UHF)":[1E-1, 1], "Radio (VHF)":[1, 10], "Radio (HF)":[10, 100], "Radio (MF)":[100, 1000], "Radio (LF)":[1E3, 10E3], "Radio (VLF)":[10E3, 100E3], "Radio (ULF)":[100E3, 1000E3], "Radio (SLF)":[1.0E6, 1.0E7], "Radio (ELF)":[1.0E7, 1.0E8]}
    
    correct_color = 'Not found'
    for color in spectrum_wavelength_bounds.keys():
        color_bounds = spectrum_wavelength_bounds[color]
        
        if color_bounds[0] <= lam_m and lam_m < color_bounds[1]:
            correct_color = color
    

    return {"lam_nm":lam_m*1.0E9, "lam_m":lam_m, "lam_a":lam_m*1.0E10, "energy_erg":energy_erg, "energy_ev":energy_ev, "energy_joule":energy_joule, "freq": freq, "correct_color":correct_color, "t_eff":t_eff}

def wavelength_from_energy(energy_joule):
    '''inputs: energy_joule-energy of light in joules
       outputs: lam_nm-wavelength of light in nanometers, lam_m-wavelength of light in meters, lam_a-wavelength of light in angstroms, energy_erg-energy of this light in ergs, energy_ev-energy of this light in ev, energy_joule-energy of this light in joules, freq-frequency of this light in hz, correct_color-string representing what band this would be (like 'radio', 'visible', 'gamma' or whatever), t_eff-effective temperature of this light, assume h*frequency = boltzmann constant*t_eff
       This function takes an energy in joules and gives back all the info of this light-including the energy in various units, frequency, effective temperature, and where on the EM spectrum it lies
    '''
    h = 6.626E-34 #planck constant in SI units
    c = 3.0E8 #speed of light in SI units
    lam_m = c/(energy_joule/h)
    
    d = energy_from_wavelength(lam_m = lam_m)
    
    return d
    
def wavelength_from_temp(temp_k): #this function will be overloaded with different measures of energy, to turn into wavelength
    '''inputs: temp_k-temperature in Kelvin
       outputs: lam_nm-wavelength of light in nanometers, lam_m-wavelength of light in meters, lam_a-wavelength of light in angstroms, energy_erg-energy of this light in ergs, energy_ev-energy of this light in ev, energy_joule-energy of this light in joules, freq-frequency of this light in hz, correct_color-string representing what band this would be (like 'radio', 'visible', 'gamma' or whatever), t_eff-effective temperature of this light, assume h*frequency = boltzmann constant*t_eff
       This function takes a temperature in Kelvin and gives back all the info of this light-including the energy in various units, frequency, effective temperature, and where on the EM spectrum it lies
    '''
    kB_SI = 1.381E-23 #boltzmann constant in SI units
    d = wavelength_from_energy(energy_joule = kB_SI*temp_k)
    
    
    return d
    
    
def clean_number(num_to_clean, num_digits = 5):
    '''inputs: num_to_clean-the number to be cleaned and be made to look pretty, clean_arg-default to 5, if passed can set how many digits to have in the number
       outputs: cleaned_num-number that has been processed and will be returned
       this function takes in a number that may be very large (13678323450) or very small (0.00000000554) or have some extra small digit added to it (4534.00000001) and will return a pretty version of the same number as a float'''
    try:
        num_to_clean = float(num_to_clean)
    except ValueError as e:
        return num_to_clean
    
    if num_to_clean > 1:
        #looks like 457832948598 or maybe like 54783928457394.0000004 or like 548372945833.45387438532
        #either way, it needs to look like 5.6584E16 or whatever-specifically, with only clean_arg-1 digits after the decimal
        
        power = math.floor(np.log10(num_to_clean))
        C = round(num_to_clean/10.**power, num_digits-1)
        cleaned_num = float(str(C)+"E"+str(power))
    else:
        #looks like 0.0438298 or 0.33333333331 or whatever. This bit is easy, we just need to round the decimal, if the rounded version is just zero, we'll have to mess around a little
        cleaned_num = round(num_to_clean, num_digits-1)
    return cleaned_num
    
def str_clean_number(num_to_clean, num_digits = 5, max_power = 4):
    '''inputs: num_to_clean-the number to be cleaned and be made to look pretty, num_digits-default to 5, if passed can set how many digits to have in the number, max_power-default to 3, if passed can set maximum power to display without scientific notation
       outputs: cleaned_num-number that has been processed and will be returned
       this function takes in a number that may be very large (13678323450) or very small (0.00000000554) or have some extra small digit added to it (4534.00000001) and will return a pretty version of the same number as a string, so it can be displayed'''
    try:
        num_to_clean = float(num_to_clean)
    except ValueError as e:
        return num_to_clean
    
    
    if num_to_clean == 0.0:
        return "0.0"
    if num_to_clean > 1:
        #looks like 457832948598 or maybe like 54783928457394.0000004 or like 548372945833.45387438532
        #either way, it needs to look like 5.6584E16 or whatever-specifically, with only clean_arg-1 digits after the decimal
        
        power = math.floor(np.log10(num_to_clean))
        C = round(num_to_clean/10.**power, num_digits-1)
        if power <= max_power:
            cleaned_num = str(C*10.**power)
        else:
            cleaned_num = str(C)+"E"+str(power)
    else:
        #looks like 0.0438298 or 0.33333333331 or whatever. This bit is easy, we just need to round the decimal
        cleaned_num = str(round(num_to_clean, num_digits-1))
        
        
        if cleaned_num == "0.0" and num_to_clean != 0.0:
            cleaned_num = str(num_to_clean)[0:num_digits]
            power = math.floor(np.log10(num_to_clean))
            
            C = round(num_to_clean/10.**power, num_digits-1)
            
            cleaned_num = str(C)+"E"+str(power)
    return cleaned_num
    
#reddening, de-reddening
def get_reddening_extinction(R_lam, E_B_V):
    
    '''inputs: R_lam-the R value for reddening in A_lam = R_lam*E(B-V), E_B_V-the reddening value E(B-V)
       outputs: A_lam-the extinction in magnitudes using A_lam = R_lam*E(B-V) at a given wavelength lam
       This function just gives the number of magnitudes that are extinguished by reddening'''
    A_lam = R_lam*E_B_V
    
    return A_lam
    
def get_E_B_V(H_alpha_H_beta_flux_obs_ratio, H_alpha_H_beta_flux_true_ratio = 2.86, R_HBeta = 3.609, R_HAlpha = 2.535):
    '''inputs: H_alpha_H_beta_flux_obs_ratio-ratio of the flux from H-alpha to the flux from H-beta, H_alpha_H_beta_flux_true_ratio-the true value of this flux, this is a constant known to be 2.86, R_HBeta-the R value for H_Beta (default to 3.609, the avg for Milky Way), R_HAlpha-the R value for H_Alpha (default to 2.535, the avg for Milky Way)
       outputs: E_B_V-value of E(B-V) for the given flux ratio
       this function assumes Milky Way values for reddening of both H-Alpha and H-Beta and gives the corresponding E_B_V'''
       
    E_B_V = 2.5*(np.log10(H_alpha_H_beta_flux_obs_ratio)-np.log10(H_alpha_H_beta_flux_true_ratio))/(R_HBeta-R_HAlpha)
    
    return E_B_V

def deredden_flux(f_lam_obs, R_lam, H_alpha_H_beta_flux_obs_ratio, H_alpha_H_beta_flux_true_ratio = 2.86, R_HBeta = 3.609, R_HAlpha = 2.535):
    '''inputs: f_lam_obs-observed (reddened) flux at a given wavelength, H_alpha_H_beta_flux_obs_ratio-ratio of the flux from H-alpha to the flux from H-beta, H_alpha_H_beta_flux_true_ratio-the true value of this flux, this is a constant known to be 2.86, R_HBeta-the R value for H_Beta (default to 3.609, the avg for Milky Way), R_HAlpha-the R value for H_Alpha (default to 2.535, the avg for Milky Way)
       outputs: f_lam_true-true (dereddened) value of f_lam
       this function assumes Milky Way values for reddening of both H-Alpha and H-Beta and finds the de-reddened flux for a given observed flux'''
    
    E_B_V = get_E_B_V(H_alpha_H_beta_flux_obs_ratio = H_alpha_H_beta_flux_obs_ratio, H_alpha_H_beta_flux_true_ratio = 2.86, R_HBeta = 3.609, R_HAlpha = 2.535)
    
    
    f_lam_true = f_lam_obs*10.**(R_lam*E_B_V/2.5)
    return f_lam_true
    
       
#density diagnostics

def density_diagnostic(f_3729, f_3726):
    '''inputs: f_3729-flux in 3729 angstroms, f_3726-flux in 3726 angstroms
       outputs: density-gives density in units of Ne/sqrt(Te/10^4), judgement-string that says it's low, medium, or high
       This function uses the figure in slide deck 15b from Robin's lecture series given the two Oxygen line fluxes to get whether this density is low, medium, or high, if it's low, we can use temperature diagnostics'''
    
    density_data = np.loadtxt("Figures/density_diagnostic_3729_3726_ratio.txt")
    dens = density_data[:,0]
    intensity = data[:,1]
    idx = find_nearest(wv = intensity, value = f_3729/f_3726)
    
    density = dens[idx]
    
    judgement = ""
    if dens < 100:
        judgement = "low"
    elif dens < 3.0E3:
        judgement = "med"
    else:
        judgement = "high"
    
    return {"density": density, "judgement": judgement}
    
#temperature diagnostics
def temp_diagnostic_low_density(f_5007, f_4959, f_4363):
    '''inputs: f_5007-flux in 5007 angstroms, f_4959-flux in 4959 angstroms, f_4363-flux in 4363 angstroms
       outputs: temp-temperature found from plot
       This function uses the figure in slide deck 15b from Robin's lecture series to get temperature given these three Oxygen line fluxes. Remember that this only works if the Ne density is low-check that first!'''
       
    temp_data = np.loadtxt("Figures/temp_diagnostic_5007_4959_4363.txt")
    temps = temp_data[:,0]
    ratios = temp_data[:,1]
    idx = find_nearest(wv=ratios, value = (f_5007+f_4959)/f_4363)
    
    temp = temps[idx]
    
    return temp


def get_temp_density_diagnostic(f_5007, f_4959, f_4363, f_3729, f_3726):
    '''inputs: f_3729-flux in 3729 angstroms, f_3726-flux in 3726 angstroms, f_5007-flux in 5007 angstroms, f_4959-flux in 4959 angstroms, f_4363-flux in 4363 angstroms
       outputs: temp-temperature found, density-Ne in cm^-3
       This function uses the figures in slide deck 15b from Robin's lecture series to get density using the two oxygen forbidden lines 3729 and 3726, and if density is low enough, then uses the Oxygen forbidden lines 5007, 4959, and 4363 to get temperature'''
    dd = density_diagnostic(f_3729 = f_3729, f_3726 = f_3726)
    density_odd_units = dd['density']
    judgement = dd['judgement']
    
    if judgement == 'low':
        temp_K = temp_diagnostic_low_density(f_5007 = f_5007, f_4959 = f_4959, f_4363 = f_4363)
        density_cm3 = density_odd_units*np.sqrt(temp_K/10000)
        return {"temp_K": temp_K, "density_cm3": density_cm3}
    
    return {"density_odd_units": density_odd_units}

#cosmological distance measurements

def physical_dist_from_comoving_dist(comoving_dist, scale_factor):
    '''inputs: comoving_dist-distance to object in comoving coordinates, scale_factor-scale factor of universe
       outputs: physical_dist-physical distance to object (units are same as units in comoving_dist
       this function takes in a comoving distance and gives a physical distance to an object. The units of the answer are the same as would come from comoving_dist'''
       
    
    return comoving_dist*scale_factor
    
def z_from_vel(vel_ms):
    '''inputs: vel_ms-velocity of object in meters per second
       outputs: z-redshift corresponding to this velocity
       This function gives back redshift as a function of velocity of an object. We have the formula 1+z = (1+v/c)*gamma. I just want to point out: gamma is independent of the direction of the velocity. So, gamma is always 1/sqrt(1-beta^2), but the classical aspect of the doppler shift has a cos \theta in it because we need to see what velocity the object is moving away from us. So, the formula is 1+z = (1+v*cos(theta)/c)/(sqrt(1-beta^2)). Theta is zero when the object is moving directly away from us, which is the assumption this function makes. More on this at https://en.wikipedia.org/wiki/Redshift#Extragalactic_observations'''
       
    z = np.sqrt((1+vel_ms/3.0E8)/(1-vel_ms/3.0E8))-1.0
    
    return z
    
       
       
       
    
    
def z_from_scale_factor(a, a_0=1.0):
    '''inputs: a-scale factor of universe at some time t, a_0-(optional) scale factor of universe at time 0-defaults to t = now and a_0=1.0
       outputs: z-redshift measurement of the object at time t, when seen from time 0
       This function gets the redhisft measurement of an object given the scale factor of the universe at that time. This function is flexible enough to get the redshift that one object in the past would see another object in the further past. We use the equation (1+z) = a_0/a'''
    
    return 1-(a_0/a)
    
def scale_factor_from_z(z, a_0 = 1.0):
    '''inputs: z-redshift of an object in the distance, a_0-(optional) scale factor of the universe when the object was observed-defaults to a_0=1.0 corresponding to now
       outputs: a-scale factor of the universe when this object emitted the light being observed
       This function turns a redshift into a scale factor using the equation (1+z) = a_0/a'''

    return a_0/(1+z)
    
def length_conversions(d_cm = None, d_m = None, d_Km = None, d_AU = None, d_ly = None, d_pc = None, d_Kpc = None, d_Mpc = None, d_Gpc = None):
    '''inputs: d_SI-distance in SI units, d_MPC-distance in megaparsec
       outputs: dist_dic-dictionary of distance conversions of d_SI into a bunch of distance units
       This function is a one-stop place to get all distance conversions. We go from meters into cm, m, km, AU, ly, pc, kpc, MPC, GPC. We can also take in whatever unit of distance and give it back in all of these units'''
    
    conv_dic = {"cm": 1.0E-2, "m":1.0, "Km": 1000., "AU": 1.496E11, "ly": 9.461E15, "pc": 3.086E16, "Kpc": 3.086E19, "Mpc": 3.086E22, "Gpc": 3.086E25}
    
    inputs = {"cm": d_cm, "m":d_m, "Km": d_Km, "AU": d_AU, "ly": d_ly, "pc": d_pc, "Kpc": d_Kpc, "Mpc": d_Mpc, "Gpc": d_Gpc}
    
    outputs = {"cm": None, "m":None, "Km": None, "AU": None, "ly": None, "pc": None, "Kpc": None, "Mpc": None, "Gpc": None}
    
    input_unit = None
    input_val = None
    for unit in inputs.keys():
        if inputs[unit] is not None:
            input_unit = unit
            input_val = inputs[unit]
            
    #get into meters first, then into <unit> units
    in_meters = input_val*conv_dic[input_unit]
    for unit in outputs.keys():
        outputs[unit] = in_meters/conv_dic[unit]
        
    dist_mod = 5*np.log10(outputs['pc']/10)
    outputs["dist_mod"] = dist_mod
    return outputs
    

    
def dist_from_redshift(z, hubble_const = H_0_SI):
    '''inputs: z-redshift of an object, hubble_const-hubble constant to be used (defaults to value set above)
       outputs: dist-distance to the object, given its redshift.
       This uses the formula D_p = \frac{2c}{H0}[1-(1+z)^-0.5] to get distance, given redshift to an object. Note that this formula only works for small redshifts'''
       
    if z > 1:
        raise Exception("Redshift passed in: "+str(z)+" is too large to use approximate formula for physical distance.")
    elif z > 0.1:
        warnings.warn("Warning: Redshift passed into distance formula might be too large.")
        
    dist_m = (2*c_ms/hubble_const)*(1-(1+z)**-0.5)
    d = length_conversions(d_m = dist_m)
    
    return d
    
    
def z_from_dist(dist_m, hubble_const = H_0_SI):
    '''inputs: dist_m-distance to the object, given its redshift, hubble_const-hubble constant to be used (defaults to value set above)
       outputs: z-redshift of an object
       This uses the formula D_p = \frac{2c}{H0}[1-(1+z)^-0.5] to get redshift, given distance to an object. Note that this formula only works for small redshifts'''
       
       
    
    z = (1-dist_m*hubble_const/(2*c_ms))**-2.0-1.0
    
    if z > 1:
        raise Exception("Redshift passed in: "+str(z)+" is too large to use approximate formula for physical distance.")
    elif z > 0.1:
        warnings.warn("Warning: Redshift passed into distance formula might be too large.")
        
    
    return z
    
#lookback time as a function of z, as a function of a

def conv_time_sec(time_sec):
    '''inputs: time_sec-time measured in seconds
       outputs: time_dic-dictionary of time_sec converted into various units, including minutes, hours, solar days, years, Myr, and Gyr
       This function converts time in seconds to a bunch of other units, returning all of them in a dictionary'''
    
    time_min = time_sec/60.0
    time_hr = time_sec/3600.0
    time_day = time_sec/86400.0
    time_yr = time_day/365.24219
    time_Myr = time_yr/1.0E6
    time_Gyr = time_yr/1.0E9
    
    return {"time_sec":time_sec,"time_min":time_min, "time_hr":time_hr, "time_day":time_day, "time_yr": time_yr, "time_Myr":time_Myr, "time_Gyr":time_Gyr}



def milne_lookback_time(z, hubble_const = H_0_SI):
    '''inputs: z-redshift of some distant object, hubble_const-(optional) current hubble constant of the universe-defaults to the value defined above
       outputs: delta_t-dictionary of lookback time measurment from now (t_0) to some object at redshift z in various units
       This function gets the lookback time to an object at redshift z using the equation t = (1/H_0)(1-1/(1+z)). Note that this assumes the Milne cosmology, which says that \rho = 0 and that the universe is empty of matter'''
    
    lookback_time_sec = (1.0/hubble_const)*(1.0-1.0/(1.0+z))
    
    time_dic = conv_time_sec(time_sec = lookback_time_sec)
    lookback_time_yr = time_dic['time_yr']
    lookback_time_Gyr = time_dic['time_Gyr']
    return {"lookback_time_sec":lookback_time_sec, "lookback_time_yr":lookback_time_yr, "lookback_time_Gyr":lookback_time_Gyr}



def de_sitter_lookback_time(z, hubble_const = H_0_SI):
    '''inputs: z-redshift of some distant object, hubble_const-(optional) current hubble constant of the universe-defaults to the value defined above
       outputs: delta_t-dictionary of lookback time measurment from now (t_0) to some object at redshift z in various units
       This function gets the lookback time to an object at redshift z using the equation t = (2/3)(1/H_0)(1-(1+z)^(-3/2)). Note that this assumes the Einstein-de Sitter cosmology, which says that E = 0 and the universe is a critical universe, where total energy is zero'''
       
    lookback_time_sec = (2.0/3.0)*(1.0/hubble_const)*(1.0-(1.0+z)**(-1.5))
    
    time_dic = conv_time_sec(time_sec = lookback_time_sec)
    lookback_time_yr = time_dic['time_yr']
    lookback_time_Gyr = time_dic['time_Gyr']
    return {"lookback_time_sec":lookback_time_sec, "lookback_time_yr":lookback_time_yr, "lookback_time_Gyr":lookback_time_Gyr}

def radians_from_HHMMSS(angle_HHMMSS):
    '''inputs: angle_HHMMSS-angle in units of HH:MM:SS
       outputs: angle_rad-the angle passed in units of radians
       This function turns an angle in units of hours, minutes, seconds into an angle in radians'''
       
    g = angle_HHMMSS.split(":")
    HH = float(g[0])
    
    
    MM = float(g[1])
    SS = float(g[2])
    
    if HH < 0:
        MM *= -1.0
        SS *= -1.0
    
    #print("HH: "+str(HH))
    #print("MM: "+str(MM))
    #print("SS: "+str(SS))
    
    angle_frac = HH/24. + MM/(24.*60.) + SS/(24.*60.*60.)
    angle_rad = angle_frac * 2.*np.pi
    
    return angle_rad
       
    


def HHMMSS_from_radians(angle_rad):
    '''inputs: angle_rad-the angle passed in units of radians
       outputs: angle_HHMMSS-angle in units of HH:MM:SS
       This function turns an angle in radians into an angle in units of hours, minutes, seconds'''


    if angle_rad < 0:
        angle_rad += 2*np.pi
        
        
    angle_frac = angle_rad/(2.*np.pi)
    
    hours = angle_frac*24.
    frac_hours = hours - math.floor(hours)
    
    mins = frac_hours * 60.
    frac_mins = mins - math.floor(mins)
    
    secs = frac_mins * 60.
    
    return str(math.floor(hours))+":"+str(math.floor(mins))+":"+str(round(secs, 5))


def equatorial_from_horizontal(altitude, azimuth, latitude_deg, time_of_day, days_after_vernal_equinox, longitude_deg, UTC_offset):
    '''inputs: altitude-altitude angle of some object in radians, azimuth-azimuth angle of some object in radians, latitude-latitude of location on Earth in radians
       outputs: RA-right ascension of the object in HH:MM:SS, dec-declination of the object in radians
       This function takes altitude, azimuth, and latitude to give the right ascension and declination of some object measured on the sky plane.'''
       
    
    latitude = latitude_deg * np.pi/180.
    #print("altitude: "+str(altitude))
    #print("azimuth: "+str(azimuth))
    #print("lat: "+str(latitude))
    #print("time of day: "+str(time_of_day))
    #print("days_after..: "+str(days_after_vernal_equinox))
    
    current_LST_hrs = get_LST(local_time_of_day = time_of_day, days_after_vernal_equinox = days_after_vernal_equinox, longitude_deg=longitude_deg, UTC_offset = UTC_offset)['LST_hrs']
    
    
    #print("current_LST_hrs: "+str(current_LST_hrs))
    
    current_LST_radians = current_LST_hrs/24.*2*np.pi
    
    #print("current_LST_radians: "+str(current_LST_radians))
    
    #we'll use the equation sin(dec) = -cos(azimuth)*cos(altitude)*cos(lat) + sin(altitude)*sin(lat) to get declination first
    
    dec = math.asin(math.cos(azimuth)*math.cos(altitude)*math.cos(latitude)+math.sin(altitude)*math.sin(latitude)) #declination can only be between -np.pi/2 to np.pi/2, so there's no degeneracy in this function
    
    #print("dec: "+str(dec))
    #print("")
    

    #I need to use another equation to tell which of h1 or h2 is the correct right ascension
    #I can use sin(h)*cos(dec) = -sin(azimuth)*cos(altitude) to get h3, the asin(..) solution
    
    h3 = math.asin(-1.0*math.sin(azimuth)*math.cos(altitude)/math.cos(dec))

    
    h_radians = h3
    #print("found h_radians: "+str(h_radians)+" but I wanted to use: "+str(A[combo[0]]))
    
    RA_radians = current_LST_radians-h_radians
    RA_HHMMSS = HHMMSS_from_radians(RA_radians)
    
    dec_rad = dec
    dec_HHMMSS = HHMMSS_from_radians(dec_rad)
    dec_deg = dec_rad*180./np.pi
    dec_deg_m_s = deg_min_sec_from_degree(dec_deg)
    
    RA = RA_HHMMSS
    dec = dec_deg_m_s
    
    return {"RA": RA, "dec": dec}
    
    
    
    
def get_UTC_offset(longitude_deg):
    '''inputs: longitude_deg-longitude of location in degrees (negative for west, positive for east)
       outputs: UTC_offset-an integer telling how many hours to add to local time at this longitude to get UTC time
       This function gives how many hours behind the longitude is from UTC time'''
    from datetime import datetime
    datetime.utcnow()
    
    hours_offset = math.round(longitude/15)
    
    return hours_offset #negative 1 because we want how many hours <<to add>>
    

def get_LST(local_time_of_day, days_after_vernal_equinox, longitude_deg, UTC_offset = -45):
    '''inputs: time_of_day-the time of day in format HH:MM:SS as a string, days_after_vernal_equinox-days since vernal equinox for the date in question, longitude_deg-the longitude of the location on Earth measured in degrees (negative for West and positive for East), UTC_offset-(optional) the number of hours to add to UTC to get local time this defaults to -45 and if it is set to -45 this will attempt to calculate UTC time using the longitude
       outputs: LST-corresponding local sidereal time of this time measure
       This function gets local sidereal time for a given time measurement on a given day. The equation I am using is LST = T + 12hours + (days after vernal equinox)*4 min'''
    
    GG_TOD = local_time_of_day.split(":")
    local_time = float(GG_TOD[0])+float(GG_TOD[1])/60.+float(GG_TOD[2])/3600.
    
    #print("local time: "+str(local_time))
    #1) go from local time to UTC
    
    if UTC_offset == -45:
        UTC_offset = get_UTC_offset(longitude_deg = longitude_deg)
    
    UTC_time = local_time-UTC_offset
    
    #print("UTC time: "+str(UTC_time))
    
    
    #2) get GST-sidereal time at the Greenwich meridian
    
    GST_hrs = (UTC_time + 12 + days_after_vernal_equinox*(3./60.+55./3600.))%24
    #print("GST: "+str(GST_hrs))
    
    
    #3) convert from GST to LST using the local longitude in HHMMSS format
    
    longitude_hrs = longitude_deg/360. * 24.
    LST_hrs = GST_hrs+longitude_hrs
    
    
    
    #LST_hrs = time_of_day_hrs + 12 + days_after_vernal_equinox*4./60.
    #LST_hrs = (time_of_day_hrs + 12 + days_after_vernal_equinox*3.9333/60.)%24.
    
    LST_hr_floor = math.floor(LST_hrs)
    
    LST_mins = (LST_hrs - LST_hr_floor)*60.
    LST_mins_floor = math.floor(LST_mins)
    
    LST_secs = round((LST_mins - LST_mins_floor)*60.,4)
    
    LST_HHMMSS = str(LST_hr_floor)+":"+str(LST_mins_floor)+":"+str(LST_secs)
    
    return {"LST_HHMMSS": LST_HHMMSS, "LST_hrs": LST_hrs}
    
    
def get_LST_good(local_time, JD, longitude_deg, UTC_offset):
    
    GG_TOD = local_time.split(":")
    local_time = float(GG_TOD[0])+float(GG_TOD[1])/60.+float(GG_TOD[2])/3600.
    
    #print("local time: "+str(local_time))
    #1) go from local time to UTC
    
    if UTC_offset == -45:
        UTC_offset = get_UTC_offset(longitude_deg = longitude_deg)
    
    UTC_time = local_time-UTC_offset
    
    T = (JD - 2451545)/36525
    theta = 280.46061837 + 360.98564736629 * (JD - 2451545.0) + (0.000387933 * T * T) - (T * T * T / 38710000.0)
    
    A = theta*360/(2*np.pi) + longitude_deg
    print("A: "+str(A*24/360 % 24))
    
    return theta%24
    
    

def horizontal_from_equatorial(declination, RA, latitude_deg, time_of_day, days_after_vernal_equinox, longitude_deg, UTC_offset):
    #a is altitude, A is azimuth
    
    latitude = latitude_deg * np.pi/180.
    
    current_LST_hrs = get_LST(local_time_of_day = time_of_day, days_after_vernal_equinox = days_after_vernal_equinox, longitude_deg=longitude_deg, UTC_offset = UTC_offset)['LST_hrs']
    
    current_LST_radians = (current_LST_hrs/24.)*2.*np.pi
    #print("declination: "+str(declination))
    #print("RA: "+str(RA))
    #print("current_LST_hrs: "+str(current_LST_hrs))
    
    h = current_LST_radians-RA #this works a good bit better if I subtract 0.5 from h
    
    
    
    
    #print("current_LST_radians:" +str(current_LST_radians))
    #print("RA: " +str(RA))
    #print("h: "+str(h))
    
    
    a = math.asin(math.cos(h)*math.cos(declination)*math.cos(latitude)+math.sin(declination)*math.sin(latitude))
    
    #use cos(a)*sin(A) = -cos(delta)*sin(H)
    A = math.asin(-1.0*math.cos(declination)*math.sin(h)/math.cos(a))
    
    if A < 0:
        A += 2.*np.pi
    
    return{"alt": deg_min_sec_from_degree(a*180./np.pi), "az":deg_min_sec_from_degree(A*180./np.pi)}


def equatorial_from_ecliptic(beta_deg, Lambda_deg):
    '''inputs: beta-beta coordinate in degrees of an object in ecliptic coordinates, Lambda_deg-lambda coordinate in degrees of an object in ecliptic coordinates
       outputs: RA_HMS-Right Ascension in units of hours, minutes, seconds, dec_deg_m_s-declination of object in degrees, minutes, and seconds
       This function takes the ecliptic coordinates of an object and gives its equatorial coordinates'''
       
    obliquity_deg = 23.5 #obliquity of the Earth's orbit wrt ecliptic plane
    #sin delta = cos (beta)*sin(lambda)*sin(obliquity)+sin(beta)*cos(obliquity)
    eps = obliquity_deg*np.pi/180
    beta = beta_deg *np.pi/180
    lam = Lambda_deg*np.pi/180
    
    
    dec = math.asin(math.cos(beta)*math.sin(lam)*math.sin(eps)+math.sin(beta)*math.cos(eps))
    
    alpha = math.acos(math.cos(beta)*math.sin(lam)/math.cos(dec))
    
    
    dec_deg_m_s = deg_min_sec_from_degree(dec)
    RA_HMS = HHMMSS_from_radians(alpha)
    
    return {"dec_deg_m_s":dec_deg_m_s, "RA_HMS":RA_HMS}
    
    
def ecliptic_from_equatorial(dec_deg_m_s, RA_HMS):
    '''inputs: RA_HMS-Right Ascension in units of hours, minutes, seconds, dec_deg_m_s-declination of object in degrees, minutes, and seconds
       outputs: beta-beta coordinate in degrees of an object in ecliptic coordinates, Lambda_deg-lambda coordinate in degrees of an object in ecliptic coordinates
       This function takes the equatorial coordinates of an object and gives its ecliptic coordinates'''
       
    obliquity_deg = 23.5 #obliquity of the Earth's orbit wrt ecliptic plane
    
    eps = obliquity_deg*np.pi/180
    
    dec = rad_from_degree_min_sec(dec_deg_m_s)
    alpha = radians_from_HHMMSS(RA_HMS)
    
    
    #sin beta = sin(dec)*cos(obliquity)-cos(delta)*sin(alpha)*sin(eps)
    beta = math.asin(math.sin(dec)*math.cos(eps)-math.cos(dec)*math.sin(alpha)*math.sin(eps))
    
    #cos(delta)*cos(alpha) = cos(beta)*sin(lambda)
    lam = math.asin(math.cos(dec)*math.cos(alpha)/math.cos(beta))
    
    
    
    beta_deg_m_s = deg_min_sec_from_degree(beta)
    
    
    
    lam_HMS = HHMMSS_from_radians(lam)
    
    return {"beta_deg_m_s":beta_deg_m_s, "lam_HMS":lam_HMS}
    
    
def equatorial_from_galactic(b_deg, l_deg, l_0_deg = 33, dec_0_deg = 62.6, alpha_0_deg = 282.25):
    '''inputs: b_deg-b coordinate in degrees, l_deg-l coordinate in degrees, l_0_deg-l_0 constant in degrees (default to 33), dec_0_deg-dec_0 constant in degrees (default to 1950 value of 62.6), alpha_0_deg-alpha_0 constant in degrees (default to 1950 value of 282.25 degrees)
       outputs: RA_HMS-right ascension value in units if hours minutes seconds, dec_deg_m_s-declination of object in degrees minutes and seconds
       This function goes from galactic coordinates to equatorial RA-dec coordinates. It defaults to 1950 values for the RA and dec constants and for the longitude of the center of the galaxy'''
       
    b = b_deg*np.pi/180.
    l = l_deg*np.pi/180.
    l_0 = l_0_deg*np.pi/180.
    dec_0 = dec_0_deg*np.pi/180.
    alpha_0 = alpha_0_deg*np.pi/180.
    
    #sin dec = cos b * sin(l-l_0)sin dec_0 + sin b * cos dec_0
    
    dec = math.asin(math.cos(b)*math.sin(l-l_0)*math.sin(dec_0)+math.sin(b)*math.cos(dec_0))
    
    #cos b cos (l-l_0) = cos dec * cos (alpha-alpha_0)
    
    alpha = (math.acos(math.cos(b)*math.cos(l-l_0)/math.cos(dec))+alpha_0)%(2*np.pi)
    
    dec_deg_m_s = deg_min_sec_from_degree(dec*180./np.pi)
    
    print("dec_deg: "+str(dec*180./np.pi))
    
    RA_HMS = HHMMSS_from_radians(alpha)
    
    return {"dec_deg_m_s":dec_deg_m_s, "RA_HMS":RA_HMS}
    


def galactic_from_equatorial(dec_deg_m_s, RA_HMS, l_0_deg = 33, dec_0_deg = 62.6, alpha_0_deg = 282.25):
    '''inputs: RA_HMS-right ascension value in units if hours minutes seconds, dec_deg_m_s-declination of object in degrees minutes and seconds, l_0_deg-l_0 constant in degrees (default to 33), dec_0_deg-dec_0 constant in degrees (default to 1950 value of 62.6), alpha_0_deg-alpha_0 constant in degrees (default to 1950 value of 282.25 degrees)
       outputs: b_deg-b coordinate in degrees, l_deg-l coordinate in degrees
       This function goes from galactic coordinates to equatorial RA-dec coordinates. It defaults to 1950 values for the RA and dec constants and for the longitude of the center of the galaxy'''
       
    dec = rad_from_degree_min_sec(dec_deg_m_s)
    alpha = radians_from_HHMMSS(RA_HMS)
    l_0 = l_0_deg*np.pi/180.
    dec_0 = dec_0_deg*np.pi/180.
    alpha_0 = alpha_0_deg*np.pi/180.
    
    #sin b = sin dec * cos(dec_0)-cos dec * sin(alpha-alpha_0)*sin dec_0
    
    b = math.asin(math.sin(dec)*math.cos(dec_0)-math.cos(dec)*math.sin(alpha-alpha_0)*math.sin(dec_0))
    
    #cos b * cos(l-l_0) = cos dec * cos(alpha-alpha_0)
    
    l = math.acos(math.cos(dec)*math.cos(alpha-alpha_0)/math.cos(b)) + l_0
    
    
    
    b_deg_m_s = deg_min_sec_from_degree(b*180./np.pi)
    l_HMS = HHMMSS_from_radians(l)
    
    return {"b_deg_m_s":b_deg_m_s, "l_HMS":l_HMS}
    
    
    
def precess_ecliptic_coord_fwd(num_years_to_precess, lam_old_deg, beta_old_deg, del_lam_arcsec_yr = 50.29, del_beta_arcsec_yr = 0.0):
    '''inputs: num_years_to_precess-number of years that have passed since lam_old and beta_old measurements, lam_old_deg-old ecliptic longitude that needs to be precessed forward in degrees, beta_old-old ecliptic latitude that needs to be processed forward in degrees, del_lam_arcsec_yr-precession rate of lambda in arcsec/yr, del_beta_arcsec_yr-precession rate of beta in arcsec/yr
       outputs: beta_new-new ecliptic latitude after precessing beta_old forward in degrees, lam_new-new ecliptic longitude after precessing lam_old forward in degrees
       This function precesses ecliptic coordinates forward in time, as a result of the Earth's spin axis precession'''
    
    lam_old_rad = lam_old_deg*np.pi/180.
    beta_old_rad = beta_old_deg*np.pi/180.
    
    lam_old_arcsec = lam_old_rad*206265.
    beta_old_arcsec = beta_old_rad*206265.
    
    lam_new_arcsec = lam_old_arcsec+del_lam_arcsec_yr*num_years_to_precess
    beta_new_arcsec = beta_old_arcsec+del_beta_arcsec_yr*num_years_to_precess
    
    lam_new_deg = lam_new_arcsec*180./(np.pi*206265)
    beta_new_deg = beta_new_arcsec*180./(np.pi*206265)
    
    print('lam_new_deg: '+str(lam_new_deg))
    print('beta_new_deg: '+str(beta_new_deg))
    
    beta_new_deg_m_s = deg_min_sec_from_degree(beta_new_deg)
    lam_HMS = HHMMSS_from_radians(lam_new_deg*np.pi/180.)
    
    return {"lam_new_deg":lam_new_deg, "beta_new_deg":beta_new_deg}
    
    
def precess_ecliptic_coord_back(num_years_to_precess, lam_new_deg, beta_new_deg, del_lam_arcsec_yr = 50.29, del_beta_arcsec_yr = 0.0):
    '''inputs: num_years_to_precess-number of years to precess back to, lam_new_deg-new ecliptic longitude that needs to be precessed backwards in degrees, beta_new-new ecliptic latitude that needs to be processed forward in degrees, del_lam_arcsec_yr-precession rate of lambda in arcsec/yr, del_beta_arcsec_yr-precession rate of beta in arcsec/yr
       outputs: beta_old-old ecliptic latitude after precessing beta_new backwards in degrees, lam_od-old ecliptic longitude after precessing lam_new backewards in degrees
       This function precesses ecliptic coordinates backwards in time, as a result of the Earth's spin axis precession'''
    
    lam_new_rad = lam_new_deg*np.pi/180.
    beta_new_rad = beta_new_deg*np.pi/180.
    
    lam_new_arcsec = lam_new_rad*206265.
    beta_new_arcsec = beta_new_rad*206265.
    
    lam_old_arcsec = lam_new_arcsec-del_lam_arcsec_yr*num_years_to_precess
    beta_old_arcsec = beta_new_arcsec-del_beta_arcsec_yr*num_years_to_precess
    
    lam_old_deg = lam_old_arcsec*180./(np.pi*206265)
    beta_old_deg = beta_old_arcsec*180./(np.pi*206265)
    
    beta_old_deg_m_s = deg_min_sec_from_degree(beta_old_deg)
    lam_HMS = HHMMSS_from_radians(lam_old_deg*np.pi/180.)
    
    return {"lam_old_deg":lam_old_deg, "beta_old_deg":beta_old_deg}



def precess_equatorial_coordinates_fwd(num_years_to_precess, alpha_old_deg, del_old_deg):
    '''inputs: num_years_to_precess-number of years that have passed since alpha_old_deg and del_old_deg measurements, alpha_old_deg-old right ascension that needs to be precessed forward in degrees, del_old_deg-old declination that needs to be processed forward in degrees
    outputs: del_new_deg-new declination after precessing del_old_deg forward in degrees, alpha_new_deg-new right ascension after precessing alpha_old_deg forward in degrees
    This function precesses equatorial coordinates forward in time, as a result of the Earth's spin axis precession'''
    alpha_old_rad = alpha_old_deg*np.pi/180.
    del_old_rad = del_old_deg*np.pi/180.



    del_alpha_sec = (3.07234+1.3365*math.sin(alpha_old_rad)*math.tan(del_old_rad))

    del_alpha_rad = 2.*np.pi*del_alpha_sec/(3600.*24.)
    alpha_new_rad = alpha_old_rad +del_alpha_rad*num_years_to_precess

    del_new_rad = del_old_rad + ((20.0468*math.cos(alpha_old_rad))*num_years_to_precess)/206265.

    alpha_new_deg = alpha_new_rad * 180./np.pi
    del_new_deg = del_new_rad*180./np.pi

    return {"alpha_new_deg":alpha_new_deg, "del_new_deg":del_new_deg}
       
       
    
    
    
    
    
    
    

    
    
    
    

    
    
