
class Observatory(astropy.coordinates.EarthLocation):
    def __init__(self, name, dec_range, loc):
        self.name = name
        self.min_dec, self.max_dec = dec_range # in degrees        
        self.min_dec_rad = (90-self.min_dec)*np.pi/180
        self.max_dec_rad = (90-self.max_dec)*np.pi/180

        self.pupil = np.loadtxt('hetpix.dat')

        # set EarthLocation object
        super().__init__(
            lat=loc[1]*u.deg,
            lon=loc[0]*u.deg,
            height=loc[2]*u.m
        )
    
        
    