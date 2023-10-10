import numpy as np
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord



def make_phaseii(lstfile, savedir = ''):
    common = {
                'PROGRAM':'HET23-2-400',
                'VIFU':'047',
                'EXP':'1080',
                'CRSPLIT':'3',
                'INSTRUMENT':'VIRUS',
                'GMAG':'22',
                'SKYBRIGHT_G':'18.0',
                'SEEING': '3.0',
                'SKYTRANS': 'S',
                'SKYCALS': 'Y',
                'PRI':'0',
                'SETUPMETHOD':'DirectGuider',
                'DITHER':'Y',
                'PMRA':'0',
                'PMDEC':'0',
                'COMMENT':'"Usual Dither, look for new object in target IFU"',
                }
    GraceID = lstfile.split('_')[1].split('.')[0]
    with open(savedir+'submission_to_HET.tsl','w') as f:
        f.write('COMMON\n')
        for key,value in common.items():
            f.write('\t{}\t{}\n'.format(key,value))
        f.write('TRACK_LIST\n')
        f.write(' OBJECT\tRA\tDEC\tPIPRI\n')
        targets = np.loadtxt(lstfile,skiprows=1,dtype=str)
        targets = np.atleast_2d(targets)
        c = SkyCoord(ra=np.asarray(targets[:,1], dtype=float) \
            * u.degree, dec = np.asarray(targets[:,2], \
                dtype=float)* u.degree, frame='icrs')
        c = c.to_string('hmsdms')
        for i,target in enumerate(targets):
            
            #process ra into format:
            ra = c[i].split(' ')[0]
            hour = "{:2.0f}".format(float(ra[:ra.index('h')]))
            if hour[0] == ' ':
                hour = '0'+hour[1:]
            min = "{:2.0f}".format(float(ra[ra.index('h')+1:ra.index('m')]))
            if min[0] == ' ':
                min = '0'+min[1:]
            sec = "{:2.2f}".format(float(ra[ra.index('m')+1:ra.index('s')]))
            if sec[0] == ' ':
                sec = '0'+sec[1:]
            if len(sec) == 4:
                sec = '0'+sec
            ra = hour+":"+min+":"+sec
            
            #processing dec into format
            dec = c[i].split(' ')[1]
            
            pos_neg = dec[0]
            deg = "{:2.0f}".format(float(dec[1:dec.index('d')]))
            if deg[0] == ' ':
                deg='0'+deg[1:]
            min = "{:2.0f}".format(float(dec[dec.index('d')+1:dec.index('m')]))
            if min[0] == ' ':
                min='0'+min[1:]
            sec = "{:2.2f}".format(float(dec[dec.index('m')+1:dec.index('s')]))
            if sec[0] == ' ':
                sec = '0'+sec[1:]
            if len(sec) == 4:
                sec = '0'+sec
            dec = pos_neg+deg+":"+min+":"+sec
            f.write('Target{}\t{}\t{}\t{}\n'.format(target[0],ra,dec,str(int(target[0]))))
def main():
    make_phaseii(sys.argv[1])
if __name__=='__main__':
    main()
