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
                'FLUX': 'Y',
                'PRI':'0',
                'SETUPMETHOD':'DirectGuider',
                'DITHER':'Y',
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
            ra = c[i].split(' ')[0].replace('h',':').replace('m',':').replace('s','')
            dec = c[i].split(' ')[1].replace('d',':').replace('m',':').replace('s','')
            f.write('GW{}\t{}\t{}\t{}\n'.format(target[1]+target[2],ra,dec,str(int(target[0]))))
def main():
    make_phaseii(sys.argv[1])
if __name__=='__main__':
    main()
