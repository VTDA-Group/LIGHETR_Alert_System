import gcn
from tabulate import tabulate as tb
import healpy as hp
from prob_obs import prob_observable
from astropy.time import Time
import astropy
from astropy.utils.data import download_file
import email_ip, get_galaxies, get_LST, make_phaseii
import os, urllib
import argparse
import urllib.request
import lxml
import numpy as np
import json
import matplotlib.pyplot as plt
import lxml.etree
import sys
import datetime


def process_fits(recipients, params = None):
        skymap,header = hp.read_map(params['skymap_fits'], h=True, verbose=False)

        # Print and save some values from the FITS header.
        header = dict(header)
        params['time'] = Time(header['DATE-OBS'],format='isot',scale='utc')
        time = Time.now()
        params['Distance'] = str(header['DISTMEAN']) + ' +/- ' + str(header['DISTSTD'])
        header['GraceID'] = params['GraceID']
        with open('./'+params['GraceID']+'.dat','w') as f:
            data_p = []
            tableheaders = ['PARAMETER','VALUE']
            for prm in interesting_parameters:
                if prm in list(params.keys()):
                    data_p.append((prm,params[prm]))
            print(tb(data_p,headers=tableheaders,tablefmt='fancy_grid'))
            f.write(tb(data_p,headers=tableheaders,tablefmt='html'))

        # Making a pie chart of the type of event for the email
        labels = ['BNS', 'NSBH', 'BBH', 'MassGap', 'Terrestrial']
        sizes = [float(params[label])*100 for label in labels]
        labels = ['%s (%.1f %%)'%(lab,pct) for lab,pct in zip(labels,sizes)]
        fig1, ax1 = plt.subplots()
        patches,texts = ax1.pie(sizes, startangle=90)
        ax1.legend(patches, labels, loc="best")
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.savefig('piechart_%s.png'%params['GraceID'])

        prob, probfull, timetill90, m = prob_observable(skymap, header, time, plot=plotting)

        params['skymap_array'] = m

        if timetill90 ==-99:
            print("HET can't observe the source.")
            return
        else:
            print("Source has a {:.1f}% chance of being observable now.".format(int(round(100 * prob))))
            print("Integrated probability over 24 hours (ignoring the sun) is {:.1f}%".format(int(round(100 * probfull))))
            print('{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90))
            send_notifications(params,timetill90,text=True,email=False)
            get_galaxies.write_catalog(params,'MANGROVE')
            mincontour = get_LST.get_LST(targf = 'galaxies%s_%s.dat'%('MANGROVE',params['GraceID']))
            make_phaseii.make_phaseii('LSTs_{}.out'.format(params['GraceID']))
            send_notifications(params,timetill90)

