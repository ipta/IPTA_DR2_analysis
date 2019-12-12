#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 17:12:59 2019

@author: dreardon

Script for generating cut datasets for those pulsars with known band/system noise, mainly:

    J1713 (legacy)
    J0437 (legacy)
    J1939 (legacy)
    J2145 (legacy)
    J1939, J0437, J1643, J1600, low frequency cut

"""

from __future__ import division, print_function, unicode_literals

import os, glob
import libstempo as t2

import dr2lite_utils as dr2u

DR2DATA = os.path.abspath('/Users/dreardon/Dropbox/git/DR2/')  # path to data local usage
datadir = os.path.join(DR2DATA, 'release/VersionB')

outdir = 'data/partim_classic'
os.system('mkdir -p {}'.format(outdir));

parfiles = glob.glob(datadir + '/J*/*IPTADR2.par')

psr_names = []
for p in parfiles:
    name = p.split('/')[-2]
    psr_names.append(name)
    outfile = os.path.join(outdir, '{}.par'.format(name))
    dr2u.clean_par(p, outfile)

timfiles = glob.glob(datadir + '/J*/*IPTADR2.tim')

for t in timfiles:
    name = t.split('/')[-2]
    outfile = os.path.join(outdir, '{}.tim'.format(name))
    dr2u.combine_tim(t, outfile)

psr_names.sort()

psrfile = os.path.join(outdir, 'psrlist_classic.txt')
with open(psrfile, 'w') as f:
    for pname in psr_names:
        f.write("{:s}\n".format(pname))

filt = {'group':
                ['327_ASP', '430_ASP', 'L-wide_ASP', 'S-wide_ASP',
                 '327_PUPPI', '430_PUPPI', 'L-wide_PUPPI',  'S-wide_PUPPI',
                 'Rcvr_800_GASP', 'Rcvr1_2_GASP',
                 'Rcvr_800_GUPPI', 'Rcvr1_2_GUPPI',
                 'PDFB_10CM', 'PDFB_20CM', 'PDFB_40CM',
                 'CPSR2_20CM', 'CPSR2_50CM',
                 'WBCORR_10CM', 'WBCORR_20CM',
                 'EFF.EBPP.1360', 'EFF.EBPP.1410', 'EFF.EBPP.2639',
                 'JBO.DFB.1400', 'JBO.DFB.1520', 'JBO.DFB.5000',
                 'NRT.BON.1400', 'NRT.BON.1600', 'NRT.BON.2000',
                 'WSRT.P1.328', 'WSRT.P1.328.C', 'WSRT.P1.323.C',
                 'WSRT.P1.382', 'WSRT.P1.382.C', 'WSRT.P1.367.C',
                 'WSRT.P1.840', 'WSRT.P1.840.C',
                 'WSRT.P1.1380', 'WSRT.P1.1380.C',
                 'WSRT.P1.1380.1',
                 'WSRT.P1.1380.2', 'WSRT.P1.1380.2.C',
                 'WSRT.P1.2273.C',
                ]
        }  # list of all non-legacy backends (is this complete?)

psrdict = {}
for p in psr_names:
    psrdict[p] = filt

dr2u.make_dataset(psrdict, indir='data/partim_classic',
                  outdir='data/partim_cut_IPTA', tmin=2, min_toas=10,
                  frequency_filter=False, bw=0.0, dt=100)


datadir = 'data/partim_cut_IPTA'
parfiles = glob.glob('data/partim_cut_IPTA/*.par')
psrlist = []
for p in parfiles:
    name = p.split('/')[-1]
    psrlist.append(name.split('.')[0])
psrlist.sort()

list_file = os.path.join(datadir, 'psrlist_cut_IPTA.txt')
with open(list_file, 'w') as f:
    for pname in psrlist:
        f.write("{:s}\n".format(pname))
