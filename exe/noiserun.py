from __future__ import (absolute_import, division,
                        print_function) # , unicode_literals)
import numpy as np
import os, argparse

from enterprise_extensions import models
from enterprise_extensions import model_utils

from enterprise.pulsar import Pulsar
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

### ARG PARSER
parser = argparse.ArgumentParser(
          description='run noise analysis with enterprise')

parser.add_argument('-p', '--psr',
                    dest='psr_name', default=None,
                    action='store',
                    help="pulsar to analyze")

parser.add_argument('-e', '--ephem',
                    dest='ephem', default='DE436',
                    action='store',
                    help="JPL ephemeris version to use")

parser.add_argument('-f', '--freespec',
                    dest='fs', default=False,
                    action='store_true',
                    help="use free spectral RN model")

parser.add_argument('-d', '--datadir',
                    dest='datadir', default='~/data/',
                    action='store',
                    help="location of par/tim files")

parser.add_argument('-o', '--outdir',
                    dest='outdir', default='~/ipta/noiseruns/{psr:s}/',
                    action='store',
                    help="location to write output")

parser.add_argument('-t', '--thin', type=int,
                    dest='thin', default=10,
                    action='store',
                    help="thinning factor (keep every [thin]th sample)")

parser.add_argument('-N', '--Nsamp', type=int,
                    dest='N', default=int(5.0e+04),
                    action='store',
                    help="number of samples to collect (after thinning!!)")

args = parser.parse_args()

outdir = os.path.abspath(args.outdir)
os.system('mkdir -p {}'.format(outdir))

# adjust Nsamp for existing chain
chfile = os.path.join(outdir, 'chain_1.txt')
if os.path.exists(chfile):
    ct = sum(1 for i in open(chfile, 'rb'))
    if ct >= args.N:
        print("{:s} has {:d} samples... exiting".format(chfile, ct))
        exit(0)
    else:
        args.N -= ct

# read in data from .par / .tim
par = os.path.join(args.datadir, '{}.par'.format(args.psr_name))
tim = os.path.join(args.datadir, '{}.tim'.format(args.psr_name))
psr = Pulsar(par, tim, ephem=args.ephem, timing_package='tempo2')

# use dm_expdip for J1713+0747
use_dmdip=False
if 'J1713+0747' == psr.name:
    use_dmdip = True

#################
##  pta model  ##
#################
red_psd = 'spectrum' if args.fs else 'powerlaw'
pta = models.model_singlepsr_noise(
        psr,
        tm_svd=True,
        white_vary=True,
        red_var=True, psd=red_psd, components=30,
        dm_var=True, dm_type='gp', dmgp_kernel='diag', dm_psd='powerlaw',
        dm_annual=True, dm_expdip=use_dmdip,
        upper_limit=False
)

outfile = os.path.join(args.outdir, 'params.txt')
with open(outfile, 'w') as f:
    for pname in pta.param_names:
        f.write(pname+'\n')

###############
##  sampler  ##
###############
sampler = model_utils.setup_sampler(pta, outdir=args.outdir, resume=True)

x0 = np.hstack(p.sample() for p in pta.params)  # initial point

# SAMPLE!!
Nsamp = args.N * args.thin
sampler.sample(x0, Nsamp,
               SCAMweight=30, AMweight=20, DEweight=50,
               burn=int(2e4), thin=args.thin)
