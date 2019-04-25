from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os, json
import numpy as np

def make_noise_file(psrname, chain, pars, outdir='partim'):
    x = {}
    for ct, par in enumerate(pars):
        x[par] = np.median(chain[:, ct])

    out_file = os.path.join(outdir, '{:s}_noise.json'.format(psrname))
    with open(out_file, 'w') as fout:
        json.dump(x, fout, sort_keys=True, indent=4, separators=(',', ': '))


### ARG PARSER
parser = argparse.ArgumentParser(
          description='run noise analysis with enterprise')

parser.add_argument('-p', '--psr',
                    dest='psr_name', default=None,
                    action='store',
                    help="pulsar to analyze")

parser.add_argument('-d', '--datadir',
                    dest='datadir', default='~/data/',
                    action='store',
                    help="location of par/tim files")

parser.add_argument('-o', '--rundir',
                    dest='rundir', default='~/ipta/noiseruns/{psr:s}/',
                    action='store',
                    help="location of chain")

args = parser.parse_args()

rundir = os.path.abspath(args.rundir)
datadir = os.path.abspath(args.datadir)


burnfrac = 0.25
thin = 10

param_file = os.path.join(args.rundir, 'params.txt')
with open(param_file, 'r') as f:
    params = [line.rstrip() for line in f]

chain_file = os.path.join(args.rundir, 'chain_1.txt')
chain_raw = np.loadtxt(chain_file)

burn = int(burnfrac * len(chain_raw))
chain = chain_raw[burn::thin]

make_noise_file(psr, chain, params, outdir=datadir)
