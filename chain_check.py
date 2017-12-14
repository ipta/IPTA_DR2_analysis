#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--burn', dest='burn', action='store',
                    type=int, default='1000', help='Chain Directory')
parser.add_argument('--chaindir', dest='chaindir', action='store',
                    type=str, default='./', help='Chain Directory')
parser.add_argument('--saveplot', dest='saveplot', action='store_true',
                    help='Save the plot in Chain Directory')
parser.add_argument('--showplot', dest='saveplot', action='store_false',
                    help='Show the plot, do not save it')

args = parser.parse_args()

chains = np.loadtxt(args.chaindir+'/chain_1.txt')

fig = plt.figure(figsize=[10,5])
upper = 10**np.percentile(chains[args.burn:, -5], q=95)
plt.subplot(121)
plt.title('Chain Iterations: {0}'.format(chains.shape[0]))
plt.plot(chains[:,-5], zorder=0, label='chain')
plt.vlines(args.burn, min(chains[:, -5]), max(chains[:, -5]), linestyle='--', zorder=1, label='burn')
plt.legend(loc='best')

plt.subplot(122)
plt.hist(chains[args.burn:,-5], bins=50, histtype='step')
plt.axvline(np.log10(upper),linestyle='--')
plt.title('Upper Limit= {0:1.2e}, Burn:{1}'.format(upper, args.burn))

plt.tight_layout()
if args.saveplot:
    plt.savefig(args.chaindir + '/chaincheck.pdf')
else:
    plt.show()
