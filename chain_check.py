#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--burn', dest='burn', action='store',
                    type=int, default='1000', help='Chain Directory')
parser.add_argument('--chaindir', dest='chaindir', action='store',
                    type=str, default='./', help='Chain Directory')

args = parser.parse_args()

chains = np.loadtxt(args.chaindir+'/chain_1.txt')

fig = plt.figure(figsize=[10,5])
upper = 10**np.percentile(chains[args.burn:, -5], q=95)
plt.subplot(121)
plt.title('Chain Iterations: {0}'.format(chains.shape[0]))
plt.plot(chains[:,-5])

plt.subplot(122)
plt.hist(chains[args.burn:,-5], bins=50, histtype='step')
plt.axvline(np.log10(upper),linestyle='--')
plt.title('Upper Limit= {0:1.2e}, Burn:{1}'.format(upper,args.burn))

plt.tight_layout()
plt.show()
