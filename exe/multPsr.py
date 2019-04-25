from __future__ import (absolute_import, division,
                        print_function) # , unicode_literals)
import numpy as np
import os, argparse, json

from enterprise_extensions import models
from enterprise_extensions import model_utils as modu

from enterprise.pulsar import Pulsar
from enterprise.signals import parameter
from enterprise.signals import selections
from enterprise.signals import utils
from enterprise.signals import signal_base
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import deterministic_signals
import enterprise.constants as const

from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
import sample_utils as su


### ARG PARSER
parser = argparse.ArgumentParser(
          description='run single pulsar UL analysis with enterprise')

parser.add_argument('-p', '--psrlist',
                    dest='psrlist', default=None,
                    action='store',
                    help="list of pulsars to analyze")

parser.add_argument('-e', '--ephem',
                    dest='ephem', default='DE436',
                    action='store',
                    help="JPL ephemeris version to use")

parser.add_argument('-b', '--bayesephem',
                    dest='be', default=False,
                    action='store_true',
                    help="use BayesEphem SSE uncertainty model")

parser.add_argument('-u', '--upperlimit',
                    dest='ul', default=False,
                    action='store_true',
                    help="use LinearExp amplitude priors for UL")

parser.add_argument('-d', '--datadir',
                    dest='datadir', default='~/data/',
                    action='store',
                    help="location of par/tim files")

parser.add_argument('-o', '--outdir',
                    dest='outdir', default='~/ipta/gwbruns/',
                    action='store',
                    help="location to write output")

parser.add_argument('-N', '--Nsamp', type=int,
                    dest='N', default=int(1.0e+05),
                    action='store',
                    help="number of samples to collect (after thinning!!)")

parser.add_argument('-t', '--thin', type=int,
                    dest='thin', default=10,
                    action='store',
                    help="thinning factor (keep every [thin]th sample)")

parser.add_argument('-y', '--dm1yr',
                    dest='use_dm1yr', default=False,
                    action='store_true',
                    help="use the DM yearly sinusoid model")

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

# read in pulsar data
with open(args.psrlist, 'r') as f:
    psrlist = [line.strip() for line in f]

psrs = []
noise_params = {}
for psr_name in psrlist:
    pfile = os.path.join(args.datadir, '{}.par'.format(psr_name))
    tfile = os.path.join(args.datadir, '{}.tim'.format(psr_name))
    nfile = os.path.join(args.datadir, "{}_noise.json".format(psr_name))
    psrs.append(Pulsar(pfile, tfile, ephem=args.ephem, timing_package='tempo2'))
    with open(nfile, 'r') as f:
        noise_params.update(json.load(f))

# find the maximum time span to set GW frequency sampling
tmin = [p.toas.min() for p in psrs]
tmax = [p.toas.max() for p in psrs]
Tspan = np.max(tmax) - np.min(tmin)

######################
##  ipta DR2 model  ##
######################
bkend = selections.Selection(selections.by_backend)
bkend_NG = selections.Selection(selections.nanograv_backends)

# fix white noise parameters
efac = parameter.Constant()
equad = parameter.Constant()
ecorr = parameter.Constant()

# DMGP parameters and powerlaw
dm_log10_A = parameter.Uniform(-20, -11)
dm_gamma = parameter.Uniform(0, 7)
dm_pl = utils.powerlaw(log10_A=dm_log10_A, gamma=dm_gamma)
dm_basis = utils.createfourierdesignmatrix_dm(nmodes=30, Tspan=Tspan)

# red noise parameters and powerlaw
if args.ul:
    rn_log10_A = parameter.LinearExp(-20, -11)
else:
    rn_log10_A = parameter.Uniform(-20, -11)
rn_gamma = parameter.Uniform(0, 7)
rn_pl = utils.powerlaw(log10_A=rn_log10_A, gamma=rn_gamma)

# GWB parameters and powerlaw
orf = utils.hd_orf()
if args.ul:
    gwb_log10_A = parameter.LinearExp(-18, -12)('gwb_log10_A')
else:
    gwb_log10_A = parameter.Uniform(-18, -12)('gwb_log10_A')
gwb_gamma = parameter.Constant(13/3)('gwb_gamma')
gwb_pl = utils.powerlaw(log10_A=gwb_log10_A, gamma=gwb_gamma)

# signals
ef = white_signals.MeasurementNoise(efac=efac, selection=bkend)
eq = white_signals.EquadNoise(log10_equad=equad, selection=bkend)
ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=bkend_NG)

dmgp = gp_signals.BasisGP(dm_pl, dm_basis, name='dm_gp')
dmexp = models.dm_exponential_dip(tmin=54500, tmax=55000, name='dmexp_1')
dm1yr = models.dm_annual_signal()

rn = gp_signals.FourierBasisGP(spectrum=rn_pl, components=30, Tspan=Tspan)
gwb = gp_signals.FourierBasisCommonGP(gwb_pl, orf, components=30, Tspan=Tspan,
                                      name='gw')

tm = gp_signals.TimingModel(use_svd=True)
be = deterministic_signals.PhysicalEphemerisSignal(  # widen prior on jup orbit
    use_epoch_toas=True,
    jup_orb_elements=parameter.Uniform(-0.1,0.1,size=6)('jup_orb_elements')
)

# full model
mod_root = tm + ef + eq + dmgp + rn + gwb
if args.be: mod_root += be
if args.use_dm1yr: mod_root += dm1yr
mod_J1713 = mod_root + dmexp
mod_ecorr = mod_root + ec

mods = []
for psr in psrs:
    if 'J1713+0747' == psr.name:
        if 'NANOGrav' in psr.flags['pta']:
            mod_J1713 += ec
        mods.append(mod_J1713(psr))
    elif 'NANOGrav' in psr.flags['pta']:
        mods.append(mod_ecorr(psr))
    else:
        mods.append(mod_root(psr))

pta = signal_base.PTA(mods)
pta.set_default_params(noise_params)

outfile = os.path.join(args.outdir, 'params.txt')
with open(outfile, 'w') as f:
    for pname in pta.param_names:
        f.write(pname+'\n')

###############
##  sampler  ##
###############
x0 = np.hstack([noise_params[p.name] if p.name in noise_params.keys()
                else p.sample() for p in pta.params])  # initial point
ndim = len(x0)

# set initial cov stdev to (starting order of magnitude)/10
stdev = np.array([10**np.floor(np.log10(abs(x)))/10 for x in x0])
cov = np.diag(stdev**2)

# generate custom sampling groups
groups = [list(range(ndim))]

# pulsar noise groups (RN + DM)
for psr in psrs:
    this_group = [pta.param_names.index(par)
                  for par in pta.param_names if psr.name in par]
    groups.append(this_group)

groups.append([pta.param_names.index('gwb_log10_A')])

if args.be:
    # all BE params
    this_group = [pta.param_names.index(par)
                  for par in pta.param_names
                  if 'jup_orb' in par or 'mass' in par or 'frame_drift' in par]
    groups.append(this_group)

    # jup_orb elements + GWB
    this_group = [pta.param_names.index(par)
                  for par in pta.param_names if 'jup_orb' in par]
    this_group.append(pta.param_names.index('gwb_log10_A'))
    groups.append(this_group)

sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, groups=groups,
                 outDir=outdir, resume=True)

# additional proposals
full_prior = su.build_prior_draw(pta, pta.param_names, name='full_prior')
sampler.addProposalToCycle(full_prior, 10)

RNA_params = [pname for pname in pta.param_names if 'red_noise_log10_A' in pname]
RN_loguni = su.build_loguni_draw(pta, RNA_params, (-20,-11), name='RN_loguni')
sampler.addProposalToCycle(RN_loguni, 5)

GWB_loguni = su.build_loguni_draw(pta, 'gwb_log10_A', (-18,-12), name='GWB_loguni')
sampler.addProposalToCycle(GWB_loguni, 5)

if args.be:
    BE_params = [pta.param_names[ii] for ii in groups[6]]
    BE_prior = su.build_prior_draw(pta, BE_params, name='BE_prior')
    sampler.addProposalToCycle(BE_prior, 5)

# SAMPLE!!
Nsamp = args.N * args.thin
sampler.sample(x0, Nsamp,
               SCAMweight=30, AMweight=20, DEweight=50,
               burn=int(5e4), thin=args.thin)
