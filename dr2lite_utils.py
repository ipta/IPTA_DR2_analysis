# dr2lite_utils.py
"""utilities for manipulating IPTA .par and .tim files to construct DR2-lite datasets
"""

from __future__ import division, print_function

import os, re

import numpy as np
import matplotlib.pyplot as plt

import libstempo as t2


# parameters to cut from .par files
_cut = [
    'T2EFAC', 'T2EQUAD', 'ECORR',
    'TNEF', 'TNEQ', 'TNECORR',
    'TNDMAmp', 'TNDMGam', 'TNDMC', 'TNSubtractDM',
    'TNRedAmp', 'TNRedGam', 'TNRedC',
    'DMMODEL', '_DM', '_CM', 'CONSTRAIN', 'DMOFF',
    'START', 'FINISH',
    'TZRSITE', 'TZRMJD', 'TZRFRQ',
]

def clean_par(infile, outfile):
    """setup parfile for GW analysis

    Remove tempo fitted noise and dm parameters.  Add DM1, DM2.

    :param infile: input .par file
    :param outfile: output .par file
    """
    if infile == outfile:
        raise ValueError("outfile must be different than infile")

    with open(infile, 'r') as fin:
        lines = fin.readlines()

    hasDM1 = False
    if any([line.startswith('DM1') for line in lines]): hasDM1 = True

    with open(outfile, 'w') as fout:
        for line in lines:
            if not any([line.startswith(flag) for flag in _cut]):
                fout.write('%s'%line)
                try:
                    if not hasDM1 and line.split()[0] == 'DM':
                        fout.write('DM1 0 1\n')
                        fout.write('DM2 0 1\n')
                except IndexError:
                    pass

def combine_tim(infile, outfile):
    """combine all .tim files into one super .tim file
    
    :param infiles: DR2 .tim file that links to other actual .tim files
    :param outfile: output .tim file
    """
    if infile == outfile:
        raise ValueError("outfile must be different than infile")

    psr_dir = os.path.dirname(infile)

    # get all included .tim files
    tims = []
    with open(infile, 'r') as fin:
        for line in fin:
            if line.startswith('INCLUDE'):
                this_tim = line.strip().split()[-1]
                tims.append(os.path.join(psr_dir, this_tim))

    # write contents of all .tim to common outfile
    with open(outfile, 'w') as fout:
        fout.write('FORMAT 1\n')
        fout.write('MODE 1\n')
        for t in tims:
            with open(t, 'r') as fin:
                for line in fin:
                    if not (line.startswith('FORMAT') or
                            line.startswith('MODE') or
                            line.startswith('C')):
                        # ensure every line ends in newline char
                        line = line.strip() + '\n'
                        fout.write(line)


def remove_addsat(timfile):
    """remove `-addsat` flags from TOAs

    When .tim files are read by `tempo`/`libstempo`, it applies site
    arrival time offsets directly to the TOAs based on the `addsat` flags.
    If .tim files are rewritten to disk the modified TOAs and the flags
    are written.  When that file is read, the `addsat`s are applied again.
    To avoid this we must manually remove `-addsat` flags after we write
    .tim files using `libstempo.tempo2pulsar.savetim()`.

    This is a dangerous operation! NEVER modify the base .tim file in the
    DR2 directory!
    """
    pattern = re.compile(r" \-addsat [+-]\d+")
    # ' \-addsat ' match -addsat and spaces around it
    # '[+-]' match plus or minus
    # '\d+' match one or more digits (addsats are integer seconds)

    with open(timfile, 'r') as f:
        tim_lines = f.read()

    tim_lines = pattern.sub("", tim_lines)

    with open(timfile, 'w') as f:
        f.write(tim_lines)

def fix_jumps(psr, verbose=True):
    """set new reference backend for jumps if original reference is filtered

    After a new reference is set that backend will have JUMP=0.  All other
    jumps are refit using libstempo's default fitter.

    :param psr: libstempo.tempopulsar object
    :param verbose: boolean flag.  If True, excecution notes are printed
    """

    # Don't use deleted points
    mask = psr.deleted[:] == 0

    # find JUMP parameter indices
    idx = np.array([1+ct for ct, p in enumerate(psr.pars()) if 'JUMP' in p])

    # get list of TOAs in each jump
    M = psr.designmatrix()[mask, :]
    jbool = [(M[:, ix]!=0.0).astype(int) for ix in idx]

    # check if all TOAs are jumped
    if np.sum(np.sum(jbool, axis=0) != 0) >= len(psr.toas()[mask]):
        if verbose: print('All TOAs are being jumped!')

        # find JUMP group with lowest weighted variance
        maxind = np.argmax([(1/psr.toaerrs[mask][np.flatnonzero(M[:, ix])]**2).sum() for ix in idx])
        jpar = psr.pars()[idx[maxind]-1]
        if verbose: print('Setting {} as reference jump.'.format(jpar))
        psr[jpar].fit = False
        psr[jpar].val = 0

        # run libstempo fitter to refit jumps relative to new reference
        try:
            _ = psr.fit()
        except np.linalg.LinAlgError:
            print("LinAlgError in libstempo.tempopulsar.fit(), skipping refit")


def get_dm_bins(toas, dt=7):
    """determine which toas belong in each DM bin

    :param toas: shape (1,) array of TOAs
    :param dt: size of DM bin (days)

    :return: list of boolean arrays for each time bin corresponding
             to TOAs within that bin
    """
    bins = int(np.ceil((toas.max() - toas.min()) / (86400*dt)))
    tmin = toas.min() - 1
    tmax = toas.max() + 1
    _, xedges = np.histogram(toas, bins=bins, range=[tmin, tmax])
    return [np.logical_and(toas >= xedges[ct], toas <= xedges[ct+1])
            for ct in range(len(xedges)-1)]


def filter_psr(psr, bw=1.1, dt=7, filter_dict=None, min_toas=10,
               frequency_filter=True, fmax=3000, verbose=True, plot=False):
    """apply frequency coverage, PTA, and/or other flag filters to pulsar

    :param psr:
        libstempo.tempopulsar object
    :param bw:
        min bandwidth ratio fmax/fmin for DM estimation
    :param dt:
        DM bin width (days)
    :param filter_dict:
        dictionary of filters to apply
        The dictionary has keys that correspond to TOA flags. The values are
        a single or list of acceptable flagvals, e.g.

        {'pta':'NANOGrav', 'group':['Rcvr_800_GUPPI', 'Rcvr1_2_GUPPI']}

        or

        {'pta':['PPTA', 'EPTA']}

    :param min_toas:
        integer, minimum number of TOAs for a given backend before dropping
    :param frequency_filter:
        boolean flag, if true apply frequency coverage filter
    :param fmax:
        high frequency threshold (MHz) to ignore filter
        If TOA freq is greater than this, always keep
    :param verbose:
        boolean flag.  If True, print notes.
    :param plot:
        boolean flag.  If True, generate diagnostic plots

    :return: libstempo.tempopulsar object with filters applied
    """
    psr.deleted[:] = 1
    print('Working on PSR {}'.format(psr.name))

    # Flag filtering
    flag_keep = []
    if filter_dict:
        for key, val in filter_dict.items():
            if verbose: print('Keeping TOAs corresponding to {} {}'
                              .format(key, val))
            if type(val) is not list:
                val = [val]
            flag_conds = [psr.flagvals(key)==v for v in val]
            # if TOA has ANY acceptable value for this flag
            flag_keep.append(np.any(flag_conds, axis=0))

    # if TOA satisfies all flags
    idx_flag = np.flatnonzero(np.alltrue(flag_keep, axis=0))

    # filter for frequency coverage
    if frequency_filter:
        if verbose: print("Running multi-frequency filter")
        bins = get_dm_bins(psr.toas()*86400, dt=dt)
        idx_freq = []
        for bn in bins:
            if sum(bn) > 1:
                ix = list(filter(lambda x: x in idx_flag, np.flatnonzero(bn)))
                if len(ix) > 0:
                    if psr.freqs[ix].max() / psr.freqs[ix].min() >= bw:
                        idx_freq.append(ix)
                    elif psr.freqs[ix].max() >= fmax:
                        idx_freq.append(ix)

        # check for empty list (i.e. there is no multi-frequency data)
        if not idx_freq:
            print("No multi-frequency data, returning original psr")
            return psr

        # delete
        idx = np.unique(np.concatenate(idx_freq))
    else:
        idx = idx_flag
    psr.deleted[idx] = 0  # mark filtered TOAs as "deleted"

    # check for "orphan" backends (less than min_toas obsv.)
    orphans = []
    for gr in np.unique(psr.flagvals('group')):
        in_group = [gr == b for b in psr.flagvals('group')]
        mask = np.logical_and(in_group, ~psr.deletedmask())
        N = np.sum(mask)
        if N>0 and N<min_toas:
            psr.deleted[mask] = True
            orphans.append([gr, N])
    if verbose: print("backends marked as 'orphan': {}".format(ophans))

    # filter design matrix
    mask = np.logical_not(psr.deleted)
    M = psr.designmatrix()[mask, :]
    dpars = []
    for ct, (par, val) in enumerate(zip(psr.pars(), M.sum(axis=0)[1:])):
        if val == 0:
            dpars.append(par)
            psr[par].fit = False
            psr[par].val = 0.0

    if verbose:
        print('Cutting {} TOAs'.format(np.sum(~mask)))
        if len(dpars): print('Turning off fit for {}'.format(dpars))

    fix_jumps(psr)

    if plot:
        plt.figure(figsize=(8,3))
        for pta in np.unique(psr.flagvals('pta')):
            nix = psr.flagvals('pta') == pta
            plt.plot(psr.toas()[nix], psr.freqs[nix], '.', label=pta)
        plt.plot(psr.toas()[~psr.deletedmask()], psr.freqs[~psr.deletedmask()], '.',
                 color='C3', alpha=0.3, label='filtered')
        plt.legend(loc='best', frameon=False)
        plt.title(psr.name)

    return psr


def make_dataset(psrdict, indir, outdir='partim_filtered',
                 frequency_filter=True, bw=1.1, dt=7, fmax=3000,
                 tmin=0,
                 plot=False, verbose=True):
    """make a filtered, DR2-lite style dataset of .par and .tim files

    :param psrdict:
        dictionary of dictionaries...
        keys are pulsar names, values are filter_dict to use in filter_psr()
    :param indir:
        string, input directory containing original .par and .tim files
    :param outdir:
        string, output directory to write new .par and .tim files
        if this directory exists, it will be overwritten
    :param frequency_filter:
        boolean flag or dict with pulsar name keys and boolean values.
        If True, apply frequency coverage filter to all / that pulsar.
    :param bw:
        min bandwidth ratio fmax/fmin for DM estimation
    :param dt:
        DM bin width (days)
    :param fmax:
        high frequency threshold (MHz) to ignore filter
        If TOA freq is greater than this, always keep
    :param tmin:
        minimum obstervation time to keep (yrs). No .par/.tim file is saved
        for pulsars with shorter observation times.
    :param verbose:
        boolean flag.  If True, print notes.
    :param plot:
        boolean flag.  If True, generate diagnostic plots
    """

    os.system('rm -rf {}'.format(outdir))
    os.system('mkdir -p {}'.format(outdir))
    for pname, filters in sorted(psrdict.items()):
        parfile = os.path.join(indir, '{}.par'.format(pname))
        timfile = os.path.join(indir, '{}.tim'.format(pname))
        psr = t2.tempopulsar(parfile, timfile, maxobs=30000)
        if isinstance(frequency_filter, dict):
            ff = frequency_filter[pname]
        else:
            ff = frequency_filter
        psr = filter_psr(psr, filter_dict=filters, frequency_filter=ff,
                         bw=bw, dt=dt, fmax=fmax,
                         plot=plot, verbose=verbose)
        toas_keep = psr.toas()[~psr.deletedmask()]
        try:
            Tobs = (toas_keep.max() - toas_keep.min())/365.25 # yrs
        except ValueError:
            Tobs = 0
        if Tobs > tmin:
            newpar = os.path.join(outdir, '{}.par'.format(pname))
            newtim = os.path.join(outdir, '{}.tim'.format(pname))
            psr.savepar(newpar)
            psr.savetim(newtim)
            remove_addsat(newtim)
            if verbose:
                print("filtered Tobs = {:.2f} yrs".format(Tobs))
        else:
            if verbose:
                print("skipping PSR {:}, filtered Tobs = {:.2f} yrs"
                      .format(psr.name, Tobs))
        print('\n')
        del psr
