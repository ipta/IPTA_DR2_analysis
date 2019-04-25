# sampling utilities for IPTA DR2 runs
#  at some point should merge into enterprise extensions

import numpy as np

class UserDraw(object):
    """object for user specified proposal distributions
    """
    def __init__(self, idxs, samplers, log_qs=None, name=None):
        """
        :param idxs: list of parameter indices to use for this jump
        :param samplers: dict of callable samplers
            keys should include all idxs
        :param lqxys: dict of callable log proposal distributions
            keys should include all idxs
            for symmetric proposals set `log_qs=None`, then `log_qxy=0`
        :param name: name for PTMCMC bookkeeping
        """
        #TODO check all idxs in keys!
        self.idxs = idxs
        self.samplers = samplers
        self.log_qs = log_qs

        if name is None:
            namestr = 'draw'
            for ii in samplers.keys():
                namestr += '_{}'.format(ii)
            self.__name__ = namestr
        else:
            self.__name__ = name

    def __call__(self, x, iter, beta):
        """proposal from parameter prior distribution
        """
        y = x.copy()

        # draw parameter from idxs
        ii = np.random.choice(self.idxs)

        try: # vector parameter
            y[ii] = self.samplers[ii]()[0]
        except (IndexError, TypeError) as e:
            y[ii] = self.samplers[ii]()

        if self.log_qs is None:
            lqxy = 0
        else:
            lqxy = self.log_qs[ii](x[ii]) - self.log_qs[ii](y[ii])

        return y, lqxy


def build_prior_draw(pta, parlist, name=None):
    """create a callable object to perfom a prior draw
    :param pta:
        instantiated PTA object
    :param parlist:
        single string or list of strings of parameter name(s) to
        use for this jump.
    :param name:
        display name for PTMCMCSampler bookkeeping
    """
    if not isinstance(parlist, list):
        parlist = [parlist]
    idxs = [pta.param_names.index(par) for par in parlist]

    # parameter map
    pmap = []
    ct = 0
    for ii, pp in enumerate(pta.params):
        size = pp.size or 1
        for nn in range(size):
            pmap.append(ii)
        ct += size

    sampler = {ii: pta.params[pmap[ii]].sample for ii in idxs}
    log_q = {ii: pta.params[pmap[ii]].get_logpdf for ii in idxs}

    return UserDraw(idxs, sampler, log_q, name=name)

def build_loguni_draw(pta, parlist, bounds, name=None):
    """create a callable object to perfom a log-uniform draw
    :param pta:
        instantiated PTA object
    :param parlist:
        single string or list of strings of parameter name(s) to
        use for this jump.
    :param bounds:
        tuple of (pmin, pmax) for draw
    :param name:
        display name for PTMCMCSampler bookkeeping
    """
    if not isinstance(parlist, list):
        parlist = [parlist]
    idxs = [pta.param_names.index(par) for par in parlist]
    pmin, pmax = bounds
    sampler = {ii: (lambda : np.random.uniform(pmin,pmax)) for ii in idxs}

    return UserDraw(idxs, sampler, None, name=name)
