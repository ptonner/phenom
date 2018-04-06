import numpy as np
import pandas as pd
import patsy
import re
from copy import copy

from itertools import product


class Design(object):

    """Abstract class defining the interface for various designs.

    Designs take as input the metadata for a given dataset, and construct a
    design matrix for use by phenom. The design also includes information about
    how latex functional variables are grouped by different priors, names for
    each generated functional variable, and properties of the design space (e.g.
    number of variables and samples)."""

    def __init__(self, name='design', *args, **kwargs):
        self.name = name
        
    def __call__(self, meta):
        return self.process(meta)

    def process(self, meta):
        raise NotImplemented()

    def matrix(self, meta):
        """The design matrix."""
        return self._matrix(meta)

    def _matrix(self, meta):
        raise NotImplemented('no matrix defined for this design!')

    @property
    def priors(self):
        """The priors for each variable in the design."""
        return np.array(self._priors()).astype(int)

    def _priors(self):
        raise NotImplemented('priors are not defined for this design!')

    @property
    def k(self):
        """The number of variables in design (e.g. columns in the design matrix)."""
        return self.matrix.shape[1]

    @property
    def n(self):
        """The number of samples in the design."""
        return self.meta.shape[0]

    @property
    def L(self):
        """The number of priors needed for the design."""
        return max(self.priors) + 1

    @property
    def names(self):
        """The names for each variable in the design."""
        return self._names()

    def _names(self):
        raise NotImplemented()

    @property
    def frame(self):
        """A dataframe of the design matrix."""
        return pd.DataFrame(self.matrix, columns=self.names)

    def __mul__(self, other):
        """multiple two designs"""
        return Kron(self, other)

    def __add__(self, other):
        """add two designs"""
        return Add(self, other)


class Add(Design):

    def __init__(self, d1, d2, name='add', *args, **kwargs):
        super(Add, self).__init__(name, *args, **kwargs)

        self.d1 = d1
        self.d2 = d2

        if self.d2.name == self.d1.name:
            self.d2.name += '_2'

    def process(self, meta):
        dm1, priors1, ap1, lsp1 = self.d1.process(meta)
        dm2, priors2, ap1, lsp2 = self.d2.process(meta)

        dm = np.concatenate((dm1, dm2), 1)
        priors = priors1 + priors2
        alphaPriors = ap1 + ap2
        lengthscalePriors = ls1 + ls2

        return dm, priors, alphaPriors, lengthscalePriors

    def _k(self):
        return self.d1.k + self.d2.k

    def _names(self):
        return self.d1.names + self.d2.names

    def _matrix(self):
        return np.concatenate((self.d1.matrix, self.d2.matrix), 1)

    def _priors(self):

        p1 = self.d1.priors
        p2 = self.d2.priors

        return p1.tolist() + (p2 + max(p1) + 1).tolist()


class Kron(Design):

    def __init__(self, d1, d2, name='kron', combineMethod = 'average', *args, **kwargs):
        super(Kron, self).__init__(name, *args, **kwargs)
        self.d1 = d1
        self.d2 = d2
        self.combineMethod = combineMethod

    def process(self, meta):

        dm1, priors1, ap1, lsp1 = self.d1.process(meta)
        dm2, priors2, ap1, lsp2 = self.d2.process(meta)
        
        dm = []
        n = dm1.shape[0]
        for i in range(n):
            dm.append(np.kron(dm1[i,:], m2[i,:]))
        dm = np.array(dm)

        priors = []
        for v1, v2 in product(priors1, priors2):
            priors.append(v1 + v2 * (max(priors1) + 1))

        def combine(l1, l2, method):
            ret = []
            for e1, e2 in product(l1, l2):
                if method == 'average':
                    ret.append([(e1[0] + e2[0]) / 2, (e1[1] + e2[1]) / 2])

        alphaPriors = combine(ap1, ap2, self.combineMethod)
        lengthscalePriors = combine(lsp1, lsp2, self.combineMethod)
        
        return dm, priors, alphaPriors, lengthscalePriors

    def _matrix(self):

        ret = []
        m1, m2 = self.d1.matrix, self.d2.matrix

        for i in range(self.n):
            ret.append(np.kron(m1[i,:], m2[i,:]))

        return np.array(ret)

    def _priors(self):

        p1 = self.d1.priors
        p2 = self.d2.priors

        priors = []

        for v1, v2 in product(p1, p2):
            priors.append(v1 + v2 * (max(p1) + 1))

        return priors

    def _names(self):

        return ['(%s: %s)x(%s: %s)' % (self.d1.name, v1, self.d2.name, v2) for v1, v2 in product(self.d1.names, self.d2.names)]

