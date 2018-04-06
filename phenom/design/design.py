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

    def __init__(self, meta, name='design', *args, **kwargs):
        self.meta = meta
        self.name = name

    @property
    def matrix(self):
        """The design matrix."""
        return self._matrix()

    def _matrix(self):
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
        super(Add, self).__init__(d1.meta, name, *args, **kwargs)

        self.d1 = d1
        self.d2 = d2

        if self.d2.name == self.d1.name:
            self.d2.name += '_2'

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

    def __init__(self, d1, d2, name='kron', *args, **kwargs):
        super(Kron, self).__init__(d1.meta, name, *args, **kwargs)
        self.d1 = d1
        self.d2 = d2

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


class Formula(Design):

    def __init__(self, meta, form, name='formula', *args, **kwargs):
        super(Formula, self).__init__(meta, name, *args, **kwargs)
        self.form = form
        self.d = patsy.dmatrix(self.form, self.meta)

    def _matrix(self):
        return np.array(self.d)

    def _names(self):
        # return self.d.design_info.column_names
        cols = copy(self.d.design_info.column_names)

        # convert all categorical factors to something pretty
        pat = 'C\((?P<factor>[a-zA-Z0-9]+)(?P<ignore>, Treatment\([a-zA-Z0-9.]+\))?\)\[T?\.?(?P<level>[a-zA-Z0-9. ]+)\]'
        comp = re.compile(pat)
        found = map(comp.findall, cols)

        for i in range(len(cols)):
            if len(found[i]) > 0:
                cols[i] = ', '.join(['%s=%s' % (f, l) for f, _, l in found[i]])

        return cols

    def _priors(self):
        priors = -1 * np.ones(self.k)

        for t in self.d.design_info.term_names:
            priors[self.d.design_info.term_name_slices[t]] = max(priors) + 1

        return priors
