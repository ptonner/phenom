import patsy
import re
import attr

import numpy as np
import pandas as pd
from copy import copy
from abc import ABCMeta, abstractmethod
from itertools import product

@attr.s
class Design(metaclass=ABCMeta):

    meta = attr.ib()

    @property
    def k(self):
        return self.matrix.shape[1]

    @property
    def n(self):
        return self.meta.shape[0]

    @property
    def L(self):
        return max(self.priors) + 1

    @property
    def names(self):
        return self._names()

    @abstractmethod
    def _names(self):
        pass

    @property
    def matrix(self):
        return self._matrix()

    @abstractmethod
    def _matrix(self):
        pass

    @property
    def frame(self):
        return pd.DataFrame(self.matrix, columns=self.names)

    @property
    def priors(self):
        return np.array(self._priors()).astype(int)

    @abstractmethod
    def _priors(self):
        pass

    def __mul__(self, other):
        return Kron(self.meta, self, other)

    def __add__(self, other):
        return Add(self.meta, self, other)


@attr.s
class Add(Design):

    d1 = attr.ib()
    d2 = attr.ib()

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


@attr.s
class Kron(Design):

    d1 = attr.ib()
    d2 = attr.ib()

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

        # return ['(%s: %s)x(%s: %s)' % (self.d1.name, v1, self.d2.name, v2) for v1, v2 in product(self.d1.names, self.d2.names)]
        return ['(%s)x(%s)' % (v1, v2) for v1, v2 in product(self.d1.names, self.d2.names)]


@attr.s
class Formula(Design):

    form = attr.ib()

    def __attrs_post_init__(self):
        self.d = patsy.dmatrix(self.form, self.meta)

    def _matrix(self):
        return np.array(self.d)

    def _names(self):
        cols = copy(self.d.design_info.column_names)

        # convert all categorical factors to something pretty
        pat = 'C\((?P<factor>[a-zA-Z0-9]+)(?P<ignore>, Treatment\([a-zA-Z0-9.]+\))?\)\[T?\.?(?P<level>[a-zA-Z0-9. ]+)\]'
        comp = re.compile(pat)
        found = map(comp.findall, cols)

        for i, f in enumerate(found): #range(len(cols)):
            if len(f) > 0:
                cols[i] = ', '.join(['%s=%s' % (ff, l) for ff, _, l in f])

        return cols

    def _priors(self):
        priors = -1 * np.ones(self.k)

        for t in self.d.design_info.term_names:
            priors[self.d.design_info.term_name_slices[t]] = max(priors) + 1

        return priors
