import numpy as np
import pandas as pd
import patsy

from itertools import product

class Design(object):

    def __init__(self, meta, name='design', *args, **kwargs):
        self.meta = meta
        self.name = name

    @property
    def k(self):
        return self._k()

    def _k(self):
        raise NotImplemented()

    @property
    def n(self):
        return self.meta.shape[0]

    @property
    def names(self):
        return self._names()

    def _names(self):
        raise NotImplemented()

    @property
    def matrix(self):
        return self._matrix()

    def _matrix(self):
        raise NotImplemented()

    @property
    def frame(self):
        return pd.DataFrame(self.matrix, columns = self.names)

    def __mul__(self, other):
        return Kron(self, other)

    def __add__(self, other):
        return Add(self, other)

class Add(Design):

    def __init__(self, d1, d2):
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


class Kron(Design):

    def __init__(self, d1, d2):
        self.d1 = d1
        self.d2 = d2

    def _matrix(self):

        ret = []
        m1, m2 = self.d1.matrix, self.d2.matrix

        for i in range(self.n):
            ret.append(np.kron(m1[i, :], m2[i, :]))

        return np.array(ret)

    def _names(self):

        return ['(%s: %s)x(%s: %s)'%(self.d1.name, v2, self.d2.name, v2) for v1, v2 in product(self.d1.names, self.d2.names)]

class Formula(Design):

    def __init__(self, meta, form, *args, **kwargs):
        super(Formula, self).__init__(meta, *args, **kwargs)
        self.form = form
        self.d = patsy.dmatrix(self.form, self.meta)

    def _matrix(self):
        return np.array(self.d)

    def _names(self):
        return self.d.design_info.column_names
