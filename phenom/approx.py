import attr
import GPy
import numpy as np


@attr.s
class Approx():

    phenotype = attr.ib()
    x_kernel = attr.ib(default=GPy.kern.RBF)

    def kernels(self):
        kx = self.x_kernel(1)

        kdesign = None
        for i in range(self.phenotype.design.L):

            # select indices that belong to this prior
            ind, = np.where(self.phenotype.design.priors == i)

            if kdesign is None:
                kdesign = GPy.kern.Linear(
                    ind.shape[0], active_dims=ind, ARD=False)
            else:
                ktmp = GPy.kern.Linear(
                    ind.shape[0], active_dims=ind, ARD=False)
                kdesign = kdesign + ktmp

        return kx, kdesign

    def model(self, optimize=True):
        kx, kdesign = self.kernels()

        cfg = self.phenotype.config()

        m = GPy.models.GPKroneckerGaussianRegression(
            cfg['x'][:, None], cfg['design'], cfg['y'].T, kx, kdesign
        )

        if optimize:
            m.optimize()

        return m

    def posterior(self, m, p):
        design = self.phenotype.design
        priors = design.priors

        ind, = np.where(priors == p)

        u = np.unique(design.matrix[:, priors == p], axis=0)

        mus = []
        vars = []

        for i in range(u.shape[0]):
            predx2 = np.zeros((1, design.matrix.shape[1]))
            predx2[0, ind] = u[i, :]

            mu, var = m.predict(m.X1, predx2)

            mus.append(mu[:, 0])
            vars.append(var[:, 0])

        return np.column_stack(mus), np.column_stack(var)
