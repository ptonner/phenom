import numpy as np
import pystan
import os
import pickle
import pandas as pd


class Phenotype(object):

    def __init__(self, data, design, model='phenom.stan',
                 alpha_priors=None, lengthscale_priors=None, sigma_prior=[1, 1],
                 maxExpectedCross=100, minExpectedCross=.01):

        self.data = data
        self.design = design

        self.alpha_priors = alpha_priors
        if self.alpha_priors is None:
            self.alpha_priors = [[1, 1]] * self.design.L

        self.lengthscale_priors = lengthscale_priors
        if self.lengthscale_priors is None:
            self.lengthscale_priors = [[1, 1]] * self.design.L

        self.sigma_prior = sigma_prior

        self.maxExpectedCross = maxExpectedCross
        self.minExpectedCross = minExpectedCross

        self.posteriors = []
        self._model = None
        self._modelFile = model

    @property
    def model(self):
        d, _ = os.path.split(__file__)
        
        if self._model is None:
            self._model = pystan.StanModel(
                file = os.path.join(d, 'stan', self._modelFile)
            )
        return self._model

    def config(self):

        x, y, _, _ = self.normalize()
        dm = self.design.matrix
        priors = self.design.priors + 1
        k = self.design.k
        L = max(priors)

        cfg = {
            'N': x.shape[0],
            'P': y.shape[1],
            'K': dm.shape[1],
            'L': L,
            'prior': priors,
            'design': dm,
            'x': x,
            'y': y.T,
            'alpha_prior': self.alpha_priors,
            'lengthscale_prior': self.lengthscale_priors,
            'sigma_prior': self.sigma_prior,
        }

        # if self.model in ['mrep', 'mfull']:
        #     cfg['marginal_alpha_prior'] = [
        #         self.marginalEffect_variance_alpha, self.marginalEffect_variance_beta]
        #     cfg['marginal_lengthscale_prior'] = [
        # self.marginalEffect_lengthscale_alpha,
        # self.marginalEffect_lengthscale_beta]

        # expected number of origin crossings = 1/(pi * lengthscale)
        cfg['ls_min'] = 1. / np.pi / self.maxExpectedCross
        cfg['ls_max'] = 1. / np.pi / self.minExpectedCross

        return cfg

    def normalize(self):

        x = self.data.index.values

        xnorm = (x.min(), x.max())
        x = (x - x.min()) / (x.max() - x.min())

        y = self.data.values
        ynorm = (y.mean(), y.std())
        y = (y - y.mean()) / y.std()

        return x, y, xnorm, ynorm

    def save(self, d):
        if not os.path.exists(os.path.join(d, 'samples')):
            os.makedirs(os.path.join(d, 'samples'))

        self.data.to_csv(os.path.join(d, 'data.csv'))
        self.design.meta.to_csv(os.path.join(d, 'meta.csv'))
        self.design.frame.to_csv(os.path.join(d, 'design.csv'))

        x, y, xnorm, ynorm = self.normalize()
        pd.DataFrame(y, index=x).to_csv(os.path.join(d, 'data-normalized.csv'))

        # pickle.dump(self, open(os.path.join(d, 'phenotype.pkl'), 'wb'))

        for i, p in enumerate(self.posteriors):
            samp = p.extract()

            # unnormalized
            x, y, xnorm, ynorm = self.normalize()
            ymean, ystd = ynorm
            xmin, xmax = xnorm

            samp['f-native'] = samp['f'] * ystd
            samp['f-native'][:, 0] += ymean

            samp['df-native'] = samp['df'] * ystd / (xmax - xmin)

            pickle.dump(
                samp, open(os.path.join(d, 'samples', 'posterior_%d.pkl' % i), 'wb'))

            summary = p.summary()
            summary = pd.DataFrame(
                summary['summary'],
                index=summary['summary_rownames'],
                columns=summary['summary_colnames'])
            summary.to_csv(os.path.join(d, 'samples', 'posterior_%d.csv' % i))

    def samples(self, *args, **kwargs):
        cfg = self.config()
        samp = self.model.sampling(data=cfg, *args, **kwargs)

        self.posteriors.append(samp)
        return samp
