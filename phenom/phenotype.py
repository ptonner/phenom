import numpy as np
import pystan
import os
import pickle
import pandas as pd


class Phenotype(object):

    def __init__(self, dataset, design, model='phenom.stan',
                 sigma_prior=[1, 1], maxExpectedCross=100, minExpectedCross=.01):

        self.dataset = dataset
        self.design = design

        self.x, self.y = self.dataset.standardize()
        self.dm, self.priors, self.alphaPriors, self.lengthscalePriors = self.design(self.dataset.meta)
        self.priors = self.priors + 1
        self.L = max(priors)

        self.sigma_prior = sigma_prior

        self.maxExpectedCross = maxExpectedCross
        self.minExpectedCross = minExpectedCross

        self.posteriors = []
        self._model = None
        self._modelFile = model

    @property
    def model(self):
        if self._model is None:
            self._model = pystan.StanModel(
                file='phenom/stan/' + self._modelFile)
        return self._model

    def config(self):

        cfg = {
            'N': self.x.shape[0],
            'P': self.y.shape[1],
            'K': self.dm.shape[1],
            'L': self.L,
            'prior': self.priors,
            'design': self.dm,
            'x': self.x,
            'y': self.y.T,
            'alpha_prior': self.alphaPriors,
            'lengthscale_prior': self.lengthscalePriors,
            'sigma_prior': self.sigma_prior,
        }

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
