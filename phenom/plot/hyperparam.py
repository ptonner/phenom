import scipy
import matplotlib.pyplot as plt
import numpy as np


def kde(sample, lim=1.0, norm=True, ** kwargs):
    l, h = sample.min(), sample.max()
    m = (l + h) / 2
    r = h - l
    l, h = m - lim * r / 2, m + lim * r / 2

    z = np.linspace(l, h)
    kde = scipy.stats.gaussian_kde(sample)
    p = kde(z)
    if norm:
        p = p / p.max()
    plt.plot(z, p, **kwargs)
