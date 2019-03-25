import matplotlib.pyplot as plt
import scipy.stats


def interval(x, f, color="C0", alpha=.4, conf=.95, **kwargs):
    plt.plot(x, (f).mean(0), color=color, **kwargs)

    # bounds from standard normal
    conf = 1 - conf
    lw = scipy.stats.norm.ppf(conf/2)
    hi = scipy.stats.norm.ppf(1 - conf/2)

    plt.fill_between(x,
                     (f).mean(0) + lw * (f).std(0),
                     (f).mean(0) + hi * (f).std(0), alpha=alpha, color=color)
