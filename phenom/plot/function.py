import matplotlib.pyplot as plt


def interval(x, f, **kwargs):
    plt.plot(x, (f).mean(0), **kwargs)
    plt.fill_between(x,
                     (f).mean(0) - 2 * (f).std(0),
                     (f).mean(0) + 2 * (f).std(0), alpha=.4)
