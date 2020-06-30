import phenom
import os
import scipy
import pandas as pd
import patsy
import numpy as np


def load_data(condition, DATA_DIR="../data"):

    if condition == "standard":
        DATA_DIR = os.path.join(DATA_DIR, "standard")
    else:
        DATA_DIR = os.path.join(DATA_DIR, "{}-oxidative".format(condition))

    # load data
    ds = None

    for dr in os.listdir(DATA_DIR):

        if not os.path.isdir(os.path.join(DATA_DIR, dr)):
            continue

        if ds is None:
            ds = phenom.dataset.DataSet.fromDirectory(os.path.join(DATA_DIR, dr))
        else:
            ds = ds.concat(
                phenom.dataset.DataSet.fromDirectory(os.path.join(DATA_DIR, dr))
            )

    ds.filter()
    ds.meta["mMPQ"] = ds.meta["mM PQ"]
    ds.data = np.log2(ds.data)
    ds.data = ds.data.iloc[5::3, :]
    return ds


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("condition", choices=["standard", "low", "hi"])
    parser.add_argument("model", choices=["mnull", "mbatch", "mfull"])
    parser.add_argument("--adapt_delta", type=float, default=0.8)
    parser.add_argument("--max_treedepth", type=int, default=10)

    args = parser.parse_args()

    ds = load_data(args.condition)

    # designs
    mnull = phenom.design.Formula(ds.meta, "C(mMPQ, Sum)")

    hierarchy = phenom.design.Formula(ds.meta, "1") + phenom.design.Formula(
        ds.meta, "C(plate) + 0"
    )
    mbatch = mnull * hierarchy

    # phenotype
    if args.model == "mnull":
        phen = phenom.phenotype.Phenotype(ds.data, mnull, model="phenom.stan")
    elif args.model == "mbatch":
        phen = phenom.phenotype.Phenotype(ds.data, mbatch, model="phenom_deriv.stan")
    else:
        phen = phenom.phenotype.Phenotype(ds.data, mbatch, model="phenom_marginal.stan")

    # sampling
    samples = phen.samples(
        control=dict(adapt_delta=args.adapt_delta, max_treedepth=args.max_treedepth)
    )
    phen.save("hsalinarum/{}/{}".format(args.condition, args.model))