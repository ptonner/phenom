import phenom

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("acid", choices=["benzoate", "citric", "malic"])
    parser.add_argument("model", choices=["mnull", "mbatch", "mfull"])
    parser.add_argument("--adapt_delta", type=float, default=0.8)
    parser.add_argument("--max_treedepth", type=int, default=10)
    parser.add_argument("--dataset", default=None)

    args = parser.parse_args()

    ds = phenom.dataset.DataSet.fromDirectory(
        "../data/pseudomonas/PA01-{}/".format(args.acid)
    )
    ds.meta["mMAcid"] = ds.meta["mM-acid"]
    ds.data = ds.data.iloc[5::3, :]
    ds.log()

    # designs
    formula = """1 + C(mMAcid, Sum, levels=[5, 10, 20, 0]) 
               + C(pH, Sum, levels=[6.5, 6.0, 5.5, 5, 7]) 
               + C(mMAcid, Sum, levels=[5, 10, 20, 0]):C(pH, Sum, levels=[6.5, 6.0, 5.5, 5, 7])"""

    mnull = phenom.design.Formula(ds.meta, formula)

    hierarchy = phenom.design.Formula(ds.meta, "1") + phenom.design.Formula(
        ds.meta, "C(plate) + 0"
    )
    mbatch = mnull * hierarchy

    # phenotype
    if args.model == "mnull":
        phen = phenom.phenotype.Phenotype(ds.data, mnull, model="phenom.stan")
    elif args.model == "mbatch":
        phen = phenom.phenotype.Phenotype(
            ds.data,
            mbatch,
            model="phenom_deriv.stan",
            lengthscale_priors=[[2.5, 0.45]] * 4 + [[100, 20]] * 4,
            alpha_priors=[[10, 10]] * 4 + [[7, 10]] * 4,
            minExpectedCross=0.01,
            maxExpectedCross=10,
            sigma_prior=[0.02, 20],
        )
    else:
        phen = phenom.phenotype.Phenotype(
            ds.data,
            mbatch,
            model="phenom_marginal.stan",
            lengthscale_priors=[[2.5, 0.45]] * 4 + [[100, 20]] * 4,
            alpha_priors=[[10, 10]] * 4 + [[7, 10]] * 4,
            minExpectedCross=0.01,
            maxExpectedCross=10,
            sigma_prior=[0.02, 20],
            marginal_lengthscale_prior=[4, 1],
            marginal_alpha_prior=[2, 20],
        )

    # sampling
    samples = phen.samples(
        control=dict(adapt_delta=args.adapt_delta, max_treedepth=args.max_treedepth)
    )

    # save
    if args.dataset is None:
        phen.save("paeruginosa/combined/{}/{}".format(args.acid, args.model))
    else:
        phen.save("paeruginosa/individual/{}/{}".format(args.acid, args.dataset))
