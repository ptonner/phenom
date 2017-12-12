
*phenom*: A hierarchical non-parametric microbial phenotype model

## design example

    from phenom.design import Formula

    treatment = Formula(meta, 'C(treatment1)')
    noise = Formula(meta, 'C(batch)')

    design = treatment * noise
