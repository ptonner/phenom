
*phenom*: A hierarchical non-parametric microbial phenotype model

# building designs

*phenom* models require a design specifying the relationship between metadata and the latent functions to be estimated. Currently, design construction is supported through the use of [patsy](https://patsy.readthedocs.io) formulas to convert metadata into a design matrix.

For example, to create a design for data with the following metadata:

| strain  | condition |
| ------------- | ------------- |
| parent  | standard  |
| mutant  | standard |
| parent  | stress |
| mutant  | stress |

A design can be constructed as:

```python
from phenom.design import Formula

treatment = Formula(meta, 'C(strain) + C(condition) + C(strain):C(condition)')
```
`C(strain)` and `C(condition)` specify categorical variables, and `C(strain):C(condition)` specifies an interaction between strain and condition effects. For more details on equation formatting see [the patsy docs](https://patsy.readthedocs.io/en/latest/formulas.html#the-formula-language).

The output of this design (`treatment.frame`) is:

| mean  | strain=mutant  | condition=stress |
| ------------- | ------------- | ------------- |
1 | 0  | 0  |
1 | 1  | 0 |
1 | 0  | 1 |
1 | 1  | 1 |

## composing more complicated designs

*patsy* supports the use of compositional operations for combining designs. These composition operations are:
* addition (`d1 + d2`): corresponds to concatenating the columns of two designs
* multiplication (`d1 * d2`): corresponds to the kronecker product of design columns. useful for repeating a design at multiple hierarchical levels (see below)

## modeling batch effects

To model batch effects, consider metadata of the form

| strain  | condition | batch |
| ------------- | ------------- | ------------- |
| parent  | standard  | 1
| mutant  | standard | 1
| parent  | stress | 1
| mutant  | stress | 1
| parent  | standard  | 2
| mutant  | standard | 2
| parent  | stress | 2
| mutant  | stress | 2

We combine the treatment design described above with a design corresponding to batch effects to make our complete design:
```python

# this is the design to be modeled both at the global and batch level
treatment = Formula(meta, 'C(strain) + C(condition)')

# base phenotype common to all observations, design is a column of 1's
base = Formula(meta, '1')

# batch effects
# the '+0' in the formula is necessary to prevent patsy from creating an un-desired intercept column
batch = Formula(meta, 'C(batch) + 0')

# hierarchy is a combination of global phenotype and batch effects
hierarchy = base + batch

# the full design replicates treatment design across hierarchy
design = treatment * hierarchy
```
