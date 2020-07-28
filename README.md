
*phenom*: A hierarchical non-parametric microbial phenotype model

# Requirements

## OS requirements
This software has been run on *macOS* and *linux*, specifically:
* ubuntu 16.04
* osX 10.14.1
  * note, special steps may be needed to install *pystan* on osX
    operating systems. [See the pystan wiki for more details.](https://github.com/stan-dev/pystan/wiki/PyStan-and-OS-X)
	
## Python

* python (2.7 or 3.6)
* major dependencies:
  * pystan
  * numpy
  * patsy
  * matplotlib
  * GPy (for running examples)

# Installation

Download the phenom repository. It is recommended to use a python
virtual environment for installation. Run the install (this should
take under five minutes):

	python setup.py install
	
# Usage

## Example

An example notebook is provided in the `notebooks` folder. 

## Phenotype

The core interface to *phenom* models is through the `phenotype`
object. `phenotype`'s requires two major components:
1. A dataset
2. and a design
The following sections outline these components in detail.

## Dataset

A dataset consists of raw growth data and associated metadata for each
well in the dataset. *phenom* expects both of these to be provided as
pandas dataframes when creating a new dataset:
```python
# load data and meta as pandas dataframes here

from phenom.dataset import DataSet
ds = DataSet(data, meta)
```

The shapes between raw growth data and metadata must match. A data
file with have `NxK` data points of `N` timepoints and `K` individual
growth curves. The metadata will then be shape `KxM` with `M` metadata elements.

A useful data storage pattern is to save the data and metadata as csv
files in the same directory:
```
path/
  to/
    folder/
      data.csv
      meta.csv
```
*phenom* then provides a convience function to load these files
together for a dataset (note that data.csv and meta.csv are the
assumed filenames):
```python
ds = DataSet.fromDirectory("path/to/folder")
```

## Building designs

*phenom* models require a design specifying the relationship between
metadata and the latent functions to be estimated. Currently, design
construction is supported through the use of
[patsy](https://patsy.readthedocs.io) formulas to convert metadata
into a design matrix.

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

### composing more complicated designs

*patsy* supports the use of compositional operations for combining designs. These composition operations are:
* addition (`d1 + d2`): corresponds to concatenating the columns of two designs
* multiplication (`d1 * d2`): corresponds to the kronecker product of design columns. useful for repeating a design at multiple hierarchical levels (see below)

### modeling batch effects

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

# License

This project is covered under the **Apache 2.0 License**
