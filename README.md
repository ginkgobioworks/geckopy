# geckopy 
[![CI-CD](https://github.com/ginkgobioworks/geckopy/actions/workflows/main.yml/badge.svg)](https://github.com/ginkgobioworks/geckopy/actions/workflows/main.yml)

**G**enome-scale model **E**nzyme **C**onstraints, using **K**inetics and 
**O**mics in **py**thon.

By combining kcats and proteomics measurement, geckopy allows for improving
the modeling capabilities in genome-scale models.

> Based on [Sanchez et al., 2017](https://dx.doi.org/10.15252/msb.20167411).

Check [https://github.com/SysBioChalmers/GECKO](https://github.com/SysBioChalmers/GECKO)
for the matlab counterpart.


## Getting started
Load a model.
```python
import geckopy

model = geckopy.io.read_sbml_ec_model("tests/data/eciML1515.xml.gz")
model.optimize()
```

Add copy number experimental data.
```python
import pandas as pd
from geckopy.experiment import from_copy_number

raw_proteomics = pd.read_csv("tests/data/ecoli_proteomics_schmidt2016S5.tsv")
exp_model = from_copy_number(
    model,
    index=raw_proteomics["uniprot"],
    cell_copies=raw_proteomics["copies_per_cell"],
    stdev=raw_proteomics["stdev"],
    vol=2.3,
    dens=1.105e-12,
    water=0.3,
)
exp_model.optimize()
```

Add pool constraint.
```python
# add some molecular weights to the proteins if the model does not have them
for prot in ec_model.proteins:
    prot.mw = 330
exp_model.constrain_pool(
    p_measured=12.,
    sigma_saturation_factor=0.8,
    fn_mass_fraction_unmeasured_matched=0.2,
)
print(exp_model.optimize())
print(exp_model.protein_pool_exchange)
```
