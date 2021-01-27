---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.9.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Proof of concept

```python
from os.path import join, pardir
import pandas as pd

from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis

ROOT = pardir
DATA = join(ROOT, "data")
```

```python
# need a name
import ecgem
```

```python
model = read_sbml_model(join(DATA, "eciML1515.xml.gz"))
```

## 1. Load the model
- Final data from $\frac{\text{mmol}}{\text{gDW}}$.
- Usual dataset coming from fermentation: copies per cell?
- Package could provide the tranformation from copies per cell (or other relevant choices):

\begin{align}
  \frac{\text{mmol}}{\text{cell}} &= \frac{\bf{molecules}}{\bf{cell}} \frac{10^3\text{mol}}{\text{molecules}} \\
  \frac{\text{mmol}}{\text{gDW}} &= \frac{\bf{mmol}}{\bf{cell}} \frac{\text{cell}}{fL} \frac{fL}{g}\frac{g}{\text{gDW}}
\end{align}

```python
raw_proteomics = pd.read_csv(join(DATA, "copy_number.csv"))
# user needs to select relevant columns + the model
# it could additionally accept just the dataframe correctly named columns00
ec_model = ecgem.from_copy_number(
    model,
    index=raw_proteomics["uniprot"],
    cell_copies = raw_proteomics["copies_per_cell"], 
    stdev = raw_proteomics["stdev"],
    vol=2.3, dens=1.105e-12, water=0.3
)
```

```python
processed_proteomics = pd.read_csv(join(DATA, "mmol_gdW_protemics.csv"))
# alternatively, from already prepared data (mmol/gDW)
# the index is has the ids and the first column has the values
# it could potentially be a pd.Series or merely a dictionary
ec_model = ecgem.from_mmol_gDW(
    model,
    processed_proteomics
)
```

## 2. Write the model
- Can the proteomics data ($\frac{mmol}{gDW}$) be written to the SBML specification?

```python
ecgem.write_sbml_model(ec_model, "eciML1515_with_prots.xml.gz")
```

If that's the case, reading directly from the SBML should be possible:

```python
ec_model = ecgem.read_sbml_model("eciML1515_with_prots.xml.gz")
```

### SBML Proposition
Proteins are [Species](https://docs.rs/rust_sbml/0.5.2/rust_sbml/struct.Species.html) with:
* `substanceUnits`: A custom unit: prot mmol/gDW 
* `initalAmount` (optional): value (double)

Protein "Exchanges" can be one of:
* [Reactions](https://docs.rs/rust_sbml/0.5.2/rust_sbml/struct.Reaction.html), compatible
with current ec models.
* Constraints with message "protein":
```xml
<model>
    ...
    <listOfConstraints>
        <constraint>
            <math xmlns="http://www.w3.org/1998/Math/MathML"35xmlns:sbml="http://www.sbml.org/sbml/level3/version2/core">
                <apply>
                    <and/>
                        <apply> <lt/> <ci> S1 </ci> <cn sbml:units="mole"> 100 </cn>
                        </apply>
                </apply>
            </math>
            <message><p xmlns="http://www.w3.org/1999/xhtml">protein</p></message>
        </constraint>
    </listOfConstraints>
    ...
</model>
```
* Be ignored and generated on read from `initialAmount`.


## 3. Optimize the model
The model should expose the same API as a normal [cobrapy](https://cobrapy.readthedocs.io/) model.

```python
ec_model.slim_optimize()
```

```python
fva_solution = flux_variability_analysis(ec_model)
fva_solution
```

- Should the `solution` (pd.Dataframe) contain the protein fluxes?
- At least the fva_solution should, since the values cannot be read from a `flux` field.


### API Proposition
* `ecgem.Model` extends `cobra.Model`.
* `ecgem.Reaction`?
* `Reaction` has an attribute `.protein` of class `ecgem.Protein`.
    * Consider `ecgem.Reaction`?
    * Consider `.proteins` instead of `.protein`; i.e., __are there EC models that assign more than one protein to 
    one reaction?__
* `ecgem.Protein` has attributes:
    * `.id`: `str`
    * `.flux`: `float`, as in `cobra.Reaction.flux`.
    * `.concentration`: `float`.
    * `.reactions`: `frozenset[cobra.Reaction]`, emulating cobrapy metabolites <-> reactions
    behaviour.
* Solution flux as field of solution.

```python
solution = ec_model.optimize()
```

```python
solution
```

```python
ec_model.reactions.EX_acgam_e.protein.id
ec_model.reactions.EX_acgam_e.protein.flux
```

```python
ec_model.proteins.Q59385.flux
```

```python
prot_fluxes = {protein.id: protein.flux for protein in ec_model.proteins}
```

### 4. Flexibiization
- As a MILP problem.
- As the caffeine implementation (iterative removal of biggest contributor to shadow prices).

```python
ecgem.flexibilize(ec_model, experimental_growth=0.57)
```
