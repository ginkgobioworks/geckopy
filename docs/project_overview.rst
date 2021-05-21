Project overview
================

.. code:: shell

    ┌── __init__.py
    ├── model.py
    ├── protein.py
    ├── reaction.py
    ├── io.py
    ├── flux_analysis.py
    │ ┌── __init__.py
    │ ├── relaxation.py
    │ ├── pytfa.py
    ├─┴ integration
    │ ┌── __init__.py
    │ ├── experimental.py                   
    │ ├── molecular_weights.py
    │ ├── relaxation.py
    ├─┴ experimental
  ┌─┴ geckopy


geckopy/model.py
~~~~~~~~~~~~~~~~

The :class:`~geckopy.Model` class, which is derived from 
:doc:`cobrapy:autoapi/cobra/core/model/index`.

* All methods that directly mention :class:`~geckopy.Reaction` were overwritten.
  This would be also required if any component of cobrapy's `Model` is derived
  in geckopy.
* Proteins were added to the population of the model.
* Pool constraint methods.

geckopy/protein.py
~~~~~~~~~~~~~~~~~~

The :class:`~geckopy.Protein` class is both represented as a constraint (the metabolite)
and as a variable (the exchange pseudorreaction). Thus, it contains methods both 
from cobrapy's `Metabolite` and `Reaction`.

geckopy/reaction.py
~~~~~~~~~~~~~~~~~~~

The class :class:`~geckopy.Reaction` merely extends cobrapy's `Reaction` to
point to its proteins.

geckopy/io.py
~~~~~~~~~~~~~

Serialization and deserialization to SBML (adapting `cobrapy:autoapi/cobra/io/sbml/index`).

geckopy/flux_analysis.py
~~~~~~~~~~~~~~~~~~~~~~~~

Cobrapy flux analysis functionality can generally be used for geckopy models
but some care must be taken to make fair comparisons of flux variability analysis,
as pointed by `Domenzain et al., 2021`. This submodule includes that and some
utility functions to inspect and process the protein solution.

geckopy/integration/
~~~~~~~~~~~~~~~~~~~~

Integration layer for `pytfa <https://github.com/EPFL-LCSB/pytfa/>`__. It
includes relaxation algorithms for the integration of metabolomics with
proteins.

geckopy/experimental/
~~~~~~~~~~~~~~~~~~~~~

Functionality for integration of experimental data. It includes the relaxation
algorithms for proteomics experimental data and utility functions for the
annotation of molecular weights from `Uniprot accesion numbers`_.

.. _Uniprot accesion numbers: https://www.uniprot.org/help/accession_numbers
.. _Domenzain et al., 2021: https://www.biorxiv.org/content/10.1101/2021.03.05.433259v1
