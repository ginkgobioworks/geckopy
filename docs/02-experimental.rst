Experimental data
=================

-  The Final proteomics data in the model has units of
   :math:`\frac{\text{mmol}}{\text{gDW}}`.
-  geckopy provides the transformation from copies per cell.

.. math::

  \begin{align}
    \frac{\text{mmol}}{\text{cell}} &= \frac{\bf{molecules}}{\bf{cell}} \frac{10^3\text{mmol}}{\text{6.022 10^23 molecules}} \\
    \frac{\text{mmol}}{\text{gDW}} &= \frac{\bf{mmol}}{\bf{cell}} \frac{\text{cell}}{fL} \frac{fL}{g}\frac{g}{\text{gDW}}
  \end{align}

.. code:: ipython3

    from os.path import join, pardir
    
    import geckopy
    import pandas as pd
    
    ROOT = pardir
    DATA = join(ROOT, "tests", "data")
    
    ec_model = geckopy.io.read_sbml_ec_model(join(DATA, "eciML1515.xml.gz"))

.. code:: ipython3

    raw_proteomics = pd.read_csv(join(DATA, "ecoli_proteomics_schmidt2016S5.tsv"))
    
    ec_model_exp = geckopy.experimental.from_copy_number(
        ec_model,
        # the index should be the IDs of the proteins exactly as in the model!
        index=raw_proteomics["uniprot"].apply(lambda x: f"prot_{x}"),
        cell_copies = raw_proteomics["copies_per_cell"], 
        stdev = raw_proteomics["stdev"],
        vol=2.3, dens=1.105e-12, water=0.3
    )

.. code:: ipython3

    ec_model_exp.slim_optimize()




.. parsed-literal::

    nan



After the proteomics data is applied, it is often the case that the
model requires some relaxation of the experimental constraints to be
able to grow.

Relaxation
~~~~~~~~~~

Geckopy provides elastic filtering (see
:func:`~geckopy.experimental.relaxation.apply_proteomics_elastic_relaxation`)
as implemented in `Chinnek and Dravnieks,
1990 <https://pubsonline.informs.org/doi/abs/10.1287/ijoc.3.2.157>`__).

.. code:: ipython3

    from geckopy.experimental.relaxation import (
        apply_proteomics_elastic_relaxation,
        apply_proteomics_relaxation,
    )

This method returns the relaxed model with the Irreducibly Inconsistent
Set of functional constraints (IIS); that is, all the proteins that
affect the feasibility of the problem.

Please note that not all proteins variables in the ISS must be relaxed
for the model to be feasible, but all of the proteins in the ISS are
part of one set of proteins that makes the model feasible.

.. code:: ipython3

    relaxed_model, iiset = apply_proteomics_elastic_relaxation(ec_model_exp)

.. code:: ipython3

    relaxed_model.slim_optimize()




.. parsed-literal::

    0.8588931565514887



Alternatively, the relaxation can be applied to just the first found
subset of the ISS with
:func:`~geckopy.experimental.relaxation.apply_proteomics_relaxation`:

.. code:: ipython3

    relaxed_model, iset = apply_proteomics_relaxation(ec_model_exp)

.. code:: ipython3

    relaxed_model.slim_optimize()




.. parsed-literal::

    0.8588940541385824



Pool constraint
~~~~~~~~~~~~~~~

.. code:: ipython3

    ec_model = geckopy.io.read_sbml_ec_model(join(DATA, "eciML1515.xml.gz"))

A pool constraint can be applied to the :class:`~geckopy.protein.Proteins`\ s
to account for protein crowding. This is useful when there are proteins with
missing concentrations in the model but the total amount of protein that the
cell can allocate is known.

The amount of flux a protein can take from the pool is their :math:`M_w`
(in :math:`\frac{g}{mmol}`). This value can be scrapped with
:func:`~geckopy.experimental.molecular_weights.extract_proteins`

.. code:: ipython3

    from geckopy.experimental.molecular_weights import extract_proteins

.. code:: ipython3

    df = extract_proteins(ec_model)
    for row in df.itertuples():
        ec_model.proteins.get_by_id(row[2]).mw = row[3]

As explained in the `Appendix of Sánchez et al.,
2017 <https://www.embopress.org/action/downloadSupplement?doi=10.15252%2Fmsb.20167411&file=msb167411-sup-0001-Appendix.pdf>`__,

-  ``p_total`` is the protein mass in :math:`\frac{g}{g_{DW}}` of the
   proteins in the model.
-  ``sigma_saturation_factor`` is the parameter adjusting how much of a
   protein pool can take part in reactions; i.e, how much of the
   measured proteome it is participation in the metabolism.
-  ``fn_mass_fraction_unmeasured_matched`` is
   :math:`\frac{f_n}{1 - f_m}`, where :math:`f_n` is the mass fraction
   of unmeasured protein divided and :math:`f_m` is the fraction of
   proteins measured. This way, it is 1 if no protein concentration is
   known.

.. code:: ipython3

    ec_model.constrain_pool(
        p_total=0.2,
        sigma_saturation_factor=0.8, 
        fn_mass_fraction_unmeasured_matched=1
    )
    ec_model.protein_pool_exchange




.. raw:: html

    
    <table>
        <tr>
            <td><strong>Reaction identifier</strong></td><td>prot_pool_exchange</td>
        </tr><tr>
            <td><strong>Name</strong></td><td></td>
        </tr><tr>
            <td><strong>Memory address</strong></td>
            <td>0x07fe6339ac760</td>
        </tr><tr>
            <td><strong>Stoichiometry</strong></td>
            <td>
                <p style='text-align:right'>--> prot_pool</p>
                <p style='text-align:right'>--></p>
            </td>
        </tr><tr>
            <td><strong>GPR</strong></td><td></td>
        </tr><tr>
            <td><strong>Lower bound</strong></td><td>0</td>
        </tr><tr>
            <td><strong>Upper bound</strong></td><td>0.16000000000000003</td>
        </tr>
    </table>




.. code:: ipython3

    ec_model.slim_optimize()




.. parsed-literal::

    0.26126914095190934


