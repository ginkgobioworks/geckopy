Getting started
===============

Enzyme Constraint (gecko) models work as regular cobra models with:

-  Proteins as additional metabolites.
-  Protein pseudorreactions, which simulate the production of the
   protein by the cell (in :math:`mmol/g_{DW}/h`).
-  (Optional) Pool constraint: all proteins must be below the total
   quantified proteome.

The interface mirrors and depends on the
`cobrapy <https://github.com/opencobra/cobrapy/>`__ API.

.. code:: ipython3

    from os.path import join, pardir
    
    import geckopy
    
    ROOT = pardir
    DATA = join(ROOT, "tests", "data")

Import models
-------------

The models can be imported directly from SBML. Check the `design
chapter <design.html>`__ for more information on how proteins are
identified.

.. code:: ipython3

    ec_model = geckopy.io.read_sbml_ec_model(join(DATA, "eciML1515.xml.gz"))

:class:`~geckopy.Model` derives from
:doc:`cobrapy:autoapi/cobra/core/model/index`.

.. code:: ipython3

    print(f"Model has {len(ec_model.reactions)} reactions and {len(ec_model.metabolites)} metabolites.")


.. parsed-literal::

    Model has 4824 reactions and 2333 metabolites.


Models can also be written to the SBML file.

.. code:: ipython3

    geckopy.io.write_sbml_ec_model(ec_model, join(DATA, "model_copy.xml"))

Python API: Protein isolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Proteins are separated from metabolites and reactions in a different
container, :attr:`~geckopy.model.Model.proteins`.

.. code:: ipython3

    ec_model.slim_optimize()




.. parsed-literal::

    0.876997313396373



.. code:: ipython3

    ec_model.proteins.prot_P0A796




.. raw:: html

    
    <table>
        <tr>
            <td><strong>Protein identifier</strong></td><td>prot_P0A796</td>
        </tr><tr>
            <td><strong>Name</strong></td><td>prot_P0A796 [cytosol]
            </td>
        </tr><tr>
            <td><strong>Memory address</strong></td>
            <td>0x0%x140634108785232</td>
        </tr><tr>
        </tr><tr>
            <td><strong>Concentration</strong></td><td>nan</td>
        </tr><tr>
            <td><strong>Upper bound</strong></td><td>1000.0</td>
        </tr><tr>
        </tr><tr>
            <td><strong>Mw</strong></td><td>0.0</td>
        </tr><tr>
            <td><strong>In 3 reaction(s)</strong></td><td>
                PFKNo1 (62.00), PFK_2No1 (62.00), PFK_3No1 (62.00)
            </td>
        </tr>
    </table>




Reactions are aware of their proteins.

.. code:: ipython3

    ec_model.reactions.PFKNo1.proteins




.. parsed-literal::

    {<Protein prot_P0A796 at 0x7fe7ee18f250>: -4.4803e-06}



Proteins are aware of their reactions.

.. code:: ipython3

    ec_model.proteins.prot_P0A796.reactions




.. parsed-literal::

    frozenset({<Reaction PFKNo1 at 0x7fe7ecc97b80>,
               <Reaction PFK_2No1 at 0x7fe7ed25ce20>,
               <Reaction PFK_3No1 at 0x7fe7ed413e50>})



Analogous to
:ref:`cobrapy:autoapi/cobra/index.html#cobra.Reaction.flux`.

.. code:: ipython3

    ec_model.proteins.prot_P0A796.contribution




.. parsed-literal::

    2.7670111370232322e-05



Fluxes are separated in the solution dataframe to avoid regexing:

.. code:: ipython3

    # Fluxes are separated in the solution dataframe to avoid regexing:
    solution_rxn, solution_prot = ec_model.optimize()

.. code:: ipython3

    solution_rxn




.. raw:: html

    <strong><em>Optimal</em> solution with objective value 0.877</strong><br><div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>fluxes</th>
          <th>reduced_costs</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>EX_acgam_e</th>
          <td>0.000000</td>
          <td>-2.821258e-01</td>
        </tr>
        <tr>
          <th>EX_cellb_e</th>
          <td>0.000000</td>
          <td>-3.630018e-01</td>
        </tr>
        <tr>
          <th>EX_chol_e</th>
          <td>0.000000</td>
          <td>-2.633174e-02</td>
        </tr>
        <tr>
          <th>EX_pi_e</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>EX_h_e</th>
          <td>8.058201</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>PUACGAMS_REVNo1</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>ARHGDx_REVNo1</th>
          <td>0.000000</td>
          <td>5.312591e-17</td>
        </tr>
        <tr>
          <th>UDPGPT_REVNo1</th>
          <td>0.000000</td>
          <td>-5.551115e-17</td>
        </tr>
        <tr>
          <th>4HTHRA_REVNo1</th>
          <td>0.000587</td>
          <td>0.000000e+00</td>
        </tr>
        <tr>
          <th>RHMND_REVNo1</th>
          <td>0.000000</td>
          <td>0.000000e+00</td>
        </tr>
      </tbody>
    </table>
    <p>4824 rows × 2 columns</p>
    </div>



.. code:: ipython3

    solution_prot




.. raw:: html

    <strong><em>Optimal</em> solution with objective value 0.877</strong><br><div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>fluxes</th>
          <th>reduced_costs</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>prot_P0A825</th>
          <td>4.246198e-07</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>prot_P75823</th>
          <td>0.000000e+00</td>
          <td>-0.0</td>
        </tr>
        <tr>
          <th>prot_P0AEA8</th>
          <td>0.000000e+00</td>
          <td>-0.0</td>
        </tr>
        <tr>
          <th>prot_P36553</th>
          <td>3.195620e-06</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>prot_P06715</th>
          <td>1.642762e-07</td>
          <td>0.0</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>prot_P77215</th>
          <td>0.000000e+00</td>
          <td>-0.0</td>
        </tr>
        <tr>
          <th>prot_P0A8Y8</th>
          <td>0.000000e+00</td>
          <td>-0.0</td>
        </tr>
        <tr>
          <th>prot_P76290</th>
          <td>0.000000e+00</td>
          <td>-0.0</td>
        </tr>
        <tr>
          <th>prot_P16691</th>
          <td>0.000000e+00</td>
          <td>-0.0</td>
        </tr>
        <tr>
          <th>prot_P32138</th>
          <td>0.000000e+00</td>
          <td>-0.0</td>
        </tr>
      </tbody>
    </table>
    <p>1259 rows × 2 columns</p>
    </div>



Kcats
-----

The kcats can be inspected and manipulated from the
:class:`~geckopy.protein.Protein` object, as a regular dictionary.

These kcats are individual for every protein-reaction pair and
correspond to the stoichiometric coefficient of the protein
pseudometabolite in the reaction.

-  The units of the input are in :math:`\frac{1}{s}`.
-  This input is translated to :math:`h` in the stoichiometric
   coefficient.

.. code:: ipython3

    ec_model.proteins.prot_P0A796.kcats




.. parsed-literal::

    {<Reaction PFKNo1 at 0x7fe7ecc97b80>: 61.99981648054322, <Reaction PFK_2No1 at 0x7fe7ed25ce20>: 61.99981648054322, <Reaction PFK_3No1 at 0x7fe7ed413e50>: 61.99981648054322}



.. code:: ipython3

    ec_model.proteins.prot_P0A796.kcats["PFKNo1"]




.. parsed-literal::

    61.99981648054322



.. code:: ipython3

    ec_model.reactions.PFKNo1.metabolites[ec_model.proteins.prot_P0A796]




.. parsed-literal::

    -4.4803e-06



.. code:: ipython3

    ec_model.proteins.prot_P0A796.kcats["PFKNo1"] = 120

.. code:: ipython3

    ec_model.reactions.PFKNo1.metabolites[ec_model.proteins.prot_P0A796]




.. parsed-literal::

    -2.3148148148148148e-06


