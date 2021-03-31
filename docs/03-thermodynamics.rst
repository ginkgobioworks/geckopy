Thermodynamics integration
==========================

Geckopy provides an integration layer with
`pytfa <https://github.com/EPFL-LCSB/pytfa/>`__ to run Thermodynamic
Flux Analysis with an enzyme constraint model.

Adapted from `pytfa
tutorials <https://github.com/EPFL-LCSB/pytfa/blob/master/tutorials/figure_paper.py>`__.

.. code:: ipython3

    import geckopy
    import pandas as pd
    import pytfa
    
    from cobra.flux_analysis.variability import flux_variability_analysis
    from geckopy.integration.pytfa import (
        adapt_gecko_to_thermo,
        translate_model_mnx_to_seed,
    )
    from geckopy.experimental import from_copy_number
    from geckopy.experimental.relaxation import apply_proteomics_relaxation
    
    from pytfa.io import load_thermoDB
    from pytfa.analysis import variability_analysis

.. code:: ipython3

    CPLEX = 'optlang-cplex'
    GUROBI = 'optlang-gurobi'
    GLPK = 'optlang-glpk'
    ROOT = pardir
    DATA = join(ROOT, "tests", "data")

Load the model
~~~~~~~~~~~~~~

Load the enzyme constraint model and transform it to a thermo model.

In this step, the proteins are addecuately added to the Thermodynamic
model with a formation energy of 0. Thus, proteins pseudometabolites are
ignored by the :math:`\Delta G_r` calculations. This is required so that they
are not interpreted as metabolites with missing thermodynamic information,
which would invalidate all reactions with proteins.

.. code:: ipython3

    ec_model = geckopy.io.read_sbml_ec_model(join(DATA, "eciML1515.xml.gz"))

.. code:: ipython3

    # Transfrom it into thermo model
    thermodb = load_thermoDB(join(DATA, "thermo_data.thermodb"))
    compartment_data = pytfa.io.read_compartment_data(
        join(DATA, "compartment_data.json")
    )
    translate_model_mnx_to_seed(
        ec_model, thermodb, join(DATA, "chem_xref_seedset.tsv")
    )
    tmodel = adapt_gecko_to_thermo(
        ec_model, thermodb, compartment_data, solver=CPLEX
    )
    
    # Info on the cobra_model
    tmodel.print_info()


.. parsed-literal::

    2021-03-19 11:32:09,278 - thermomodel_ecModel of iML1515 - INFO - # Model initialized with units kcal/mol and temperature 298.15 K
    2021-03-19 11:32:11,752 - thermomodel_ecModel of iML1515 - INFO - # Model preparation starting...
    2021-03-19 11:32:18,493 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:18,494 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:18,792 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:18,792 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:18,914 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 11:32:18,914 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 11:32:18,934 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 11:32:18,936 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 11:32:18,939 - thermomodel_ecModel of iML1515 - WARNING - Warning : Cs/Cs
    2021-03-19 11:32:18,940 - thermomodel_ecModel of iML1515 - WARNING - Warning : Cs/Cs
    2021-03-19 11:32:18,986 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:18,986 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:19,191 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:19,192 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:19,243 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 11:32:19,244 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 11:32:19,246 - thermomodel_ecModel of iML1515 - WARNING - Warning : Cs/Cs
    2021-03-19 11:32:19,247 - thermomodel_ecModel of iML1515 - WARNING - Warning : Cs/Cs
    2021-03-19 11:32:19,334 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:19,706 - thermomodel_ecModel of iML1515 - WARNING - Warning : O3Te/Te
    2021-03-19 11:32:19,707 - thermomodel_ecModel of iML1515 - WARNING - Warning : CH3O3Te/Te
    2021-03-19 11:32:19,725 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:19,752 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 11:32:19,802 - thermomodel_ecModel of iML1515 - INFO - # Model preparation done.
    2021-03-19 11:32:19,804 - thermomodel_ecModel of iML1515 - INFO - # Model conversion starting...
    2021-03-19 11:32:50,122 - thermomodel_ecModel of iML1515 - INFO - # Model conversion done.
    2021-03-19 11:32:50,123 - thermomodel_ecModel of iML1515 - INFO - # Updating cobra_model variables...
    2021-03-19 11:32:50,219 - thermomodel_ecModel of iML1515 - INFO - # cobra_model variables are up-to-date


.. parsed-literal::

                                  value
    key                                
    name             ecModel of iML1515
    description      ecModel of iML1515
    num constraints               19819
    num variables                 25248
    num metabolites                2333
    num reactions                  4824
                               value
    key                             
    num metabolites(thermo)     2264
    num reactions(thermo)        585
    pct metabolites(thermo)  97.0424
    pct reactions(thermo)    12.1269


Check that the model works.

.. code:: ipython3

    solution = tmodel.optimize()
    solution




.. raw:: html

    <strong><em>Optimal</em> solution with objective value 0.885</strong><br><div>
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
          <td>None</td>
        </tr>
        <tr>
          <th>EX_cellb_e</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>EX_chol_e</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>EX_pi_e</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>EX_h_e</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>PUACGAMS_REVNo1</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>ARHGDx_REVNo1</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>UDPGPT_REVNo1</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>4HTHRA_REVNo1</th>
          <td>0.000592</td>
          <td>None</td>
        </tr>
        <tr>
          <th>RHMND_REVNo1</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
      </tbody>
    </table>
    <p>4824 rows × 2 columns</p>
    </div>



Flux variability analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

We can use pytfa’s `variability_analysis` for the integrated model. We
will run four simulations:

1. FBA-based (no enzyme nor thermodynamic constraints).
2. TFA-based.
3. TFA + metabolite concentrations
4. TFA + metabolites concentrations + enzyme constraints.

.. code:: ipython3

    fva_fluxes = flux_variability_analysis(ec_model)
    tva_fluxes = variability_analysis(tmodel, kind='reactions')


.. parsed-literal::

    2021-03-19 11:33:02,042 - thermomodel_ecModel of iML1515 - INFO - Beginning variability analysis for variable of type reactions
    minimizing: 100%|██████████| 4824/4824 [40:45<00:00,  1.97it/s] 
    maximizing: 100%|██████████| 4824/4824 [1:30:22<00:00,  1.12s/it]


Save the results just in case.

.. code:: ipython3

    fva_fluxes.to_csv("fva_fluxes.tsv", sep="\t", index_label="reaction")
    tva_fluxes.to_csv("tva_fluxes.tsv", sep="\t", index_label="reaction")

Now, the same with specific concentration data.

.. code:: ipython3

    # Add more specific concentration data
    def apply_concentration_bound(met, lb, ub):
        the_conc_var = tmodel.log_concentration.get_by_id(met)
        # Do not forget the variables in the model are logs !
        the_conc_var.ub = log(ub)
        the_conc_var.lb = log(lb)

.. code:: ipython3

    apply_concentration_bound('atp_c', lb=1e-3, ub=1e-2)
    apply_concentration_bound('adp_c', lb=4e-4, ub=7e-4)
    apply_concentration_bound('amp_c', lb=2e-4, ub=3e-4)
    
    tmodel.optimize()
    # Perform variability analysis again
    tva_fluxes_lc = variability_analysis(tmodel, kind='reactions')


.. parsed-literal::

    2021-03-19 13:44:11,972 - thermomodel_ecModel of iML1515 - INFO - Beginning variability analysis for variable of type reactions
    minimizing: 100%|██████████| 4824/4824 [54:08<00:00,  1.48it/s]  
    maximizing: 100%|██████████| 4824/4824 [1:33:59<00:00,  1.17s/it]


.. code:: ipython3

    tva_fluxes_lc.to_csv("tva_fluxes_lc.tsv", sep="\t", index_label="reaction")

And finally, thermo + concentration + proteomics.

.. code:: ipython3

    raw_proteomics = pd.read_csv(join(DATA, "ecoli_proteomics_schmidt2016S5.tsv"))
    ec_model_constrained = from_copy_number(
        ec_model,
        index=raw_proteomics["uniprot"],
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    # we need to relax the experimental assumptions!
    relaxed_model, iss = apply_proteomics_relaxation(ec_model_constrained)
    tmodel_prot = adapt_gecko_to_thermo(
        relaxed_model, thermodb, compartment_data, solver=CPLEX
    )


.. parsed-literal::

    2021-03-19 16:17:36,090 - thermomodel_ecModel of iML1515 - INFO - # Model initialized with units kcal/mol and temperature 298.15 K
    2021-03-19 16:17:38,553 - thermomodel_ecModel of iML1515 - INFO - # Model preparation starting...
    2021-03-19 16:17:45,492 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:45,493 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:45,772 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:45,773 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:45,888 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 16:17:45,889 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 16:17:45,906 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 16:17:45,908 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 16:17:45,911 - thermomodel_ecModel of iML1515 - WARNING - Warning : Cs/Cs
    2021-03-19 16:17:45,912 - thermomodel_ecModel of iML1515 - WARNING - Warning : Cs/Cs
    2021-03-19 16:17:45,953 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:45,954 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:46,156 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:46,157 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:46,204 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 16:17:46,205 - thermomodel_ecModel of iML1515 - WARNING - Warning : F/F
    2021-03-19 16:17:46,208 - thermomodel_ecModel of iML1515 - WARNING - Warning : Cs/Cs
    2021-03-19 16:17:46,209 - thermomodel_ecModel of iML1515 - WARNING - Warning : Cs/Cs
    2021-03-19 16:17:46,294 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:46,667 - thermomodel_ecModel of iML1515 - WARNING - Warning : O3Te/Te
    2021-03-19 16:17:46,668 - thermomodel_ecModel of iML1515 - WARNING - Warning : CH3O3Te/Te
    2021-03-19 16:17:46,685 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:46,711 - thermomodel_ecModel of iML1515 - WARNING - Warning : C16H12NO3ClF/F
    2021-03-19 16:17:46,761 - thermomodel_ecModel of iML1515 - INFO - # Model preparation done.
    2021-03-19 16:17:46,762 - thermomodel_ecModel of iML1515 - INFO - # Model conversion starting...
    2021-03-19 16:18:16,250 - thermomodel_ecModel of iML1515 - INFO - # Model conversion done.
    2021-03-19 16:18:16,251 - thermomodel_ecModel of iML1515 - INFO - # Updating cobra_model variables...
    2021-03-19 16:18:16,341 - thermomodel_ecModel of iML1515 - INFO - # cobra_model variables are up-to-date


.. code:: ipython3

    tva_fluxes_prot = variability_analysis(tmodel_prot, kind='reactions')


.. parsed-literal::

    2021-03-19 16:22:07,789 - thermomodel_ecModel of iML1515 - INFO - Beginning variability analysis for variable of type reactions
    minimizing: 100%|██████████| 4824/4824 [19:53<00:00,  4.04it/s]
    maximizing: 100%|██████████| 4824/4824 [1:56:38<00:00,  1.45s/it]  


.. code:: ipython3

    tva_fluxes_prot.to_csv("tva_fluxes_prot.tsv", sep="\t", index_label="reaction")

.. code:: ipython3

    tva_fluxes_prot.head()




.. raw:: html

    <div>
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
          <th>minimum</th>
          <th>maximum</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>EX_acgam_e</th>
          <td>0.0</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>EX_cellb_e</th>
          <td>0.0</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>EX_chol_e</th>
          <td>0.0</td>
          <td>0.000000</td>
        </tr>
        <tr>
          <th>EX_pi_e</th>
          <td>0.0</td>
          <td>1000.000000</td>
        </tr>
        <tr>
          <th>EX_h_e</th>
          <td>0.0</td>
          <td>165.250256</td>
        </tr>
      </tbody>
    </table>
    </div>


