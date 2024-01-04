Thermodynamics integration
==========================

Geckopy provides an integration layer with
`pytfa <https://github.com/EPFL-LCSB/pytfa/>`__ to run Thermodynamic
Flux Analysis with an enzyme constraint model.

.. code:: ipython3

    from pathlib import Path
    
    import geckopy
    import numpy as np
    import pandas as pd
    import pytfa
    
    from cobra.flux_analysis.variability import flux_variability_analysis
    from geckopy.integration.pytfa import (
        adapt_gecko_to_thermo,
        translate_model_mnx_to_seed,
    )
    from geckopy.experimental import from_copy_number
    from geckopy.experimental.relaxation import apply_proteomics_relaxation
    from geckopy.experimental.molecular_weights import extract_proteins
    
    from pytfa.io import load_thermoDB
    from pytfa.io.plotting import plot_fva_tva_comparison
    from pytfa.optim.variables import DeltaG,DeltaGstd,ThermoDisplacement
    from pytfa.analysis import  variability_analysis

.. code:: ipython3

    DATA = Path.cwd().parent / "tests" / "data"

Load the model
~~~~~~~~~~~~~~

Load the enzyme constraint model and transform it to a thermo model
(taking into account the proteins).

.. code:: ipython3

    ec_model = geckopy.io.read_sbml_ec_model(DATA / "ec_coli_core.xml", hardcoded_rev_reactions=False)

We fix the glucose exchange rate so that the FVA runs are more
comparable between each other.

.. code:: ipython3

    ec_model.slim_optimize()
    ec_model.reactions.EX_glc__D_e.bounds = ec_model.reactions.EX_glc__D_e.flux, ec_model.reactions.EX_glc__D_e.flux

Geckopy to pytfa adapter. There is a translation of ids although in this
case the model already contains seed identifiers.

.. code:: ipython3

    thermodb = load_thermoDB(DATA / "thermo_data.thermodb")
    compartment_data = pytfa.io.read_compartment_data(
        str(DATA / "compartment_data.json")
    )
    # ranslate the mnx ids in the model to seed ids
    translate_model_mnx_to_seed(
        ec_model, thermodb, DATA / "chem_xref_seedset.tsv"
    )
    tmodel = adapt_gecko_to_thermo(
        ec_model, thermodb, compartment_data, solver="optlang-glpk"
    )
    tmodel.print_info()


.. parsed-literal::

    2021-05-21 17:30:42,353 - thermomodel_ - INFO - # Model initialized with units kcal/mol and temperature 298.15 K
    2021-05-21 17:30:42,358 - thermomodel_ - INFO - # Model preparation starting...
    2021-05-21 17:30:42,406 - thermomodel_ - INFO - # Model preparation done.
    2021-05-21 17:30:42,434 - thermomodel_ - INFO - # Model conversion starting...
    2021-05-21 17:30:43,134 - thermomodel_ - INFO - # Model conversion done.
    2021-05-21 17:30:43,135 - thermomodel_ - INFO - # Updating cobra_model variables...
    2021-05-21 17:30:43,139 - thermomodel_ - INFO - # cobra_model variables are up-to-date


.. parsed-literal::

                    value
    key                  
    name                 
    description          
    num constraints   631
    num variables     763
    num metabolites    72
    num reactions      95
                               value
    key                             
    num metabolites(thermo)      127
    num reactions(thermo)         73
    pct metabolites(thermo)  176.389
    pct reactions(thermo)    76.8421


Check that the model works.

.. code:: ipython3

    solution = tmodel.optimize()
    solution




.. raw:: html

    <strong><em>Optimal</em> solution with objective value 0.374</strong><br><div>
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
          <th>PFK</th>
          <td>9.627754</td>
          <td>None</td>
        </tr>
        <tr>
          <th>PFL</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>PGI</th>
          <td>9.923283</td>
          <td>None</td>
        </tr>
        <tr>
          <th>PGK</th>
          <td>-19.005185</td>
          <td>None</td>
        </tr>
        <tr>
          <th>PGL</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>...</th>
          <td>...</td>
          <td>...</td>
        </tr>
        <tr>
          <th>NADH16</th>
          <td>30.034558</td>
          <td>None</td>
        </tr>
        <tr>
          <th>NADTRHD</th>
          <td>0.000000</td>
          <td>None</td>
        </tr>
        <tr>
          <th>NH4t</th>
          <td>2.040601</td>
          <td>None</td>
        </tr>
        <tr>
          <th>O2t</th>
          <td>15.017279</td>
          <td>None</td>
        </tr>
        <tr>
          <th>PDH</th>
          <td>16.118563</td>
          <td>None</td>
        </tr>
      </tbody>
    </table>
    <p>95 rows × 2 columns</p>
    </div>



Flux variability analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate variability analysis on all continuous variables

Plain FVA
^^^^^^^^^

We use the FVA from pytfa for consistency here (which requires the model
to have and :code:`logger` property).

Since this model does not have splitted reactions (so that proteins are
reactants in both directions), we don’t need to take additional measures
to make the different runs comparable. See
:func:`~geckopy.flux_analysis.flux_variability_analysis` that takes
that into account if that is not the case.

.. code:: ipython3

    import logging

.. code:: ipython3

    # avoid loading the screen with logs
    logging.basicConfig(filename='docs03.log', level=logging.DEBUG)
    ec_model.logger = logging.getLogger(__name__)
    fva_fluxes = variability_analysis(ec_model.copy(), kind="reaction")


.. parsed-literal::

    minimizing: 100%|██████████| 95/95 [00:00<00:00, 810.72it/s]
    maximizing: 100%|██████████| 95/95 [00:00<00:00, 859.12it/s]


Thermo FVA
^^^^^^^^^^

.. code:: ipython3

    tva_fluxes = variability_analysis(tmodel, kind="reaction")


.. parsed-literal::

    2021-05-21 17:30:45,685 - thermomodel_ - INFO - Beginning variability analysis for variable of type reaction
    minimizing: 100%|██████████| 95/95 [00:27<00:00,  3.42it/s]
    maximizing: 100%|██████████| 95/95 [00:22<00:00,  4.26it/s]


Thermo concentration proteomics FVA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can apply both proteomics and thermodynamics (which may include
metabolomics) constraints.

This will usually require a relaxation of the experimental constraints
or the thermodynamic constraint, see the
`relaxation chapter <relaxation.html>`_ for a detailed explanation.

.. code:: ipython3

    raw_proteomics = pd.read_csv(DATA / "ecoli_proteomics_schmidt2016S5.tsv")
    ec_model_constrained = from_copy_number(
        ec_model.copy(),
        # the index should be the IDs of the proteins exactly as in the model!
        index=raw_proteomics["uniprot"].apply(lambda x: f"prot_{x}"),
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    tmodel_prot = adapt_gecko_to_thermo(
        ec_model_constrained, thermodb, compartment_data, solver="optlang-glpk"
    )
    # constrain the model objective so that the feashibility relaxation recovers growth
    tmodel_prot.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = solution.objective_value
    iis, status = geckopy.integration.relaxation.relax_thermo_proteins(
        tmodel_prot, 
        prot_candidates=[prot.id for prot in tmodel_prot.proteins], 
        objective_rule=geckopy.experimental.relaxation.Objective_rule.MIN_ELASTIC_SUM
    )
    print(f"IIS (status '{status}'): {iis}")


.. parsed-literal::

    2021-05-21 17:31:36,086 - thermomodel_ - INFO - # Model initialized with units kcal/mol and temperature 298.15 K
    2021-05-21 17:31:36,090 - thermomodel_ - INFO - # Model preparation starting...
    2021-05-21 17:31:36,138 - thermomodel_ - INFO - # Model preparation done.
    2021-05-21 17:31:36,165 - thermomodel_ - INFO - # Model conversion starting...
    2021-05-21 17:31:36,867 - thermomodel_ - INFO - # Model conversion done.
    2021-05-21 17:31:36,868 - thermomodel_ - INFO - # Updating cobra_model variables...
    2021-05-21 17:31:36,871 - thermomodel_ - INFO - # cobra_model variables are up-to-date
    adding thermo slacks: 100%|██████████| 73/73 [00:00<00:00, 308.11it/s]
    adding protein slacks: 100%|██████████| 55/55 [00:00<00:00, 421.43it/s]


.. parsed-literal::

    IIS (status 'optimal'): {'prot_P0AB71', 'prot_P00370', 'prot_P0A9P0', 'prot_P0A9N4', 'prot_P25516', 'prot_P21599', 'prot_P0A9C5', 'prot_P0A9B2', 'prot_P06999', 'prot_P37689', 'prot_P06959', 'prot_P0AFG8', 'prot_P0A6P9', 'prot_P33221', 'prot_P0A867'}


Since this relaxation is an in-place operation, we need to first rebuild
the model.

In this case, all of the problematic constraints are proteins, so we
will simply remove their concentrations for illustration purposes.

.. code:: ipython3

    ec_model_constrained = from_copy_number(
        ec_model.copy(),
        # the index should be the IDs of the proteins exactly as in the model!
        index=raw_proteomics["uniprot"].apply(lambda x: f"prot_{x}"),
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    tmodel_prot = adapt_gecko_to_thermo(
        ec_model_constrained, thermodb, compartment_data, solver="optlang-glpk"
    )
    for prot_id in iis:
        tmodel_prot.proteins.get_by_id(prot_id).add_concentration(None)


.. parsed-literal::

    2021-05-21 17:31:49,759 - thermomodel_ - INFO - # Model initialized with units kcal/mol and temperature 298.15 K
    2021-05-21 17:31:49,763 - thermomodel_ - INFO - # Model preparation starting...
    2021-05-21 17:31:49,810 - thermomodel_ - INFO - # Model preparation done.
    2021-05-21 17:31:49,838 - thermomodel_ - INFO - # Model conversion starting...
    2021-05-21 17:31:50,545 - thermomodel_ - INFO - # Model conversion done.
    2021-05-21 17:31:50,546 - thermomodel_ - INFO - # Updating cobra_model variables...
    2021-05-21 17:31:50,548 - thermomodel_ - INFO - # cobra_model variables are up-to-date


.. code:: ipython3

    tmodel_prot.slim_optimize()




.. parsed-literal::

    0.3742298749331099



.. code:: ipython3

    # Perform variability analysis again
    tva_fluxes_prot = variability_analysis(tmodel_prot, kind='reactions')


.. parsed-literal::

    2021-05-21 17:31:52,185 - thermomodel_ - INFO - Beginning variability analysis for variable of type reactions
    minimizing: 100%|██████████| 95/95 [00:40<00:00,  2.35it/s]
    maximizing: 100%|██████████| 95/95 [08:00<00:00,  5.06s/it]


Plotting
~~~~~~~~

Similar to `Figure 5 B of Sánchez et al.,
2017 <https://www.embopress.org/doi/full/10.15252/msb.20167411#msb167411-fig-0005>`__.

.. code:: ipython3

    import plotly.express as px

.. code:: ipython3

    dfs = (fva_fluxes, tva_fluxes, tva_fluxes_prot)

`Tidy <https://tidyr.tidyverse.org/articles/tidy-data.html>`__ up the
data to plot it with plotly express.

.. code:: ipython3

    df_plot = pd.concat([abs(df.maximum - df.minimum).apply(lambda x: x if x < 2000 else 2000) for df in dfs], axis=1).melt(var_name="Variability method", value_name="Flux")
    df_plot.loc[df_plot["Variability method"] == 0, "Variability method"] = "FBA"
    df_plot.loc[df_plot["Variability method"] == 1, "Variability method"] = "Thermo"
    df_plot.loc[df_plot["Variability method"] == 2, "Variability method"] = "Thermo + proteins"

.. code:: ipython3

    fig = px.histogram(
        df_plot, 
        x="Flux", color="Variability method",
        cumulative=True, nbins=200, 
        color_discrete_sequence=["rgba(26, 200, 26, 0.8)", "rgba(200, 26, 80, 0.5)", "rgba(26, 26, 200, 0.5)"],
        marginal="violin", barmode="overlay",
    )
    fig.show(renderer="notebook_connected")



.. raw:: html

    <script type="text/javascript">
    window.PlotlyConfig = {MathJaxConfig: 'local'};
    if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}
    if (typeof require !== 'undefined') {
    require.undef("plotly");
    requirejs.config({
        paths: {
            'plotly': ['https://cdn.plot.ly/plotly-latest.min']
        }
    });
    require(['plotly'], function(Plotly) {
        window._Plotly = Plotly;
    });
    }
    </script>



    <div>                            <div id="2fce539d-2e7b-4f2c-89f5-a78b827c537d" class="plotly-graph-div" style="height:525px; width:100%;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("2fce539d-2e7b-4f2c-89f5-a78b827c537d")) {                    Plotly.newPlot(                        "2fce539d-2e7b-4f2c-89f5-a78b827c537d",                        [{"alignmentgroup": "True", "bingroup": "x", "cumulative": {"enabled": true}, "hovertemplate": "Variability method=FBA<br>Flux=%{x}<br>count=%{y}<extra></extra>", "legendgroup": "FBA", "marker": {"color": "rgba(26, 200, 26, 0.8)", "opacity": 0.5}, "name": "FBA", "nbinsx": 200, "offsetgroup": "FBA", "orientation": "v", "showlegend": true, "type": "histogram", "x": [176.6099999999985, 39.99999999999987, 59.99999999999996, 20.000000000000068, 60.000000000000284, 19.999999999999755, 10.000000000000087, 20.000000000000323, 3.214895047684773, 19.99999999999985, 20.000000000000007, 20.000000000000053, 166.6099999999999, 19.999999999999993, 19.999999999999996, 166.6099999999999, 166.61000000000013, 20.0, 166.60999999999993, 166.60999999999993, 19.999999999999986, 181.61000000000018, 19.99999999999998, 176.60999999999973, 0.873921506968429, 20.00000000000002, 71.10424242424335, 40.62090900687144, 20.00000000000006, 20.000000000000107, 222.14666666666713, 120.00000000000041, 20.0, 19.99999999999986, 222.14666666666722, 19.999999999999982, 1000.0, 19.999999999999975, 20.154536201071043, 333.219999999998, 20.154536201070123, 20.466372805801072, 19.999999999999993, 19.99999999999988, 20.000000000000618, 10.000000000000199, 71.10424242424286, 19.9999999999995, 39.99999999999785, 0.0, 0.0, 0.0, 1.7409062416280233e-33, 10.000000000000076, 39.99999999999728, 69.99999999999915, 19.999999999999638, 0.0, 9.9999999999999, 60.00000000000011, 3.214895047684771, 19.999999999999254, 16.384166666667042, 19.999999999999996, 166.60999999998884, 666.4399999999572, 666.4399999999966, 1000.0, 0.0, 35.04114285714235, 0.0, 59.99999999999894, 19.999999999999844, 0.0, 166.60999999999746, 0.0, 176.60999999999746, 166.60999999999763, 166.60999999999746, 9.999999999999991, 59.99999999999901, 69.99999999999976, 20.00000000000009, 20.000000000000398, 19.999999999999986, 20.000000000000398, 0.0, 108.30500000000066, 98.3050000000018, 98.30500000000166, 120.00000000000178, 378.220000000007, 10.000000000000606, 59.99999999999923, 39.999999999999716], "xaxis": "x", "yaxis": "y"}, {"alignmentgroup": "True", "hovertemplate": "Variability method=FBA<br>Flux=%{x}<extra></extra>", "legendgroup": "FBA", "marker": {"color": "rgba(26, 200, 26, 0.8)"}, "name": "FBA", "offsetgroup": "FBA", "scalegroup": "x", "showlegend": false, "type": "violin", "x": [176.6099999999985, 39.99999999999987, 59.99999999999996, 20.000000000000068, 60.000000000000284, 19.999999999999755, 10.000000000000087, 20.000000000000323, 3.214895047684773, 19.99999999999985, 20.000000000000007, 20.000000000000053, 166.6099999999999, 19.999999999999993, 19.999999999999996, 166.6099999999999, 166.61000000000013, 20.0, 166.60999999999993, 166.60999999999993, 19.999999999999986, 181.61000000000018, 19.99999999999998, 176.60999999999973, 0.873921506968429, 20.00000000000002, 71.10424242424335, 40.62090900687144, 20.00000000000006, 20.000000000000107, 222.14666666666713, 120.00000000000041, 20.0, 19.99999999999986, 222.14666666666722, 19.999999999999982, 1000.0, 19.999999999999975, 20.154536201071043, 333.219999999998, 20.154536201070123, 20.466372805801072, 19.999999999999993, 19.99999999999988, 20.000000000000618, 10.000000000000199, 71.10424242424286, 19.9999999999995, 39.99999999999785, 0.0, 0.0, 0.0, 1.7409062416280233e-33, 10.000000000000076, 39.99999999999728, 69.99999999999915, 19.999999999999638, 0.0, 9.9999999999999, 60.00000000000011, 3.214895047684771, 19.999999999999254, 16.384166666667042, 19.999999999999996, 166.60999999998884, 666.4399999999572, 666.4399999999966, 1000.0, 0.0, 35.04114285714235, 0.0, 59.99999999999894, 19.999999999999844, 0.0, 166.60999999999746, 0.0, 176.60999999999746, 166.60999999999763, 166.60999999999746, 9.999999999999991, 59.99999999999901, 69.99999999999976, 20.00000000000009, 20.000000000000398, 19.999999999999986, 20.000000000000398, 0.0, 108.30500000000066, 98.3050000000018, 98.30500000000166, 120.00000000000178, 378.220000000007, 10.000000000000606, 59.99999999999923, 39.999999999999716], "xaxis": "x2", "yaxis": "y2"}, {"alignmentgroup": "True", "bingroup": "x", "cumulative": {"enabled": true}, "hovertemplate": "Variability method=Thermo<br>Flux=%{x}<br>count=%{y}<extra></extra>", "legendgroup": "Thermo", "marker": {"color": "rgba(200, 26, 80, 0.5)", "opacity": 0.5}, "name": "Thermo", "nbinsx": 200, "offsetgroup": "Thermo", "orientation": "v", "showlegend": true, "type": "histogram", "x": [176.61000000000013, 39.99999999999993, 45.00000000000018, 14.999999999999947, 45.000000000000206, 20.00000000000001, 10.000000000000004, 14.999999999999945, 1.3766794409164191, 20.00000000000001, 20.0, 20.0, 166.61000000000018, 19.999999999999968, 19.999999999999968, 166.61000000000013, 166.6100000000004, 20.000000000000004, 31.610000000000017, 31.610000000000007, 19.999999999999932, 177.8600000000011, 20.000000000000004, 41.61000000000003, 0.3742298749331099, 20.000000000000007, 71.104242424242, 30.26899643410217, 20.0, 15.000000000000057, 222.14666666666682, 119.99999999999949, 20.0, 14.999999999999943, 222.14666666666795, 20.0, 1000.0, 20.000000000000185, 15.066949724625502, 333.2199999999969, 15.066949724625502, 15.202046709476358, 15.000000000000057, 19.999999999999538, 20.0, 10.000000000000083, 71.10424242424246, 19.999999999999183, 39.999999999999794, 0.0, 0.0, 0.0, 0.0, 10.000000000000004, 39.9999999999999, 69.99999999999974, 20.000000000000025, 0.0, 10.0, 60.00000000000103, 1.376679440916424, 19.99999999999984, 15.402500000000225, 15.000000000000059, 166.61000000000195, 666.44, 666.4400000000119, 1000.0, 0.0, 33.639999999999844, 0.0, 45.00000000000008, 14.999999999999945, 0.0, 31.609999999999026, 0.0, 41.61000000000041, 31.609999999999694, 31.610000000000404, 10.000000000000004, 44.99999999999989, 69.99999999999986, 20.000000000000004, 20.00000000000013, 20.000000000000004, 20.0, 0.0, 60.80500000000006, 30.80500000000002, 30.80500000000002, 114.99999999999912, 378.219999999999, 10.0, 60.000000000000014, 39.999999999999815], "xaxis": "x", "yaxis": "y"}, {"alignmentgroup": "True", "hovertemplate": "Variability method=Thermo<br>Flux=%{x}<extra></extra>", "legendgroup": "Thermo", "marker": {"color": "rgba(200, 26, 80, 0.5)"}, "name": "Thermo", "offsetgroup": "Thermo", "scalegroup": "x", "showlegend": false, "type": "violin", "x": [176.61000000000013, 39.99999999999993, 45.00000000000018, 14.999999999999947, 45.000000000000206, 20.00000000000001, 10.000000000000004, 14.999999999999945, 1.3766794409164191, 20.00000000000001, 20.0, 20.0, 166.61000000000018, 19.999999999999968, 19.999999999999968, 166.61000000000013, 166.6100000000004, 20.000000000000004, 31.610000000000017, 31.610000000000007, 19.999999999999932, 177.8600000000011, 20.000000000000004, 41.61000000000003, 0.3742298749331099, 20.000000000000007, 71.104242424242, 30.26899643410217, 20.0, 15.000000000000057, 222.14666666666682, 119.99999999999949, 20.0, 14.999999999999943, 222.14666666666795, 20.0, 1000.0, 20.000000000000185, 15.066949724625502, 333.2199999999969, 15.066949724625502, 15.202046709476358, 15.000000000000057, 19.999999999999538, 20.0, 10.000000000000083, 71.10424242424246, 19.999999999999183, 39.999999999999794, 0.0, 0.0, 0.0, 0.0, 10.000000000000004, 39.9999999999999, 69.99999999999974, 20.000000000000025, 0.0, 10.0, 60.00000000000103, 1.376679440916424, 19.99999999999984, 15.402500000000225, 15.000000000000059, 166.61000000000195, 666.44, 666.4400000000119, 1000.0, 0.0, 33.639999999999844, 0.0, 45.00000000000008, 14.999999999999945, 0.0, 31.609999999999026, 0.0, 41.61000000000041, 31.609999999999694, 31.610000000000404, 10.000000000000004, 44.99999999999989, 69.99999999999986, 20.000000000000004, 20.00000000000013, 20.000000000000004, 20.0, 0.0, 60.80500000000006, 30.80500000000002, 30.80500000000002, 114.99999999999912, 378.219999999999, 10.0, 60.000000000000014, 39.999999999999815], "xaxis": "x2", "yaxis": "y2"}, {"alignmentgroup": "True", "bingroup": "x", "cumulative": {"enabled": true}, "hovertemplate": "Variability method=Thermo + proteins<br>Flux=%{x}<br>count=%{y}<extra></extra>", "legendgroup": "Thermo + proteins", "marker": {"color": "rgba(26, 26, 200, 0.5)", "opacity": 0.5}, "name": "Thermo + proteins", "nbinsx": 200, "offsetgroup": "Thermo + proteins", "orientation": "v", "showlegend": true, "type": "histogram", "x": [0.42720717115602724, 20.004420030439, 0.24098761270525415, 1.0462086970222764, 0.16453670753734664, 20.000000000000004, 1.1964639973576883, 1.6041138391256347, 1.37667944091645, 7.886263421906325, 20.000000000000007, 20.000000000000007, 1.3582238719885864, 1.5588315339439938, 1.5588315339439935, 81.72764469128707, 0.6481475621172655, 20.000000000000007, 1.8103912187073592, 1.8103912187073592, 1.1964639973576878, 81.89311822513909, 20.000000000000007, 4.9541545685247375, 0.37422987493311, 20.000000000000007, 21.434443650960432, 0.37336622552655496, 1.558831533943994, 0.32290854788492795, 108.9701929217161, 40.41158687711495, 0.0015504790096054688, 1.6041138391256382, 108.97019292171602, 7.886263421906324, 0.004420030439002111, 1.1964639973576845, 0.1204708795497516, 6.618722012872128, 0.12047087954975161, 0.2528953459768036, 0.4258002783035888, 20.000000000000007, 19.999999999999996, 1.196463997357688, 21.434443650960436, 7.886263421906324, 20.004420030439004, 0.0, 0.0, 0.0, 0.0, 1.362643902427589, 39.99999999999993, 29.516894704383173, 0.0015504790096054688, 0.0, 2.3011584020574443, 20.205793438557475, 1.376679440916444, 20.00000000000001, 1.5544115035049921, 0.4258002783035888, 0.0014068928524357733, 326.91057876514816, 326.91057876514816, 1.3626439024275887, 0.0, 1.3626439024275887, 0.0, 0.16453670753734667, 1.0462086970222693, 0.0, 31.610000000000028, 0.0, 2.2463291186006584, 31.610000000000014, 0.03792665251853772, 1.362643902427589, 0.16453670753734667, 29.51689470438317, 1.3626439024275887, 0.1961876315164189, 0.0015504790096054686, 0.19618763151640517, 0.0, 1.558831533943994, 0.0889026181127603, 1.1771318057855826, 40.407166846675956, 0.2021391542815347, 2.3011584020574447, 20.205793438557475, 20.004420030439007], "xaxis": "x", "yaxis": "y"}, {"alignmentgroup": "True", "hovertemplate": "Variability method=Thermo + proteins<br>Flux=%{x}<extra></extra>", "legendgroup": "Thermo + proteins", "marker": {"color": "rgba(26, 26, 200, 0.5)"}, "name": "Thermo + proteins", "offsetgroup": "Thermo + proteins", "scalegroup": "x", "showlegend": false, "type": "violin", "x": [0.42720717115602724, 20.004420030439, 0.24098761270525415, 1.0462086970222764, 0.16453670753734664, 20.000000000000004, 1.1964639973576883, 1.6041138391256347, 1.37667944091645, 7.886263421906325, 20.000000000000007, 20.000000000000007, 1.3582238719885864, 1.5588315339439938, 1.5588315339439935, 81.72764469128707, 0.6481475621172655, 20.000000000000007, 1.8103912187073592, 1.8103912187073592, 1.1964639973576878, 81.89311822513909, 20.000000000000007, 4.9541545685247375, 0.37422987493311, 20.000000000000007, 21.434443650960432, 0.37336622552655496, 1.558831533943994, 0.32290854788492795, 108.9701929217161, 40.41158687711495, 0.0015504790096054688, 1.6041138391256382, 108.97019292171602, 7.886263421906324, 0.004420030439002111, 1.1964639973576845, 0.1204708795497516, 6.618722012872128, 0.12047087954975161, 0.2528953459768036, 0.4258002783035888, 20.000000000000007, 19.999999999999996, 1.196463997357688, 21.434443650960436, 7.886263421906324, 20.004420030439004, 0.0, 0.0, 0.0, 0.0, 1.362643902427589, 39.99999999999993, 29.516894704383173, 0.0015504790096054688, 0.0, 2.3011584020574443, 20.205793438557475, 1.376679440916444, 20.00000000000001, 1.5544115035049921, 0.4258002783035888, 0.0014068928524357733, 326.91057876514816, 326.91057876514816, 1.3626439024275887, 0.0, 1.3626439024275887, 0.0, 0.16453670753734667, 1.0462086970222693, 0.0, 31.610000000000028, 0.0, 2.2463291186006584, 31.610000000000014, 0.03792665251853772, 1.362643902427589, 0.16453670753734667, 29.51689470438317, 1.3626439024275887, 0.1961876315164189, 0.0015504790096054686, 0.19618763151640517, 0.0, 1.558831533943994, 0.0889026181127603, 1.1771318057855826, 40.407166846675956, 0.2021391542815347, 2.3011584020574447, 20.205793438557475, 20.004420030439007], "xaxis": "x2", "yaxis": "y2"}],                        {"barmode": "overlay", "legend": {"title": {"text": "Variability method"}, "tracegroupgap": 0}, "margin": {"t": 60}, "template": {"data": {"bar": [{"error_x": {"color": "#2a3f5f"}, "error_y": {"color": "#2a3f5f"}, "marker": {"line": {"color": "#E5ECF6", "width": 0.5}}, "type": "bar"}], "barpolar": [{"marker": {"line": {"color": "#E5ECF6", "width": 0.5}}, "type": "barpolar"}], "carpet": [{"aaxis": {"endlinecolor": "#2a3f5f", "gridcolor": "white", "linecolor": "white", "minorgridcolor": "white", "startlinecolor": "#2a3f5f"}, "baxis": {"endlinecolor": "#2a3f5f", "gridcolor": "white", "linecolor": "white", "minorgridcolor": "white", "startlinecolor": "#2a3f5f"}, "type": "carpet"}], "choropleth": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "type": "choropleth"}], "contour": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "contour"}], "contourcarpet": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "type": "contourcarpet"}], "heatmap": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "heatmap"}], "heatmapgl": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "heatmapgl"}], "histogram": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "histogram"}], "histogram2d": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "histogram2d"}], "histogram2dcontour": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "histogram2dcontour"}], "mesh3d": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "type": "mesh3d"}], "parcoords": [{"line": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "parcoords"}], "pie": [{"automargin": true, "type": "pie"}], "scatter": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatter"}], "scatter3d": [{"line": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatter3d"}], "scattercarpet": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scattercarpet"}], "scattergeo": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scattergeo"}], "scattergl": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scattergl"}], "scattermapbox": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scattermapbox"}], "scatterpolar": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatterpolar"}], "scatterpolargl": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatterpolargl"}], "scatterternary": [{"marker": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "type": "scatterternary"}], "surface": [{"colorbar": {"outlinewidth": 0, "ticks": ""}, "colorscale": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "type": "surface"}], "table": [{"cells": {"fill": {"color": "#EBF0F8"}, "line": {"color": "white"}}, "header": {"fill": {"color": "#C8D4E3"}, "line": {"color": "white"}}, "type": "table"}]}, "layout": {"annotationdefaults": {"arrowcolor": "#2a3f5f", "arrowhead": 0, "arrowwidth": 1}, "autotypenumbers": "strict", "coloraxis": {"colorbar": {"outlinewidth": 0, "ticks": ""}}, "colorscale": {"diverging": [[0, "#8e0152"], [0.1, "#c51b7d"], [0.2, "#de77ae"], [0.3, "#f1b6da"], [0.4, "#fde0ef"], [0.5, "#f7f7f7"], [0.6, "#e6f5d0"], [0.7, "#b8e186"], [0.8, "#7fbc41"], [0.9, "#4d9221"], [1, "#276419"]], "sequential": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]], "sequentialminus": [[0.0, "#0d0887"], [0.1111111111111111, "#46039f"], [0.2222222222222222, "#7201a8"], [0.3333333333333333, "#9c179e"], [0.4444444444444444, "#bd3786"], [0.5555555555555556, "#d8576b"], [0.6666666666666666, "#ed7953"], [0.7777777777777778, "#fb9f3a"], [0.8888888888888888, "#fdca26"], [1.0, "#f0f921"]]}, "colorway": ["#636efa", "#EF553B", "#00cc96", "#ab63fa", "#FFA15A", "#19d3f3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"], "font": {"color": "#2a3f5f"}, "geo": {"bgcolor": "white", "lakecolor": "white", "landcolor": "#E5ECF6", "showlakes": true, "showland": true, "subunitcolor": "white"}, "hoverlabel": {"align": "left"}, "hovermode": "closest", "mapbox": {"style": "light"}, "paper_bgcolor": "white", "plot_bgcolor": "#E5ECF6", "polar": {"angularaxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}, "bgcolor": "#E5ECF6", "radialaxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}}, "scene": {"xaxis": {"backgroundcolor": "#E5ECF6", "gridcolor": "white", "gridwidth": 2, "linecolor": "white", "showbackground": true, "ticks": "", "zerolinecolor": "white"}, "yaxis": {"backgroundcolor": "#E5ECF6", "gridcolor": "white", "gridwidth": 2, "linecolor": "white", "showbackground": true, "ticks": "", "zerolinecolor": "white"}, "zaxis": {"backgroundcolor": "#E5ECF6", "gridcolor": "white", "gridwidth": 2, "linecolor": "white", "showbackground": true, "ticks": "", "zerolinecolor": "white"}}, "shapedefaults": {"line": {"color": "#2a3f5f"}}, "ternary": {"aaxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}, "baxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}, "bgcolor": "#E5ECF6", "caxis": {"gridcolor": "white", "linecolor": "white", "ticks": ""}}, "title": {"x": 0.05}, "xaxis": {"automargin": true, "gridcolor": "white", "linecolor": "white", "ticks": "", "title": {"standoff": 15}, "zerolinecolor": "white", "zerolinewidth": 2}, "yaxis": {"automargin": true, "gridcolor": "white", "linecolor": "white", "ticks": "", "title": {"standoff": 15}, "zerolinecolor": "white", "zerolinewidth": 2}}}, "xaxis": {"anchor": "y", "domain": [0.0, 1.0], "title": {"text": "Flux"}}, "xaxis2": {"anchor": "y2", "domain": [0.0, 1.0], "matches": "x", "showgrid": true, "showticklabels": false}, "yaxis": {"anchor": "x", "domain": [0.0, 0.7326], "title": {"text": "count"}}, "yaxis2": {"anchor": "x2", "domain": [0.7426, 1.0], "matches": "y2", "showgrid": false, "showline": false, "showticklabels": false, "ticks": ""}},                        {"responsive": true}                    ).then(function(){
    
    var gd = document.getElementById('2fce539d-2e7b-4f2c-89f5-a78b827c537d');
    var x = new MutationObserver(function (mutations, observer) {{
            var display = window.getComputedStyle(gd).display;
            if (!display || display === 'none') {{
                console.log([gd, 'removed!']);
                Plotly.purge(gd);
                observer.disconnect();
            }}
    }});
    
    // Listen for the removal of the full notebook cells
    var notebookContainer = gd.closest('#notebook-container');
    if (notebookContainer) {{
        x.observe(notebookContainer, {childList: true});
    }}
    
    // Listen for the clearing of the current output cell
    var outputEl = gd.closest('.output');
    if (outputEl) {{
        x.observe(outputEl, {childList: true});
    }}
    
                            })                };                });            </script>        </div>

