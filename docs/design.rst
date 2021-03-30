Design principles
=================

- Proteins should be separated from metabolites and reactions at the user-level

.. code:: python

  model.metabolites
  model.reactions
  model.proteins

- Writable models that comply with community standard (SBML) **without hacks**. 
- Experimental data should be easy to load.

SBML compliance
~~~~~~~~~~~~~~~
Proteins are Species_:

* in Group_ (SBML extension) Proteins;
* with `initialAmount` (optional): concentration (double);
* and kcat value (double) in reaction stoichiometry.

Protein Exchanges:

* Not present in the SBML file.
* Generated at the time of parsing.

For backwards compatibility, with the available enzyme constrained models, 
:class:`~geckopy.io.read_sbml_ec_model` also parses ungroupped metabolites
with following the naming convention `prot_{UNIPROT_ID}`.

.. code:: xml

  ...
  <listOfSpecies>
      <species id="prot_P00363" compartment="c" initialAmount="0.12"/>
      <species id="prot_P08921" compartment="c" initialAmount="1.4"/>
  </listOfSpecies>
  <listOfGroups xmlns="http://www.sbml.org/sbml/level3/version1/groups/version1">
      <group kind="classification" id="Protein">
          <listOfMembers id="all_proteins">
              <member idRef="prot_P00363"/>
              <member idRef="prot_P08921"/>
          </listOfMembers>
      </group>
  </listOfGroups>
  ...


.. _Species: https://www.embopress.org/doi/epdf/10.15252/msb.20199110
.. _Group: http://europepmc.org/article/MED/28187406
