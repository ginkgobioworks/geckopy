Relaxation
==========

Oftentimes, appying experimental constraints such as proteomics and metabolomics
leave the model in a infeashible state. To recover from the infeashibility,
geckopy provides different relaxation methods that compute an Irreducibly Inconsistent Set (IIS); i.e., a minimal set of infeashible constraints.

The IIS may not be unique, so we need to set a criterion to select the best IIS.
This criterion - the objective of the relaxation problem - is implemented in the
form of a enum :class:`~geckopy.experimental.relaxation.ObjectiveRule`
that all of the relaxation functions of geckopy accept as a parameter
(:code:`objective_rule`):

- Objective_rule.MIN_ELASTIC_SUM: :math:`\sum_{v \in \text{elastic vars}} v_{flux}` (LP).
- Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE: :math:`\sum_{e \in \text{elastic vars}} + \text{prev objective}` (LP).
- Objective_rule.MIN_MILP_COUNT: :math:`\sum_i^{N} e_i` where e is a binary variable (MILP).

Relaxation functions on :class:`geckopy.Model`
----------------------------------------------

- :func:`~geckopy.experimental.relaxation.get_upper_relaxation` builds the
  relaxation problem (adding elastic variables to protein constraints), sets
  the objective and solves, returning the IIS and the status of the solver.
- :func:`~geckopy.experimental.relaxation.elastic_upper_relaxation` runs
  `get_upper_relaxation` iteratively, removing previous found infeashible
  variables, to compute a wider IIS (although it might still be non-unique). It
  is not an in-place operation and returns both a relaxed model and the IIS.

Relaxation functions on thermodynamic :class:`geckopy.Model`
------------------------------------------------------------

- :func:`~geckopy.integration.relaxation.relax_thermo_proteins` mimics 
  :func:`~geckopy.experimental.relaxation.get_upper_relaxation` but also adds
  :math:`\Delta G_r` relaxation variables to take into account thermodynamic
  constraints as source of the infeashibility.
- :func:`~geckopy.integration.relaxation.relax_thermo_concentrations_proteins`
  mimics :func:`~geckopy.experimental.relaxation.get_upper_relaxation` but also
  adds relaxation variables to take into account metabolomics measurements
  as source of the infeashibility.
