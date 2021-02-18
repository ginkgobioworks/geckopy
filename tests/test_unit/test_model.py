"""Tests related to `.proteins` API."""


def test_constraint_pool_changes_objective_value(ec_model):
    """Test .proteins interface."""
    dry_solution = ec_model.slim_optimize()
    # mw are not included in the model, let's say all of them are 330 Da
    for prot in ec_model.proteins:
        prot.mw = 330
    ec_model.constrain_pool(2000000000, 0.8, 1.0)
    pool_solution = ec_model.slim_optimize()
    assert dry_solution == pool_solution
