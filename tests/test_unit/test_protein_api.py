"""Tests related to `.proteins` API."""


def test_all_proteins_are_parsed(ec_model):
    """Test .proteins interface."""
    assert len(ec_model.proteins) == 1259


def test_protein_point_to_reactions(ec_model):
    """Test protein.reactions interface."""
    assert (
        ec_model.reactions.get_by_id("PUACGAMtrNo1")
        in ec_model.proteins.prot_P75905.reactions
    )
    assert len(ec_model.proteins.prot_P75905.reactions) == 3


def test_reaction_point_to_proteins(ec_model):
    """Test protein.reactions.REAC_ID.proteins interface."""
    assert (
        ec_model.proteins.get_by_id("prot_P75905")
        in ec_model.reactions.PUACGAMtrNo1.proteins
    )


def test_protein_pseudoreactions_are_not_model_reactions(ec_model):
    """Ensure that model.reactions is not polluted with protein pseudorreactions."""
    assert not ec_model.reactions.query("prot_P75905_exchange")


def test_protein_are_not_model_metabolites(ec_model):
    """Ensure that model.metabolites is not polluted with proteins."""
    assert not ec_model.metabolites.query("prot_P75905")
