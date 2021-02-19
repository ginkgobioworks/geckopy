"""Test Input/Output capabilities."""
import geckopy


def test_read_geckopy_from_file(path_eciML1515):
    """Read model directly from file."""
    model = geckopy.io.read_sbml_ec_model(path_eciML1515)
    assert len(model.proteins) == 1259


def test_copy_geckopy(ec_model):
    """Check that deepcopy works."""
    copied = ec_model.copy()
    assert len(copied.proteins) == len(ec_model.proteins)
    assert len(copied.reactions) == len(ec_model.reactions)
    assert len(copied.metabolites) == len(ec_model.metabolites)


def test_parsing_captures_naming_convention(dummy_ec_model):
    """Check proteins rely on the naming convention prot_UNIPROT are parsed."""
    assert dummy_ec_model.proteins.query("prot_P0A805")


def test_parsing_captures_protein_group(dummy_ec_model):
    """Check members of Protein group are parsed as proteins."""
    assert dummy_ec_model.groups.query("Protein")
    assert dummy_ec_model.proteins.query("prot_P0A825")
    assert dummy_ec_model.proteins.query("dummy_prot")


def test_protein_parsing_does_not_get_normal_metabolites(dummy_ec_model):
    """Check normal metabolites are not parsed as proteins."""
    assert not dummy_ec_model.proteins.query("normal_met")
    assert dummy_ec_model.metabolites.query("normal_met")
    mets = set(dummy_ec_model.metabolites)
    prots = set(dummy_ec_model.proteins)
    assert mets ^ prots == mets | prots
