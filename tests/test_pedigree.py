from repare.pedigree import Pedigree


def _add_node(
    pedigree: Pedigree,
    node_id: str,
    sex: str,
    y_haplogroup: str,
    mt_haplogroup: str,
    years_before_present: float = 0,
) -> None:
    pedigree.add_node(
        node_id=node_id,
        sex=sex,
        y_haplogroup=y_haplogroup,
        mt_haplogroup=mt_haplogroup,
        can_have_children=True,
        can_be_inbred=True,
        years_before_present=years_before_present,
    )


def test_fill_node_parents_creates_shared_placeholder_parents():
    pedigree = Pedigree()
    _add_node(pedigree, "A", "M", "Y1", "MT1")
    _add_node(pedigree, "B", "F", "", "MT2")
    pedigree.add_sibling_relation("A", "B")

    pedigree.fill_node_parents("A")

    assert pedigree.num_placeholders == 2
    assert pedigree.get_father("A") == pedigree.get_father("B") == "0"
    assert pedigree.get_mother("A") == pedigree.get_mother("B") == "1"
    assert pedigree.get_children("0") == {"A", "B"}
    assert pedigree.get_children("1") == {"A", "B"}
    assert pedigree.get_data("0")["sex"] == "M"
    assert pedigree.get_data("0")["y_haplogroup"] == "*"
    assert pedigree.get_data("1")["sex"] == "F"
    assert pedigree.get_data("1")["y_haplogroup"] == ""


def test_add_parent_relation_keeps_siblings_when_reassigning_same_parent():
    pedigree = Pedigree()
    _add_node(pedigree, "father", "M", "Y0", "MTF")
    _add_node(pedigree, "mother", "F", "", "MTM")
    _add_node(pedigree, "child1", "M", "Y1", "MT1")
    _add_node(pedigree, "child2", "F", "", "MT2")

    for parent in ("father", "mother"):
        for child in ("child1", "child2"):
            pedigree.add_parent_relation(parent, child)

    assert pedigree.get_siblings("child1") == {"child2"}
    assert pedigree.get_siblings("child2") == {"child1"}

    pedigree.add_parent_relation("father", "child1")

    assert pedigree.get_siblings("child1") == {"child2"}
    assert pedigree.get_siblings("child2") == {"child1"}
    assert pedigree.get_children("father") == {"child1", "child2"}
    assert pedigree.get_children("mother") == {"child1", "child2"}


def test_update_haplogroups_propagates_specific_lineage_information():
    pedigree = Pedigree()
    _add_node(pedigree, "father", "M", "A*", "*")
    _add_node(pedigree, "mother", "F", "", "H*")
    _add_node(pedigree, "son", "M", "A1", "H1")
    _add_node(pedigree, "daughter", "F", "", "H1")
    pedigree.add_parent_relation("father", "son")
    pedigree.add_parent_relation("mother", "son")
    pedigree.add_parent_relation("father", "daughter")
    pedigree.add_parent_relation("mother", "daughter")

    pedigree.update_haplogroups()

    assert pedigree.get_data("father")["y_haplogroup"] == "A1*"
    assert pedigree.get_data("mother")["mt_haplogroup"] == "H1*"
    assert pedigree.get_data("son")["y_haplogroup"] == "A1"
    assert pedigree.get_data("daughter")["mt_haplogroup"] == "H1"


def test_validate_years_before_present_detects_inconsistent_ancestors():
    pedigree = Pedigree()
    _add_node(pedigree, "father", "M", "Y1", "MTF", years_before_present=80)
    _add_node(pedigree, "mother", "F", "", "MTM", years_before_present=120)
    _add_node(pedigree, "child", "F", "", "MTM", years_before_present=100)
    pedigree.add_parent_relation("father", "child")
    pedigree.add_parent_relation("mother", "child")

    assert not pedigree.validate_years_before_present()

    pedigree.get_data("father")["years_before_present"] = 150
    assert pedigree.validate_years_before_present()


def test_is_relation_in_pedigree_distinguishes_maternal_and_paternal_paths():
    pedigree = Pedigree()
    # Maternal side
    _add_node(pedigree, "grandma_m", "F", "", "MTM")
    _add_node(pedigree, "grandpa_m", "M", "YM", "")
    _add_node(pedigree, "aunt", "F", "", "MTM")
    _add_node(pedigree, "mother", "F", "", "MTM")
    # Paternal side
    _add_node(pedigree, "grandma_p", "F", "", "MTP")
    _add_node(pedigree, "grandpa_p", "M", "YP", "")
    _add_node(pedigree, "paternal_aunt", "F", "", "MTP")
    _add_node(pedigree, "father", "M", "YP", "")
    # Child
    _add_node(pedigree, "child", "M", "YP", "MTM")

    for parent in ("grandma_m", "grandpa_m"):
        pedigree.add_parent_relation(parent, "aunt")
        pedigree.add_parent_relation(parent, "mother")
    for parent in ("grandma_p", "grandpa_p"):
        pedigree.add_parent_relation(parent, "paternal_aunt")
        pedigree.add_parent_relation(parent, "father")

    pedigree.add_parent_relation("mother", "child")
    pedigree.add_parent_relation("father", "child")

    assert pedigree.is_relation_in_pedigree("aunt", "child", ["maternal aunt/uncle-nephew/niece"])
    assert not pedigree.is_relation_in_pedigree("aunt", "child", ["paternal aunt/uncle-nephew/niece"])
    assert pedigree.is_relation_in_pedigree("grandma_m", "child", ["maternal grandparent-grandchild"])
    assert not pedigree.is_relation_in_pedigree("grandma_p", "child", ["maternal grandparent-grandchild"])
    assert pedigree.is_relation_in_pedigree("paternal_aunt", "child", ["paternal aunt/uncle-nephew/niece"])


def test_is_relation_in_pedigree_identifies_double_cousins():
    pedigree = Pedigree()
    # Grandparents for each sibling pair
    _add_node(pedigree, "pgf1", "M", "Y1", "")
    _add_node(pedigree, "pgm1", "F", "", "MT1")
    _add_node(pedigree, "pgf2", "M", "Y2", "")
    _add_node(pedigree, "pgm2", "F", "", "MT2")

    _add_node(pedigree, "paternal_parent", "M", "Y1", "MT1")
    _add_node(pedigree, "paternal_sibling", "F", "", "MT1")
    _add_node(pedigree, "maternal_parent", "F", "", "MT2")
    _add_node(pedigree, "maternal_sibling", "M", "Y2", "")

    _add_node(pedigree, "child_one", "M", "Y1", "MT2")
    _add_node(pedigree, "child_two", "F", "", "MT2")

    for parent in ("pgf1", "pgm1"):
        pedigree.add_parent_relation(parent, "paternal_parent")
        pedigree.add_parent_relation(parent, "paternal_sibling")
    for parent in ("pgf2", "pgm2"):
        pedigree.add_parent_relation(parent, "maternal_parent")
        pedigree.add_parent_relation(parent, "maternal_sibling")

    pedigree.add_parent_relation("paternal_parent", "child_one")
    pedigree.add_parent_relation("maternal_parent", "child_one")
    pedigree.add_parent_relation("paternal_sibling", "child_two")
    pedigree.add_parent_relation("maternal_sibling", "child_two")

    assert pedigree.is_relation_in_pedigree("child_one", "child_two", ["double cousins"])


def test_is_relation_in_pedigree_detects_maternal_half_siblings():
    pedigree = Pedigree()
    _add_node(pedigree, "shared_mother", "F", "", "MT1")
    _add_node(pedigree, "father_a", "M", "YA", "")
    _add_node(pedigree, "father_b", "M", "YB", "")
    _add_node(pedigree, "child_a", "M", "YA", "MT1")
    _add_node(pedigree, "child_b", "F", "", "MT1")

    pedigree.add_parent_relation("shared_mother", "child_a")
    pedigree.add_parent_relation("father_a", "child_a")
    pedigree.add_parent_relation("shared_mother", "child_b")
    pedigree.add_parent_relation("father_b", "child_b")

    assert pedigree.is_relation_in_pedigree("child_a", "child_b", ["maternal half-siblings"])
    assert not pedigree.is_relation_in_pedigree("child_a", "child_b", ["paternal half-siblings"])
