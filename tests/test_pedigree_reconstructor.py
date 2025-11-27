import math
from pathlib import Path
from unittest import mock

from repare.pedigree import Pedigree
from repare.pedigree_reconstructor import PedigreeReconstructor


def _write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    lines = [",".join(header)]
    lines.extend(",".join(row) for row in rows)
    path.write_text("\n".join(lines))


class _NoOpTqdm:
    def __init__(self, iterable, *_, **__):
        self._iterable = iterable

    def __iter__(self):
        return iter(self._iterable)

    def set_description(self, *_args, **_kwargs):  # pragma: no cover - trivial shim
        pass


def _run_reconstruction(reconstructor: PedigreeReconstructor) -> Pedigree:
    with mock.patch("repare.pedigree_reconstructor.tqdm", _NoOpTqdm):
        return reconstructor.find_best_pedigree()


def _build_reconstructor(
    tmp_path: Path,
    nodes_rows: list[list[str]] | None = None,
    relation_rows: list[list[str]] | None = None,
) -> PedigreeReconstructor:
    nodes_path = tmp_path / "nodes.csv"
    relations_path = tmp_path / "relations.csv"
    output_dir = tmp_path / "outputs"
    output_dir.mkdir()

    default_nodes = [
        ["Alice", "F", "", "H1*", "", "", "150"],
        ["Bob", "M", "A1", "M1", "False", "True", ""],
        ["Charlie", "M", "A1", "T1", "True", "False", "80"],
    ]
    default_relations = [
        ["Bob", "Alice", "1", "parent-child", ""],
        ["Charlie", "Alice", "2", "", ""],
        ["Charlie", "Bob", "3", "", "True"],
        ["Alice", "Bob", "2", "", ""],
    ]

    _write_csv(
        nodes_path,
        [
            "id",
            "sex",
            "y_haplogroup",
            "mt_haplogroup",
            "can_have_children",
            "can_be_inbred",
            "years_before_present",
        ],
        nodes_rows or default_nodes,
    )

    _write_csv(
        relations_path,
        ["id1", "id2", "degree", "constraints", "force_constraints"],
        relation_rows or default_relations,
    )

    return PedigreeReconstructor(
        relations_path=relations_path,
        nodes_path=nodes_path,
        outputs_dir=output_dir,
        plot=False,
        write_alternate_pedigrees=False,
    )


def _build_family_pedigree(include_maternal_sibling: bool) -> Pedigree:
    pedigree = Pedigree()
    pedigree.add_node("Alice", "F", "", "H1", True, True, 0)
    pedigree.add_node("Ben", "M", "Y1", "H2", True, True, 0)
    pedigree.add_node("Cara", "F", "", "H3", True, True, 0)
    pedigree.add_node("Dora", "F", "", "H4", True, True, 0)
    pedigree.add_parent_relation("Alice", "Cara")
    pedigree.add_parent_relation("Ben", "Cara")
    if include_maternal_sibling:
        pedigree.add_sibling_relation("Alice", "Dora")
    return pedigree


def test_process_data_normalizes_node_and_relation_inputs(tmp_path):
    reconstructor = _build_reconstructor(tmp_path)

    node_data = reconstructor._node_data.set_index("id")
    assert bool(node_data.loc["Alice", "can_have_children"]) is True
    assert bool(node_data.loc["Alice", "can_be_inbred"]) is True
    assert bool(node_data.loc["Bob", "can_have_children"]) is False
    assert bool(node_data.loc["Bob", "can_be_inbred"]) is True
    assert math.isnan(node_data.loc["Bob", "years_before_present"])
    assert node_data.loc["Alice", "years_before_present"] == 150

    relation_data = reconstructor._relation_data
    first_degree = relation_data[
        (relation_data["id1"] == "Alice") & (relation_data["id2"] == "Bob") & (relation_data["degree"] == "1")
    ].iloc[0]
    assert first_degree["constraints"] == "child-parent"
    assert bool(first_degree["force_constraints"]) is False

    second_degree = relation_data[
        (relation_data["id1"] == "Alice") & (relation_data["id2"] == "Charlie") & (relation_data["degree"] == "2")
    ].iloc[0]
    assert second_degree["constraints"] == reconstructor._DEFAULT_CONSTRAINTS["2"]
    assert bool(second_degree["force_constraints"]) is False

    third_degree = relation_data[
        (relation_data["id1"] == "Bob") & (relation_data["id2"] == "Charlie") & (relation_data["degree"] == "3")
    ].iloc[0]
    assert bool(third_degree["force_constraints"]) is True


def test_pair_to_constraints_prioritizes_specific_relations(tmp_path):
    reconstructor = _build_reconstructor(tmp_path)
    pair_constraints = reconstructor._pair_to_constraints[("Alice", "Bob")]
    lengths = [len(constraints) for constraints in pair_constraints]
    assert lengths == sorted(lengths)
    assert pair_constraints[0] == ("child-parent",)
    assert pair_constraints[1] == tuple(reconstructor._DEFAULT_CONSTRAINTS["2"].split(";"))


def test_check_parent_child_haplogroups_enforces_lineage_consistency():
    pedigree = Pedigree()
    pedigree.add_node("father", "M", "A1", "H1", True, True, 0)
    pedigree.add_node("son", "M", "B2", "H1", True, True, 0)
    pedigree.add_parent_relation("father", "son")
    assert not PedigreeReconstructor._check_parent_child_haplogroups(pedigree, "father", "son")

    pedigree.add_node("mother", "F", "", "H2", True, True, 0)
    pedigree.add_node("daughter", "F", "", "L2", True, True, 0)
    pedigree.add_parent_relation("mother", "daughter")
    assert not PedigreeReconstructor._check_parent_child_haplogroups(pedigree, "mother", "daughter")

    pedigree.add_node("compatible_father", "M", "Q1*", "H3", True, True, 0)
    pedigree.add_node("compatible_son", "M", "Q1A", "H4", True, True, 0)
    pedigree.add_parent_relation("compatible_father", "compatible_son")
    assert PedigreeReconstructor._check_parent_child_haplogroups(pedigree, "compatible_father", "compatible_son")

    pedigree.add_node("compatible_mother", "F", "", "J1*", True, True, 0)
    pedigree.add_node("compatible_daughter", "F", "", "J1B", True, True, 0)
    pedigree.add_parent_relation("compatible_mother", "compatible_daughter")
    assert PedigreeReconstructor._check_parent_child_haplogroups(pedigree, "compatible_mother", "compatible_daughter")


def test_find_best_pedigree_reconstructs_forced_relations(tmp_path):
    nodes = [
        ["Ivy", "F", "", "H1", "True", "True", "120"],
        ["Leo", "M", "Q1", "T2", "True", "True", "118"],
        ["Mia", "F", "", "H1", "True", "True", "80"],
    ]
    relations = [
        ["Ivy", "Mia", "1", "parent-child", "True"],
        ["Leo", "Mia", "1", "parent-child", "True"],
    ]
    reconstructor = _build_reconstructor(tmp_path, nodes_rows=nodes, relation_rows=relations)

    pedigree = _run_reconstruction(reconstructor)

    assert pedigree.get_mother("Mia") == "Ivy"
    assert pedigree.get_father("Mia") == "Leo"
    assert reconstructor._sample_strike_count == 0

    corrected_relations = (tmp_path / "outputs" / "corrected_input_relations.csv").read_text()
    assert "# Final inconsistency count: 0" in corrected_relations

    exact_relations = (tmp_path / "outputs" / "reconstructed_exact_relations.csv").read_text().splitlines()
    assert "Ivy,Mia,parent-child" in exact_relations[1:]
    assert "Leo,Mia,parent-child" in exact_relations[1:]


def test_find_best_pedigree_logs_removed_conflicting_relations(tmp_path):
    nodes = [
        ["Ada", "F", "", "H1", "True", "True", "120"],
        ["Ben", "M", "Q1", "H2", "True", "True", "118"],
        ["Cara", "F", "", "H1", "True", "False", "80"],
    ]
    relations = [
        ["Ada", "Cara", "1", "parent-child", "True"],
        ["Ben", "Cara", "1", "parent-child", "True"],
        ["Ada", "Ben", "2", "", ""],
    ]
    reconstructor = _build_reconstructor(tmp_path, nodes_rows=nodes, relation_rows=relations)

    pedigree = _run_reconstruction(reconstructor)

    assert pedigree.get_mother("Cara") == "Ada"
    assert pedigree.get_father("Cara") == "Ben"
    assert reconstructor._sample_strike_count == 1

    corrected_relations = (tmp_path / "outputs" / "corrected_input_relations.csv").read_text()
    assert "# Ada,Ben,2," in corrected_relations

    exact_relations = (tmp_path / "outputs" / "reconstructed_exact_relations.csv").read_text()
    assert "Ada,Cara,parent-child" in exact_relations
    assert "Ben,Cara,parent-child" in exact_relations


def test_write_constant_relations_emits_shared_pairs_only(tmp_path):
    reconstructor = _build_reconstructor(tmp_path)
    reconstructor._final_pedigrees = [
        _build_family_pedigree(include_maternal_sibling=True),
        _build_family_pedigree(include_maternal_sibling=False),
    ]

    constants_path = tmp_path / "outputs" / "constant_relations.csv"
    reconstructor._write_constant_relations(constants_path)

    output_lines = constants_path.read_text().splitlines()
    assert output_lines[1] == "node1,node2,relation"
    assert output_lines[2:] == [
        "Alice,Cara,parent-child",
        "Ben,Cara,parent-child",
    ]
    assert "Alice,Dora,siblings" not in output_lines


def test_connect_parent_relation_assigns_known_parent():
    pedigree = Pedigree()
    pedigree.add_node("mom", "F", "", "H1", True, True, 0)
    pedigree.add_node("dad", "M", "Y1", "", True, True, 0)
    pedigree.add_node("child", "M", "Y1", "H1", True, True, 0)

    result = PedigreeReconstructor._connect_parent_relation(pedigree, "mom", "child")

    assert len(result) == 1
    updated = result[0]
    assert updated.get_mother("child") == "mom"
    assert updated.get_children("mom") == {"child"}


def test_connect_sibling_relation_merges_parent_placeholders():
    pedigree = Pedigree()
    pedigree.add_node("alex", "M", "Y1", "H1", True, True, 0)
    pedigree.add_node("blair", "F", "", "H1", True, True, 0)

    result = PedigreeReconstructor._connect_sibling_relation(pedigree, "alex", "blair")

    assert len(result) == 1
    updated = result[0]
    assert updated.get_father("alex") == updated.get_father("blair")
    assert updated.get_mother("alex") == updated.get_mother("blair")
    assert updated.get_siblings("alex") == {"blair"}


def test_connect_second_degree_relation_adds_maternal_aunt_path():
    pedigree = Pedigree()
    pedigree.add_node("aunt", "F", "", "H1", True, True, 0)
    pedigree.add_node("mother", "F", "", "H1", True, True, 0)
    pedigree.add_node("father", "M", "Y1", "", True, True, 0)
    pedigree.add_node("child", "M", "Y1", "H1", True, True, 0)
    pedigree.add_parent_relation("mother", "child")
    pedigree.add_parent_relation("father", "child")

    result = PedigreeReconstructor._connect_second_degree_relation(
        pedigree, "aunt", "child", "maternal aunt/uncle-nephew/niece"
    )

    assert len(result) == 1
    updated = result[0]
    assert "aunt" in updated.get_siblings("mother")
    assert updated.is_relation_in_pedigree("aunt", "child", ["maternal aunt/uncle-nephew/niece"])
