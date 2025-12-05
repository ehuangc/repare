from collections import defaultdict
from pathlib import Path

import pytest

from benchmarks.published.evaluator.pedigree_evaluator import PedigreeEvaluator


class StubPedigree:
    def __init__(self, relations: dict[tuple[str, str], dict[str, int]], inconsistency_count: int = 0):
        self._relations = {}
        nodes: set[str] = set()
        for pair, relation_counts in relations.items():
            key = tuple(sorted(pair))
            nodes.update(key)
            self._relations[key] = defaultdict(int, relation_counts)
        self.node_to_data = {node: {} for node in sorted(nodes)}
        self._inconsistency_count = inconsistency_count
        self.count_inconsistencies_calls: list[tuple[defaultdict, defaultdict, bool]] = []

    def get_relations_between_nodes(self, id1: str, id2: str, include_maternal_paternal: bool = True):
        return defaultdict(int, self._relations.get(tuple(sorted((id1, id2))), {}))

    def count_inconsistencies(
        self,
        pair_to_constraints,
        pair_to_relations_so_far,
        check_half_siblings: bool,
    ):
        self.count_inconsistencies_calls.append((pair_to_constraints, pair_to_relations_so_far, check_half_siblings))
        return self._inconsistency_count, []


def _write_csv(path: Path, header: list[str], rows: list[tuple[str, ...]]) -> None:
    lines = ["# Test data", ",".join(header)]
    for row in rows:
        lines.append(",".join(str(value) for value in row))
    path.write_text("\n".join(lines) + "\n")


def _instantiate_evaluator(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    relations: dict[tuple[str, str], dict[str, int]],
    published_rows: list[tuple[str, ...]],
    algorithm_rows: list[tuple[str, ...]],
) -> "PedigreeEvaluator":
    published_path = tmp_path / "published.csv"
    algorithm_nodes_path = tmp_path / "nodes.csv"
    algorithm_relations_path = tmp_path / "inferred.csv"

    _write_csv(published_path, ["id1", "id2", "degree", "constraints", "force_constraints"], published_rows)
    _write_csv(algorithm_relations_path, ["id1", "id2", "degree", "constraints"], algorithm_rows)

    nodes = sorted({node for pair in relations for node in pair})
    _write_csv(algorithm_nodes_path, ["id", "sex"], [(node, "M") for node in nodes])

    stub_pedigree = StubPedigree(relations)
    monkeypatch.setattr(
        PedigreeEvaluator,
        "_run_algorithm",
        staticmethod(lambda nodes_path, relations_path, _stub=stub_pedigree: _stub),
    )

    return PedigreeEvaluator(str(published_path), str(algorithm_nodes_path), str(algorithm_relations_path))


def test_sort_relation_flips_directional_relation_when_ids_are_unsorted():
    id1, id2, relation = PedigreeEvaluator._sort_relation("B", "A", "parent-child")

    assert (id1, id2, relation) == ("A", "B", "child-parent")


def test_sort_relation_leaves_symmetric_relation_unchanged():
    id1, id2, relation = PedigreeEvaluator._sort_relation("A", "B", "siblings")

    assert (id1, id2, relation) == ("A", "B", "siblings")


def test_get_metrics_reports_perfect_scores(monkeypatch, tmp_path):
    relations = {
        ("A", "B"): {"parent-child": 1},
        ("B", "C"): {"siblings": 1},
        ("A", "C"): {"maternal grandparent-grandchild": 1},
    }
    published_rows = [
        ("A", "B", "1", "parent-child", ""),
        ("B", "C", "1", "siblings", ""),
        ("A", "C", "2", "maternal grandparent-grandchild", ""),
    ]
    algorithm_rows = [
        ("A", "B", "1", "parent-child;child-parent"),
        ("B", "C", "1", "siblings"),
        ("A", "C", "2", "maternal grandparent-grandchild"),
    ]

    evaluator = _instantiate_evaluator(tmp_path, monkeypatch, relations, published_rows, algorithm_rows)

    metrics = evaluator.get_metrics()

    assert metrics["Pairwise Relation Accuracy"] == pytest.approx(1.0)
    assert metrics["Relation Precision"] == pytest.approx(1.0)
    assert metrics["Relation Recall"] == pytest.approx(1.0)
    assert metrics["Relation F1"] == pytest.approx(1.0)
    assert metrics["Pairwise Degree Accuracy"] == pytest.approx(1.0)
    assert metrics["Degree Precision"] == pytest.approx(1.0)
    assert metrics["Degree Recall"] == pytest.approx(1.0)
    assert metrics["Degree F1"] == pytest.approx(1.0)
    assert metrics["Connectivity R-squared"] == pytest.approx(1.0)


def test_get_metrics_detects_relation_and_degree_errors(monkeypatch, tmp_path):
    relations = {
        ("A", "B"): {"maternal grandparent-grandchild": 1},
        ("B", "C"): {"siblings": 1},
        ("A", "C"): {"maternal grandparent-grandchild": 1},
    }
    published_rows = [
        ("A", "B", "1", "parent-child", ""),
        ("B", "C", "1", "siblings", ""),
        ("A", "C", "2", "maternal grandparent-grandchild", ""),
    ]
    algorithm_rows = [
        ("A", "B", "2", "maternal grandparent-grandchild"),
        ("B", "C", "1", "siblings"),
        ("A", "C", "2", "maternal grandparent-grandchild"),
    ]

    evaluator = _instantiate_evaluator(tmp_path, monkeypatch, relations, published_rows, algorithm_rows)
    metrics = evaluator.get_metrics()

    two_thirds = pytest.approx(2 / 3)
    assert metrics["Pairwise Relation Accuracy"] == two_thirds
    assert metrics["Relation Precision"] == two_thirds
    assert metrics["Relation Recall"] == two_thirds
    assert metrics["Relation F1"] == two_thirds
    assert metrics["Pairwise Degree Accuracy"] == two_thirds
    assert metrics["Degree Precision"] == two_thirds
    assert metrics["Degree Recall"] == two_thirds
    assert metrics["Degree F1"] == two_thirds
    assert metrics["Connectivity R-squared"] == pytest.approx(1.0)


def test_count_input_relation_inconsistencies_reports_counts(monkeypatch, tmp_path):
    relations = {("A", "B"): {"parent-child": 1}}
    published_rows = [("A", "B", "1", "parent-child", "")]
    algorithm_rows = [("A", "B", "1", "parent-child")]

    evaluator = _instantiate_evaluator(tmp_path, monkeypatch, relations, published_rows, algorithm_rows)
    published_stub = StubPedigree({}, inconsistency_count=5)
    inferred_stub = StubPedigree({}, inconsistency_count=2)
    evaluator._published_pedigree = published_stub
    evaluator.algorithm_pedigree = inferred_stub

    metrics = evaluator.count_input_relation_inconsistencies()

    assert metrics["Published Pedigree Input Inconsistencies"] == 5
    assert metrics["Inferred Pedigree Input Inconsistencies"] == 2
    assert published_stub.count_inconsistencies_calls[0][2] is True
    assert inferred_stub.count_inconsistencies_calls[0][2] is True
