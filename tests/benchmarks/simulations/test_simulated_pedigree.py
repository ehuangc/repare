from collections import defaultdict

import pandas as pd
import pytest

from benchmarks.simulations.simulator.simulated_pedigree import SimulatedPedigree


class StubRelationsPedigree:
    def __init__(self, relations: dict[tuple[str, str], dict[str, int]]):
        self._relations = {}
        for pair, relation_counts in relations.items():
            key = tuple(sorted(pair))
            self._relations[key] = defaultdict(int, relation_counts)

    def get_relations_between_nodes(self, node1: str, node2: str, include_maternal_paternal: bool = True):
        return defaultdict(int, self._relations.get(tuple(sorted((node1, node2))), {}))


def _build_simulated_pedigree(tmp_path, truth_relations, inferred_relations):
    pedigree = SimulatedPedigree(pedigree_data_dir=tmp_path)
    nodes = sorted(
        {node for pair in truth_relations for node in pair} | {node for pair in inferred_relations for node in pair}
    )
    pedigree._final_nodes_df = pd.DataFrame({"id": nodes})
    pedigree._ground_truth_pedigree = StubRelationsPedigree(truth_relations)
    pedigree._algorithm_pedigree = StubRelationsPedigree(inferred_relations)
    pedigree._algorithm_found_pedigree = True
    return pedigree


def test_scale_error_rates_respects_scale(tmp_path):
    pedigree = SimulatedPedigree(pedigree_data_dir=tmp_path, random_seed=0)

    degree_probs, relation_probs = pedigree._scale_error_rates(scale=0.5)

    assert degree_probs["1"][0] == pytest.approx(0.995)
    assert degree_probs["1"][1] == pytest.approx(0.005)
    assert relation_probs["parent-child;child-parent"][0] == pytest.approx(0.975)
    assert relation_probs["parent-child;child-parent"][1] == pytest.approx(0.025)


def test_scale_error_rates_raises_when_probabilities_exceed_one(tmp_path):
    pedigree = SimulatedPedigree(pedigree_data_dir=tmp_path)

    with pytest.raises(ValueError):
        pedigree._scale_error_rates(scale=25)


def test_calculate_tp_fp_fn_counts():
    ground_truth = defaultdict(int, {"parent-child": 2, "siblings": 1})
    algorithm = defaultdict(int, {"parent-child": 1, "siblings": 2, "Unrelated": 1})

    tp, fp, fn = SimulatedPedigree._calculate_tp_fp_fn(ground_truth, algorithm)

    assert tp == 2
    assert fp == 2
    assert fn == 1


def test_get_metrics_reports_perfect_scores(tmp_path):
    relations = {
        ("A", "B"): {"parent-child": 1},
        ("B", "C"): {"siblings": 1},
        ("A", "C"): {"maternal grandparent-grandchild": 1},
    }
    pedigree = _build_simulated_pedigree(tmp_path, relations, relations)

    metrics = pedigree.get_metrics()

    assert metrics["Pairwise Relation Accuracy"] == pytest.approx(1.0)
    assert metrics["Relation Precision"] == pytest.approx(1.0)
    assert metrics["Relation Recall"] == pytest.approx(1.0)
    assert metrics["Relation F1"] == pytest.approx(1.0)
    assert metrics["Pairwise Degree Accuracy"] == pytest.approx(1.0)
    assert metrics["Degree Precision"] == pytest.approx(1.0)
    assert metrics["Degree Recall"] == pytest.approx(1.0)
    assert metrics["Degree F1"] == pytest.approx(1.0)
    assert metrics["Connectivity R-squared"] == pytest.approx(1.0)


def test_get_metrics_detects_relation_and_degree_errors(tmp_path):
    truth_relations = {
        ("A", "B"): {"parent-child": 1},
        ("B", "C"): {"siblings": 1},
        ("A", "C"): {"maternal grandparent-grandchild": 1},
    }
    inferred_relations = {
        ("A", "B"): {"maternal grandparent-grandchild": 1},
        ("B", "C"): {"siblings": 1},
        ("A", "C"): {"maternal grandparent-grandchild": 1},
    }
    pedigree = _build_simulated_pedigree(tmp_path, truth_relations, inferred_relations)

    metrics = pedigree.get_metrics()

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
