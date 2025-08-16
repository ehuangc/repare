import logging
import tempfile
from collections import defaultdict
from itertools import combinations

import pandas as pd
from sklearn.metrics import r2_score

from repare.pedigree import Pedigree
from repare.pedigree_reconstructor import PedigreeReconstructor

logger = logging.getLogger(__name__)


class PedigreeEvaluator:
    """
    Generates an algorithm-reconstructed pedigree and scores it against a published/ground-truth pedigree.
    """

    def __init__(self, published_relations_path: str, algorithm_nodes_path: str, algorithm_relations_path: str) -> None:
        self._published_relations_path = published_relations_path
        self._algorithm_nodes_path = algorithm_nodes_path
        self._algorithm_relations_path = algorithm_relations_path

        self._published_relation_counts: defaultdict[tuple[str, str], defaultdict[str, int]] = (
            self._load_published_relation_counts(self._published_relations_path)
        )
        self._algorithm_relation_counts: defaultdict[tuple[str, str], defaultdict[str, int]] = (
            self._load_algorithm_relation_counts(self._algorithm_nodes_path, self._algorithm_relations_path)
        )
        self._fill_uncertain_relations()

    def _load_published_relation_counts(self, path: str) -> defaultdict[tuple[str, str], defaultdict[str, int]]:
        published_relations_df = pd.read_csv(path, comment="#", dtype=str, keep_default_na=False)
        relation_counts: defaultdict[tuple[str, str], defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
        for id1, id2, degree, constraints, _ in published_relations_df.itertuples(index=False):
            if constraints:
                # Ensure only one constraint (exact relation) per node pair
                assert ";" not in constraints
                id1, id2, relation = self._sort_relation(id1, id2, constraints)
                relation_counts[(id1, id2)][relation] += 1
            else:
                id1, id2, relation = self._sort_relation(id1, id2, degree)
                relation_counts[(id1, id2)][relation] += 1
        return relation_counts

    def _load_algorithm_relation_counts(
        self, nodes_path: str, relations_path: str
    ) -> defaultdict[tuple[str, str], defaultdict[str, int]]:
        self.algorithm_pedigree: Pedigree = self._run_algorithm(nodes_path, relations_path)
        algorithm_relations: defaultdict[tuple[str, str], defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))

        for id1, id2 in combinations(self.algorithm_pedigree.node_to_data, 2):
            if not id1.isnumeric() and not id2.isnumeric():  # Skip placeholder nodes
                relations_between_nodes = self.algorithm_pedigree.get_relations_between_nodes(
                    id1, id2, include_maternal_paternal=True
                )
                for relation, count in relations_between_nodes.items():
                    id1, id2, relation = self._sort_relation(id1, id2, relation)
                    algorithm_relations[(id1, id2)][relation] += count
        return algorithm_relations

    @staticmethod
    def _sort_relation(id1: str, id2: str, relation: str) -> tuple[str, str, str]:
        flipped_relations = {
            "parent-child": "child-parent",
            "child-parent": "parent-child",
            "siblings": "siblings",  # Symmetric
            "maternal aunt/uncle-nephew/niece": "maternal nephew/niece-aunt/uncle",
            "maternal nephew/niece-aunt/uncle": "maternal aunt/uncle-nephew/niece",
            "paternal aunt/uncle-nephew/niece": "paternal nephew/niece-aunt/uncle",
            "paternal nephew/niece-aunt/uncle": "paternal aunt/uncle-nephew/niece",
            "maternal grandparent-grandchild": "maternal grandchild-grandparent",
            "maternal grandchild-grandparent": "maternal grandparent-grandchild",
            "paternal grandparent-grandchild": "paternal grandchild-grandparent",
            "paternal grandchild-grandparent": "paternal grandparent-grandchild",
            "maternal half-siblings": "maternal half-siblings",  # Symmetric
            "paternal half-siblings": "paternal half-siblings",  # Symmetric
            "1": "1",  # Symmetric
            "2": "2",  # Symmetric
        }
        if id2 < id1:
            return id2, id1, flipped_relations[relation]
        else:
            return id1, id2, relation

    @staticmethod
    def _run_algorithm(nodes_path: str, relations_path: str) -> Pedigree:
        with tempfile.TemporaryDirectory() as temp_dir:
            pedigree_reconstructor = PedigreeReconstructor(
                relations_path, nodes_path, outputs_dir=temp_dir, max_candidate_pedigrees=1000, plot=False
            )
            return pedigree_reconstructor.find_best_pedigree()

    def _fill_uncertain_relations(self) -> None:
        uncertain_to_exact_relations = {
            "1": ["parent-child", "child-parent", "siblings"],
            "2": [
                "maternal aunt/uncle-nephew/niece",
                "maternal nephew/niece-aunt/uncle",
                "paternal aunt/uncle-nephew/niece",
                "paternal nephew/niece-aunt/uncle",
                "maternal grandparent-grandchild",
                "maternal grandchild-grandparent",
                "paternal grandparent-grandchild",
                "paternal grandchild-grandparent",
                "maternal half-siblings",
                "paternal half-siblings",
            ],
        }

        for (id1, id2), relation_counts_between_nodes in self._published_relation_counts.items():
            for uncertain_relation, count in list(relation_counts_between_nodes.items()):  # Cast to list to copy items
                if uncertain_relation not in uncertain_to_exact_relations:
                    continue

                for exact_relation in uncertain_to_exact_relations[uncertain_relation]:
                    available_count = self._algorithm_relation_counts[(id1, id2)][exact_relation]
                    assign_count = min(count, available_count)
                    self._published_relation_counts[(id1, id2)][exact_relation] += assign_count
                    self._published_relation_counts[(id1, id2)][uncertain_relation] -= assign_count

                    count -= assign_count
                    if count == 0:
                        del relation_counts_between_nodes[uncertain_relation]
                        break

    def get_metrics(self) -> dict[str, float]:
        metrics: dict[str, float] = dict()
        pairwise_relation_accuracy, relation_precision, relation_recall, relation_f1 = (
            self._calculate_relation_metrics()
        )
        pairwise_degree_accuracy, degree_precision, degree_recall, degree_f1 = self._calculate_degree_metrics()

        metrics["Pairwise Relation Accuracy"] = pairwise_relation_accuracy
        metrics["Relation Precision"] = relation_precision
        metrics["Relation Recall"] = relation_recall
        metrics["Relation F1"] = relation_f1
        metrics["Pairwise Degree Accuracy"] = pairwise_degree_accuracy
        metrics["Degree Precision"] = degree_precision
        metrics["Degree Recall"] = degree_recall
        metrics["Degree F1"] = degree_f1
        metrics["Connectivity R-squared"] = self._calculate_connectivity_r_squared()
        metrics["Kinship Inference Errors"] = self._calculate_kinship_inference_errors()
        return metrics

    @staticmethod
    def _calculate_tp_fp_fn(
        published_counts: defaultdict[str, int], algorithm_counts: defaultdict[str, int], nodes: tuple[str, str]
    ) -> tuple[int, int, int]:
        tp = 0  # True positives
        fp = 0  # False positives
        fn = 0  # False negatives
        relations = published_counts.keys() | algorithm_counts.keys()
        for relation in relations:
            true_count = published_counts[relation]
            algorithm_count = algorithm_counts[relation]

            if true_count == algorithm_count:
                tp += true_count
            elif true_count > algorithm_count:
                tp += algorithm_count
                fn += true_count - algorithm_count
                logger.info(f"False Negative: {nodes[0]} - {nodes[1]}: {relation} ({true_count} > {algorithm_count})")
            else:
                tp += true_count
                fp += algorithm_count - true_count
                logger.info(f"False Positive: {nodes[0]} - {nodes[1]}: {relation} ({true_count} < {algorithm_count})")
        return tp, fp, fn

    def _calculate_relation_metrics(self) -> tuple[float, float, float, float]:
        correct_node_pairs: int = 0
        total_node_pairs: int = 0
        relation_tp: int = 0
        relation_fp: int = 0
        relation_fn: int = 0

        nodes = [node for node in self.algorithm_pedigree.node_to_data if not node.isnumeric()]
        for id1, id2 in combinations(sorted(nodes), 2):
            published_relations_between_nodes = self._published_relation_counts[(id1, id2)]
            algorithm_relations_between_nodes = self._algorithm_relation_counts[(id1, id2)]

            if published_relations_between_nodes == algorithm_relations_between_nodes:
                correct_node_pairs += 1
            total_node_pairs += 1

            tp, fp, fn = self._calculate_tp_fp_fn(
                published_relations_between_nodes, algorithm_relations_between_nodes, (id1, id2)
            )
            relation_tp += tp
            relation_fp += fp
            relation_fn += fn

        pairwise_relation_accuracy = correct_node_pairs / total_node_pairs
        relation_precision = relation_tp / (relation_tp + relation_fp)
        relation_recall = relation_tp / (relation_tp + relation_fn)
        relation_f1 = 2 * (relation_precision * relation_recall) / (relation_precision + relation_recall)
        relation_f1 = (
            (2 * relation_precision * relation_recall) / (relation_precision + relation_recall)
            if relation_precision + relation_recall > 0
            else 0
        )
        return pairwise_relation_accuracy, relation_precision, relation_recall, relation_f1

    def _calculate_degree_metrics(self) -> tuple[float, float, float, float]:
        correct_node_pairs: int = 0
        total_node_pairs: int = 0
        degree_tp: int = 0
        degree_fp: int = 0
        degree_fn: int = 0

        nodes = [node for node in self.algorithm_pedigree.node_to_data if not node.isnumeric()]
        for id1, id2 in combinations(sorted(nodes), 2):
            published_relations_between_nodes = self._published_relation_counts[(id1, id2)]
            algorithm_relations_between_nodes = self._algorithm_relation_counts[(id1, id2)]

            published_degrees_between_nodes = defaultdict(int)
            algorithm_degrees_between_nodes = defaultdict(int)
            for relation in ["parent-child", "child-parent", "siblings"]:
                published_degrees_between_nodes["1"] += published_relations_between_nodes[relation]
                algorithm_degrees_between_nodes["1"] += algorithm_relations_between_nodes[relation]

            for relation in [
                "maternal aunt/uncle-nephew/niece",
                "paternal aunt/uncle-nephew/niece",
                "maternal nephew/niece-aunt/uncle",
                "paternal nephew/niece-aunt/uncle",
                "maternal grandparent-grandchild",
                "paternal grandparent-grandchild",
                "maternal grandchild-grandparent",
                "paternal grandchild-grandparent",
                "maternal half-siblings",
                "paternal half-siblings",
            ]:
                published_degrees_between_nodes["2"] += published_relations_between_nodes[relation]
                algorithm_degrees_between_nodes["2"] += algorithm_relations_between_nodes[relation]

            if published_degrees_between_nodes == algorithm_degrees_between_nodes:
                correct_node_pairs += 1
            total_node_pairs += 1

            tp, fp, fn = self._calculate_tp_fp_fn(
                published_degrees_between_nodes, algorithm_degrees_between_nodes, (id1, id2)
            )
            degree_tp += tp
            degree_fp += fp
            degree_fn += fn

        pairwise_degree_accuracy = correct_node_pairs / total_node_pairs
        degree_precision = degree_tp / (degree_tp + degree_fp)
        degree_recall = degree_tp / (degree_tp + degree_fn)
        degree_f1 = (
            (2 * degree_precision * degree_recall) / (degree_precision + degree_recall)
            if degree_precision + degree_recall > 0
            else 0
        )
        return pairwise_degree_accuracy, degree_precision, degree_recall, degree_f1

    def _calculate_connectivity_r_squared(self) -> float:
        published_relation_counter: defaultdict[str, int] = defaultdict(int)
        algorithm_relation_counter: defaultdict[str, int] = defaultdict(int)

        nodes = [node for node in self.algorithm_pedigree.node_to_data if not node.isnumeric()]
        for node1, node2 in combinations(sorted(nodes), 2):
            published_relations: defaultdict[str, int] = self._published_relation_counts[(node1, node2)]
            algorithm_relations: defaultdict[str, int] = self._algorithm_relation_counts[(node1, node2)]

            for _relation, count in published_relations.items():
                published_relation_counter[node1] += count
                published_relation_counter[node2] += count

            for _relation, count in algorithm_relations.items():
                algorithm_relation_counter[node1] += count
                algorithm_relation_counter[node2] += count

        published_connectivities: list[int] = []
        algorithm_connectivities: list[int] = []
        for node in nodes:
            published_connectivities.append(published_relation_counter[node])
            algorithm_connectivities.append(algorithm_relation_counter[node])
        return r2_score(published_connectivities, algorithm_connectivities)

    def _calculate_kinship_inference_errors(self) -> int:
        """
        Calculate the number of node pairs that share a different inferred kinship degree than in the published pedigree
        or share a relation constraint not consistent with the published pedigree.
        """
        published_exact_relations = pd.read_csv(
            self._published_relations_path, dtype=str, comment="#", keep_default_na=False
        )
        inferred_relations = pd.read_csv(self._algorithm_relations_path, dtype=str, comment="#", keep_default_na=False)

        pair_to_published_degree = {}
        for id1, id2, degree, _, _ in published_exact_relations.itertuples(index=False):
            assert degree in ["1", "2"]
            pair_to_published_degree[tuple(sorted((id1, id2)))] = degree

        # Map constraints to their flipped value
        flipped_constraints = {
            "parent-child": "child-parent",
            "child-parent": "parent-child",
            "maternal aunt/uncle-nephew/niece": "maternal nephew/niece-aunt/uncle",
            "paternal aunt/uncle-nephew/niece": "paternal nephew/niece-aunt/uncle",
            "maternal nephew/niece-aunt/uncle": "maternal aunt/uncle-nephew/niece",
            "paternal nephew/niece-aunt/uncle": "paternal aunt/uncle-nephew/niece",
            "maternal grandparent-grandchild": "maternal grandchild-grandparent",
            "paternal grandparent-grandchild": "paternal grandchild-grandparent",
            "maternal grandchild-grandparent": "maternal grandparent-grandchild",
            "paternal grandchild-grandparent": "paternal grandparent-grandchild",
            "siblings": "siblings",  # Symmetric
            "maternal half-siblings": "maternal half-siblings",  # Symmetric
            "paternal half-siblings": "paternal half-siblings",  # Symmetric
        }

        pair_to_inferred_degree = {}
        pair_to_inferred_constraints = {}
        for id1, id2, degree, constraints in inferred_relations.itertuples(index=False):
            if degree == "1" or degree == "2":
                pair_to_inferred_degree[tuple(sorted((id1, id2)))] = degree
            if constraints:
                if id2 < id1:
                    assert sorted((id1, id2)) == [id2, id1]
                    # Split constraints and map each to its flipped value
                    constraints_list = [c.strip() for c in constraints.split(";")]
                    flipped = [flipped_constraints[c] for c in constraints_list]
                    pair_to_inferred_constraints[tuple(sorted((id1, id2)))] = set(flipped)
                else:
                    assert sorted((id1, id2)) == [id1, id2]
                    pair_to_inferred_constraints[tuple(sorted((id1, id2)))] = set(constraints.split(";"))

        # Compare the degree dicts
        kinship_inference_errors = 0
        for pair, algorithm_degree in pair_to_inferred_degree.items():
            if pair not in pair_to_published_degree:
                kinship_inference_errors += 1
                continue
            published_degree = pair_to_published_degree[pair]
            if algorithm_degree != published_degree:
                kinship_inference_errors += 1
                continue

        for pair in pair_to_published_degree:
            if pair not in pair_to_inferred_degree:
                kinship_inference_errors += 1

        # Count within-degree relation constraint inference errors
        for id1, id2, _, constraints, _ in published_exact_relations.itertuples(index=False):
            # Skip "dotted lines"
            if not constraints:
                continue
            pair = tuple(sorted((id1, id2)))
            if id2 < id1:
                constraints = flipped_constraints[constraints]
            if pair in pair_to_inferred_constraints and constraints not in pair_to_inferred_constraints[pair]:
                kinship_inference_errors += 1
        return kinship_inference_errors

    def write_relation_differences(self, path: str) -> None:
        """
        Write the differences between the published and inferred relations to a CSV file.
        """
        false_positives: defaultdict[tuple[str, str], list[str]] = defaultdict(list)
        false_negatives: defaultdict[tuple[str, str], list[str]] = defaultdict(list)
        nodes = [node for node in self.algorithm_pedigree.node_to_data if not node.isnumeric()]
        for id1, id2 in combinations(sorted(nodes), 2):
            assert id1 < id2
            published_relations_between_nodes = self._published_relation_counts[(id1, id2)]
            algorithm_relations_between_nodes = self._algorithm_relation_counts[(id1, id2)]

            relations = published_relations_between_nodes.keys() | algorithm_relations_between_nodes.keys()
            for relation in relations:
                true_count = published_relations_between_nodes[relation]
                algorithm_count = algorithm_relations_between_nodes[relation]

                if true_count > algorithm_count:
                    for _ in range(true_count - algorithm_count):
                        false_negatives[(id1, id2)].append(relation)

                elif algorithm_count > true_count:
                    for _ in range(algorithm_count - true_count):
                        false_positives[(id1, id2)].append(relation)

        # Write false positives and false negatives to CSV file
        with open(path, "w") as file:
            file.write("id1,id2,published_relation,inferred_relation\n")
            for id1, id2 in sorted(set(false_positives.keys()) | set(false_negatives.keys())):
                false_positive_relations = (
                    ";".join(false_positives[(id1, id2)]) if (id1, id2) in false_positives else "None"
                )
                false_negative_relations = (
                    ";".join(false_negatives[(id1, id2)]) if (id1, id2) in false_negatives else "None"
                )
                file.write(f"{id1},{id2},{false_positive_relations},{false_negative_relations}\n")
