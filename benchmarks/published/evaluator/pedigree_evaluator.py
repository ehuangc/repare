import logging
import tempfile
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import pandas as pd
from sklearn.metrics import r2_score

from repare.pedigree import Pedigree
from repare.pedigree_reconstructor import PedigreeReconstructor

logger = logging.getLogger(__name__)


class PedigreeEvaluator:
    """
    Generates an algorithm-reconstructed pedigree and scores it against a published/ground-truth pedigree.
    """

    def __init__(
        self,
        published_relations_path: str,
        algorithm_nodes_path: str,
        algorithm_relations_path: str,
        outputs_dir: Path | str,
    ) -> None:
        self._published_relations_path = published_relations_path
        self._algorithm_nodes_path = algorithm_nodes_path
        self._algorithm_relations_path = algorithm_relations_path
        self._outputs_dir = Path(outputs_dir)

        self._published_relation_counts: defaultdict[tuple[str, str], defaultdict[str, int]] = (
            self._load_published_relation_counts(self._published_relations_path)
        )
        self._algorithm_relation_counts: defaultdict[tuple[str, str], defaultdict[str, int]] = (
            self._load_algorithm_relation_counts(self._algorithm_nodes_path, self._algorithm_relations_path)
        )
        (
            self._input_pair_to_constraints,
            self._input_pair_to_relations,
        ) = self._load_input_relation_mappings(self._algorithm_relations_path)
        self._published_pedigree: Pedigree | None = None
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
        self.algorithm_pedigree: Pedigree = self._run_algorithm(
            nodes_path, relations_path, outputs_dir=self._outputs_dir
        )
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

    def _load_input_relation_mappings(
        self, relations_path: str
    ) -> tuple[
        defaultdict[tuple[str, str], list[tuple[str, ...]]],
        defaultdict[tuple[str, str], list[tuple[str, str, bool]]],
    ]:
        relation_data = pd.read_csv(relations_path, dtype=str, comment="#", keep_default_na=False)
        if "force_constraints" not in relation_data.columns:
            relation_data["force_constraints"] = ""
        relation_data["force_constraints"] = relation_data["force_constraints"].eq("True")

        flipped_constraints = {
            "parent-child": "child-parent",
            "child-parent": "parent-child",
            "siblings": "siblings",
            "maternal aunt/uncle-nephew/niece": "maternal nephew/niece-aunt/uncle",
            "maternal nephew/niece-aunt/uncle": "maternal aunt/uncle-nephew/niece",
            "paternal aunt/uncle-nephew/niece": "paternal nephew/niece-aunt/uncle",
            "paternal nephew/niece-aunt/uncle": "paternal aunt/uncle-nephew/niece",
            "maternal grandparent-grandchild": "maternal grandchild-grandparent",
            "maternal grandchild-grandparent": "maternal grandparent-grandchild",
            "paternal grandparent-grandchild": "paternal grandchild-grandparent",
            "paternal grandchild-grandparent": "paternal grandparent-grandchild",
            "maternal half-siblings": "maternal half-siblings",
            "paternal half-siblings": "paternal half-siblings",
            "double cousins": "double cousins",
            "half aunt/uncle-half nephew/niece": "half nephew/niece-half aunt/uncle",
            "half nephew/niece-half aunt/uncle": "half aunt/uncle-half nephew/niece",
            "greatgrandparent-greatgrandchild": "greatgrandchild-greatgrandparent",
            "greatgrandchild-greatgrandparent": "greatgrandparent-greatgrandchild",
            "grandaunt/granduncle-grandnephew/grandniece": "grandnephew/grandniece-grandaunt/granduncle",
            "grandnephew/grandniece-grandaunt/granduncle": "grandaunt/granduncle-grandnephew/grandniece",
            "first cousins": "first cousins",
        }

        def sort_nodes(row: pd.Series) -> pd.Series:
            if row["id2"] < row["id1"]:
                constraints = row["constraints"]
                if constraints:
                    constraints_list = [c.strip() for c in constraints.split(";")]
                    flipped = [flipped_constraints[c] for c in constraints_list]
                    constraints = ";".join(flipped)
                return pd.Series(
                    {
                        "id1": row["id2"],
                        "id2": row["id1"],
                        "degree": row["degree"],
                        "constraints": constraints,
                        "force_constraints": row["force_constraints"],
                    }
                )
            return row

        relation_data = relation_data.apply(sort_nodes, axis=1)

        default_constraints = {
            "1": "parent-child;child-parent;siblings",
            "2": (
                "maternal aunt/uncle-nephew/niece;"
                "maternal nephew/niece-aunt/uncle;"
                "paternal aunt/uncle-nephew/niece;"
                "paternal nephew/niece-aunt/uncle;"
                "maternal grandparent-grandchild;"
                "maternal grandchild-grandparent;"
                "paternal grandparent-grandchild;"
                "paternal grandchild-grandparent;"
                "maternal half-siblings;"
                "paternal half-siblings;"
                "double cousins"
            ),
            "3": (
                "half aunt/uncle-half nephew/niece;"
                "half nephew/niece-half aunt/uncle;"
                "greatgrandparent-greatgrandchild;"
                "greatgrandchild-greatgrandparent;"
                "grandaunt/granduncle-grandnephew/grandniece;"
                "grandnephew/grandniece-grandaunt/granduncle;"
                "first cousins"
            ),
        }
        relation_data["constraints"] = relation_data.apply(
            lambda row: row["constraints"] if row["constraints"] else default_constraints[row["degree"]], axis=1
        )

        pair_to_constraints: defaultdict[tuple[str, str], list[tuple[str, ...]]] = defaultdict(list)
        pair_to_relations: defaultdict[tuple[str, str], list[tuple[str, str, bool]]] = defaultdict(list)
        for id1, id2, degree, constraints, force_constraints in relation_data.itertuples(index=False):
            constraint_tuple = tuple(c.strip() for c in constraints.split(";"))
            pair_to_constraints[(id1, id2)].append(constraint_tuple)
            # Only 1st- and 2nd-degree relations are included in pair_to_relations(_so_far)
            if degree in ["1", "2"]:
                pair_to_relations[(id1, id2)].append((degree, constraints, bool(force_constraints)))

        for constraint_lists in pair_to_constraints.values():
            constraint_lists.sort(key=len)
        return pair_to_constraints, pair_to_relations

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
            "double cousins": "double cousins",  # Symmetric
            "1": "1",  # Symmetric
            "2": "2",  # Symmetric
        }
        if id2 < id1:
            return id2, id1, flipped_relations[relation]
        else:
            return id1, id2, relation

    @staticmethod
    def _run_algorithm(nodes_path: str, relations_path: str, outputs_dir: Path | str) -> Pedigree:
        outputs_path = Path(outputs_dir)
        outputs_path.mkdir(parents=True, exist_ok=True)
        pedigree_reconstructor = PedigreeReconstructor(
            relations_path, nodes_path, outputs_dir=outputs_path, max_candidate_pedigrees=1000, plot=False
        )
        return pedigree_reconstructor.find_best_pedigree()

    def _get_published_pedigree(self) -> Pedigree:
        if self._published_pedigree is None:
            with tempfile.TemporaryDirectory() as temp_dir:
                self._published_pedigree = self._run_algorithm(
                    self._algorithm_nodes_path, self._published_relations_path, outputs_dir=temp_dir
                )
        return self._published_pedigree

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
                "double cousins",
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
        metrics.update(self.count_input_relation_inconsistencies())
        return metrics

    def count_input_relation_inconsistencies(self) -> dict[str, int]:
        """
        Count how many input relations are inconsistent with the published pedigree and with the inferred pedigree.
        """
        published_pedigree = self._get_published_pedigree()
        published_inconsistencies, _ = published_pedigree.count_inconsistencies(
            self._input_pair_to_constraints, self._input_pair_to_relations, check_half_siblings=True
        )
        inferred_inconsistencies, _ = self.algorithm_pedigree.count_inconsistencies(
            self._input_pair_to_constraints, self._input_pair_to_relations, check_half_siblings=True
        )
        return {
            "Published Pedigree Input Inconsistencies": published_inconsistencies,
            "Inferred Pedigree Input Inconsistencies": inferred_inconsistencies,
        }

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
                "double cousins",
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

        inconsistencies = self.count_input_relation_inconsistencies()
        published_inconsistencies = inconsistencies["Published Pedigree Input Inconsistencies"]
        inferred_inconsistencies = inconsistencies["Inferred Pedigree Input Inconsistencies"]

        # Write false positives and false negatives to CSV file
        with open(path, "w") as file:
            file.write(
                f"# Published pedigree inconsistencies: {published_inconsistencies}, "
                f"inferred pedigree inconsistencies: {inferred_inconsistencies}\n"
            )
            file.write("id1,id2,published_relation,inferred_relation\n")
            for id1, id2 in sorted(set(false_positives.keys()) | set(false_negatives.keys())):
                # False negatives are published relations not inferred by the algorithm
                false_negative_relations = (
                    ";".join(false_negatives[(id1, id2)]) if (id1, id2) in false_negatives else "None"
                )
                # False positives are inferred relations not in the published pedigree
                false_positive_relations = (
                    ";".join(false_positives[(id1, id2)]) if (id1, id2) in false_positives else "None"
                )
                file.write(f"{id1},{id2},{false_negative_relations},{false_positive_relations}\n")
