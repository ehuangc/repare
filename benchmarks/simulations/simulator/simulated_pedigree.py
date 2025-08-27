import math
import os
import random
import tempfile
from collections import defaultdict
from itertools import combinations

import pandas as pd
from sklearn.metrics import r2_score

from repare.pedigree import Pedigree
from repare.pedigree_reconstructor import PedigreeReconstructor


class SimulatedPedigree:
    """
    Simulates a pedigree, masks and corrupts the data to generate algorithm inputs,
    runs the pedigree reconstruction algorithm, and calculates performance metrics.
    """

    def __init__(
        self,
        p_mask_node: float = 0.4,
        error_rate_scale: float = 1,
        max_candidate_pedigrees: int = 1000,
        epsilon: float | None = None,
        random_seed: int | None = None,
    ) -> None:
        self._ground_truth_pedigree = Pedigree()
        self._y_haplogroup_pool = ["a", "b"]
        self._mt_haplogroup_pool = ["a", "b", "c", "d", "e"]
        # "True" value, will be more conservative when writing node data
        self._p_can_have_children = 0.6
        self._mean_children_per_mate = 2
        self._sd_children_per_mate = 1
        self._num_generations = 4
        self._generation_zero_starting_size = 3

        # Probability that a node will be masked (i.e., not included in node data)
        self._p_mask_node = p_mask_node
        # Scale to apply to relation classification error rates
        self._error_rate_scale = error_rate_scale

        if p_mask_node < 0 or p_mask_node > 1:
            raise ValueError("p_mask_node must be between 0 and 1.")
        if error_rate_scale < 0:
            raise ValueError("error_rate_scale must be non-negative.")

        self._node_count = 0
        # Maps generation number to set of node IDs
        self._generation_to_nodes = defaultdict(set)
        # Maps node ID to generation number
        self._node_to_generation = {}
        self._max_candidate_pedigrees = max_candidate_pedigrees
        self._epsilon = epsilon
        self._random_seed = random_seed
        random.seed(self._random_seed)

        self._base_degree_classification_probs: dict[str, tuple[float, float, float, float]] = {
            "1": (0.99, 0.01, 0.0, 0.0),
            "2": (0.01, 0.94, 0.05, 0.0),
            "3": (0.0, 0.02, 0.88, 0.1),
            "Unrelated": (0.0, 0.0, 0.01, 0.99),
        }
        self._base_relation_classification_probs: dict[str, tuple[float, float]] = {
            "parent-child;child-parent": (0.95, 0.05),
            "siblings": (0.05, 0.95),
        }

    def create_pedigree(self) -> None:
        self._create_node(generation=0, sex="M", can_have_children=True)
        self._create_node(generation=0, sex="F", can_have_children=True)
        for _ in range(self._generation_zero_starting_size - 2):
            self._create_node(generation=0)

        for generation in range(self._num_generations - 1):
            # 1) Make copy so iterator doesn't change size, 2) Sort to ensure deterministic results
            for node in sorted(self._generation_to_nodes[generation]):
                if self._ground_truth_pedigree.node_to_data[node]["can_have_children"]:
                    self._create_children(node)

    def _create_node(
        self,
        generation: int,
        sex: str | None = None,
        y_haplogroup: str | None = None,
        mt_haplogroup: str | None = None,
        can_have_children: bool | None = None,
    ) -> str:
        """
        Create a new node in the pedigree with the given attributes. Returns the new node's ID.
        """
        node_id = "N" + str(self._node_count)
        self._node_count += 1

        if sex is None:
            sex = random.choice(["M", "F"])
        if y_haplogroup is None and sex == "M":
            y_haplogroup = random.choice(self._y_haplogroup_pool)
        if mt_haplogroup is None:
            mt_haplogroup = random.choice(self._mt_haplogroup_pool)
        if can_have_children is None:
            can_have_children = random.random() < self._p_can_have_children

        # Set can_be_inbred to True; we will update later in self._get_nodes()
        self._ground_truth_pedigree.add_node(
            node_id,
            sex,
            y_haplogroup,
            mt_haplogroup,
            can_have_children,
            can_be_inbred=True,
            years_before_present=math.nan,
        )
        self._generation_to_nodes[generation].add(node_id)
        self._node_to_generation[node_id] = generation
        return node_id

    def _create_children(self, parent1: str) -> None:
        """
        Given a parent node, select a valid mate and produce a random number of children.
        """
        parent1_sex = self._ground_truth_pedigree.get_data(parent1)["sex"]
        parent1_generation = self._node_to_generation[parent1]

        # Parent 2 is either an existing node or is an "outside" node
        parent2 = None
        parent2_sex = "F" if parent1_sex == "M" else "M"
        if random.random() < 0.25:
            # Existing mate
            potential_mates = []
            for node in self._ground_truth_pedigree.node_to_data:
                if (
                    node != parent1
                    and self._ground_truth_pedigree.get_data(node)["can_have_children"]
                    and 0 <= parent1_generation - self._node_to_generation[node] <= 1
                ):
                    if parent1_sex == "M" and self._ground_truth_pedigree.get_data(node)["sex"] == "F":
                        potential_mates.append(node)
                    if parent1_sex == "F" and self._ground_truth_pedigree.get_data(node)["sex"] == "M":
                        potential_mates.append(node)
            if len(potential_mates) > 0:
                parent2 = random.choice(potential_mates)
            else:
                parent2 = self._create_node(generation=parent1_generation, sex=parent2_sex, can_have_children=True)
        else:
            # Outside mate
            parent2 = self._create_node(generation=parent1_generation, sex=parent2_sex, can_have_children=True)

        child_generation = max(parent1_generation, self._node_to_generation[parent2]) + 1
        num_children = round(random.normalvariate(mu=self._mean_children_per_mate, sigma=self._sd_children_per_mate))
        num_children = max(1, num_children)
        for _ in range(num_children):
            child_sex = random.choice(["M", "F"])
            child_y_haplogroup = None
            child_mt_haplogroup = None
            if parent1_sex == "M":
                if child_sex == "M":
                    child_y_haplogroup = self._ground_truth_pedigree.get_data(parent1)["y_haplogroup"]
                child_mt_haplogroup = self._ground_truth_pedigree.get_data(parent2)["mt_haplogroup"]
            else:
                if child_sex == "M":
                    child_y_haplogroup = self._ground_truth_pedigree.get_data(parent2)["y_haplogroup"]
                child_mt_haplogroup = self._ground_truth_pedigree.get_data(parent1)["mt_haplogroup"]

            child = self._create_node(
                generation=child_generation,
                sex=child_sex,
                y_haplogroup=child_y_haplogroup,
                mt_haplogroup=child_mt_haplogroup,
            )
            self._ground_truth_pedigree.add_parent_relation(parent1, child)
            self._ground_truth_pedigree.add_parent_relation(parent2, child)

    def _get_nodes(self) -> pd.DataFrame:
        nodes_list: list[
            str, str, str, str, str, str
        ] = []  # id, sex, y_haplogroup, mt_haplogroup, can_have_children, can_be_inbred
        pedigree_relation_pairs: set[tuple[str, str]] = self._ground_truth_pedigree.get_related_pairs()
        for node in self._ground_truth_pedigree.node_to_data:
            sex = self._ground_truth_pedigree.get_data(node)["sex"]
            y_haplogroup = self._ground_truth_pedigree.get_data(node)["y_haplogroup"]
            mt_haplogroup = self._ground_truth_pedigree.get_data(node)["mt_haplogroup"]

            can_have_children = "True"
            # Even if node has no children, conservatively set can_have_children to False
            if not self._ground_truth_pedigree.get_data(node)["can_have_children"] and random.random() < 0.25:
                can_have_children = "False"

            can_be_inbred = "True"
            father = self._ground_truth_pedigree.get_father(node)
            mother = self._ground_truth_pedigree.get_mother(node)
            # Even if node is not inbred, conservatively set can_be_inbred to False
            if (
                (father, mother) not in pedigree_relation_pairs
                and (mother, father) not in pedigree_relation_pairs
                and random.random() < 0.5
            ):
                can_be_inbred = "False"

            # Conservatively do not set sample age
            years_before_present = math.nan
            nodes_list.append(
                (node, sex, y_haplogroup, mt_haplogroup, can_have_children, can_be_inbred, years_before_present)
            )
        return pd.DataFrame(
            nodes_list,
            columns=[
                "id",
                "sex",
                "y_haplogroup",
                "mt_haplogroup",
                "can_have_children",
                "can_be_inbred",
                "years_before_present",
            ],
        )

    def _get_first_degree_relations(self) -> pd.DataFrame:
        relations_list: list[str, str, str, str] = []  # id1, id2, degree, constraints
        for parent, child in self._ground_truth_pedigree.get_parent_child_pairs(include_placeholders=False):
            relations_list.append((parent, child, "1", "parent-child;child-parent"))

        for sibling1, sibling2 in self._ground_truth_pedigree.get_sibling_pairs(include_placeholders=False):
            relations_list.append((sibling1, sibling2, "1", "siblings"))
        return pd.DataFrame(relations_list, columns=["id1", "id2", "degree", "constraints"])

    def _get_second_degree_relations(self) -> pd.DataFrame:
        relations_list: list[str, str, str, str] = []  # id1, id2, degree, constraints (no constraints for 2nd-degree)
        for aunt_uncle, nephew_niece in self._ground_truth_pedigree.get_aunt_uncle_nephew_niece_pairs(
            include_placeholders=False
        ):
            relations_list.append((aunt_uncle, nephew_niece, "2", ""))

        for grandparent, grandchild in self._ground_truth_pedigree.get_grandparent_grandchild_pairs(
            include_placeholders=False
        ):
            relations_list.append((grandparent, grandchild, "2", ""))

        for half_sibling1, half_sibling2 in self._ground_truth_pedigree.get_half_sibling_pairs(
            include_placeholders=False
        ):
            relations_list.append((half_sibling1, half_sibling2, "2", ""))
        return pd.DataFrame(relations_list, columns=["id1", "id2", "degree", "constraints"])

    def _get_third_degree_relations(self) -> pd.DataFrame:
        relations_list: list[str, str, str, str] = []  # id1, id2, degree, constraints (no constraints for 3rd-degree)
        for half_aunt_uncle, half_nephew_niece in self._ground_truth_pedigree.get_half_aunt_uncle_nephew_niece_pairs(
            include_placeholders=False
        ):
            relations_list.append((half_aunt_uncle, half_nephew_niece, "3", ""))

        for greatgrandparent, greatgrandchild in self._ground_truth_pedigree.get_greatgrandparent_greatgrandchild_pairs(
            include_placeholders=False
        ):
            relations_list.append((greatgrandparent, greatgrandchild, "3", ""))

        for (
            grandaunt_granduncle,
            grandnephew_grandniece,
        ) in self._ground_truth_pedigree.get_grandaunt_granduncle_grandnephew_grandniece_pairs(
            include_placeholders=False
        ):
            relations_list.append((grandaunt_granduncle, grandnephew_grandniece, "3", ""))

        for first_cousin1, first_cousin2 in self._ground_truth_pedigree.get_first_cousin_pairs(
            include_placeholders=False
        ):
            relations_list.append((first_cousin1, first_cousin2, "3", ""))
        return pd.DataFrame(relations_list, columns=["id1", "id2", "degree", "constraints"])

    def _get_unrelated_relations(self) -> pd.DataFrame:
        relations_list: list[str, str, str, str] = []  # id1, id2, degree, constraints (no constraints for unrelated)
        for node1, node2 in combinations(self._ground_truth_pedigree.node_to_data.keys(), 2):
            if not self._ground_truth_pedigree.get_relations_between_nodes(node1, node2):
                relations_list.append((node1, node2, "Unrelated", ""))
        return pd.DataFrame(relations_list, columns=["id1", "id2", "degree", "constraints"])

    def _scale_error_rates(
        self, scale: float
    ) -> tuple[dict[str, tuple[float, float, float, float]], dict[str, tuple[float, float]]]:
        scaled_degree_probs = {}
        for degree, probs in self._base_degree_classification_probs.items():
            correct_prob_idx = None
            if degree == "1":
                correct_prob_idx = 0
            elif degree == "2":
                correct_prob_idx = 1
            elif degree == "3":
                correct_prob_idx = 2
            else:
                correct_prob_idx = 3

            scaled_probs = []
            for idx, prob in enumerate(probs):
                if idx == correct_prob_idx:
                    scaled_probs.append(1 - ((1 - prob) * scale))
                else:
                    scaled_probs.append(prob * scale)
            scaled_degree_probs[degree] = tuple(scaled_probs)

        scaled_relation_probs = {}
        for relation, probs in self._base_relation_classification_probs.items():
            correct_prob_idx = 0 if relation == "parent-child;child-parent" else 1
            scaled_probs = []
            for idx, prob in enumerate(probs):
                if idx == correct_prob_idx:
                    scaled_probs.append(1 - ((1 - prob) * scale))
                else:
                    scaled_probs.append(prob * scale)
            scaled_relation_probs[relation] = tuple(scaled_probs)

        for _degree, probs in scaled_degree_probs.items():
            if any(prob > 1 for prob in probs):
                raise ValueError("Scale is too high. Error rates exceed 1.")
        for _relation, probs in scaled_relation_probs.items():
            if any(prob > 1 for prob in probs):
                raise ValueError("Scale is too high. Error rates exceed 1.")
        return scaled_degree_probs, scaled_relation_probs

    def mask_and_corrupt_data(self) -> None:
        nodes_df = self._get_nodes()
        relations_df = pd.concat(
            [
                self._get_first_degree_relations(),
                self._get_second_degree_relations(),
                self._get_third_degree_relations(),
                self._get_unrelated_relations(),
            ]
        )
        self._ground_truth_pedigree.clean_data()
        self._ground_truth_nodes_df = nodes_df.copy()
        # Remove unrelated entries
        self._ground_truth_relations_df = relations_df[relations_df["degree"] != "Unrelated"].copy()

        nodes_to_mask = [node for node in nodes_df["id"] if random.random() < self._p_mask_node]
        nodes_df = nodes_df.loc[~nodes_df["id"].isin(nodes_to_mask)]
        relations_df = relations_df.loc[
            ~relations_df["id1"].isin(nodes_to_mask) & ~relations_df["id2"].isin(nodes_to_mask)
        ]

        # Ensure only the first (closest) relation is kept for each pair of nodes to simulate e.g., READv2 outputs
        relations_df["pair"] = relations_df.apply(lambda row: tuple(sorted([row["id1"], row["id2"]])), axis=1)
        relations_df = (
            relations_df.sort_values(by=["pair", "degree"])
            .drop_duplicates(subset="pair", keep="first")
            .drop(columns=["pair"])
        )

        degree_classification_probs, relation_classification_probs = self._scale_error_rates(
            scale=self._error_rate_scale
        )

        def corrupt_relation(row: pd.Series) -> pd.Series:
            node1, node2, degree, constraints = row
            new_degree_probs = degree_classification_probs[degree]
            new_degree = random.choices(population=["1", "2", "3", "Unrelated"], weights=new_degree_probs, k=1)[0]

            new_constraints = ""
            if constraints:
                assert degree == "1"
                if new_degree == "1":
                    new_constraints_probs = relation_classification_probs[constraints]
                    new_constraints = random.choices(
                        population=["parent-child;child-parent", "siblings"], weights=new_constraints_probs, k=1
                    )[0]
                else:
                    new_constraints = ""
            return pd.Series({"id1": node1, "id2": node2, "degree": new_degree, "constraints": new_constraints})

        relations_df = relations_df.apply(corrupt_relation, axis=1)
        # Remove unrelated entries
        relations_df = relations_df[relations_df["degree"] != "Unrelated"]

        empty_iters = 0
        # Ensure at least 2 nodes and 1 relation in input data
        while len(nodes_df) < 2 or len(relations_df) < 1:
            relations_df = self._ground_truth_relations_df.sample(n=1, axis=0, random_state=empty_iters)
            relations_df = relations_df.apply(corrupt_relation, axis=1)
            relations_df = relations_df[relations_df["degree"] != "Unrelated"]
            nodes = relations_df["id1"].tolist() + relations_df["id2"].tolist()
            nodes_df = self._ground_truth_nodes_df.loc[self._ground_truth_nodes_df["id"].isin(nodes)]
            empty_iters += 1

        self._final_nodes_df = nodes_df.copy()
        self._final_relations_df = relations_df.copy()

    def run_algorithm(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            nodes_path = os.path.join(temp_dir, "nodes.csv")
            relations_path = os.path.join(temp_dir, "relations.csv")
            self._final_nodes_df.to_csv(nodes_path, index=False)
            self._final_relations_df.to_csv(relations_path, index=False)

            outputs_dir = os.path.join(temp_dir, "outputs")
            os.makedirs(outputs_dir, exist_ok=True)
            if self._epsilon is not None:
                pedigree_reconstructor = PedigreeReconstructor(
                    relations_path,
                    nodes_path,
                    outputs_dir,
                    max_candidate_pedigrees=self._max_candidate_pedigrees,
                    epsilon=self._epsilon,
                    random_seed=self._random_seed,
                    plot=False,
                )
            else:
                pedigree_reconstructor = PedigreeReconstructor(
                    relations_path,
                    nodes_path,
                    outputs_dir,
                    max_candidate_pedigrees=self._max_candidate_pedigrees,
                    random_seed=self._random_seed,
                    plot=False,
                )

            try:
                self._algorithm_pedigree = pedigree_reconstructor.find_best_pedigree()
                self._algorithm_found_pedigree = True
            except RuntimeError:
                self._algorithm_found_pedigree = False

    def get_pedigree_statistics(self) -> dict[str, int | float]:
        statistics: dict[str, int | float] = dict()
        statistics["Total Node Count"] = len(self._ground_truth_pedigree.node_to_data)
        statistics["Sampled Node Count"] = len(self._final_nodes_df)
        statistics["Proportion of Inbred Nodes"] = self._calculate_inbred_proportion()
        statistics["Proportion of Non-Final-Generation Nodes with Children"] = self._calculate_has_children_proportion()
        statistics["Mean Children Count per Parent"] = self._calculate_mean_children_per_node()
        return statistics

    def _calculate_inbred_proportion(self) -> float:
        num_nodes_with_parents = 0
        num_nodes_with_related_parents = 0

        related_pairs = set(self._ground_truth_relations_df[["id1", "id2"]].itertuples(index=False, name=None))
        for node in self._ground_truth_pedigree.node_to_data:
            father = self._ground_truth_pedigree.get_father(node)
            mother = self._ground_truth_pedigree.get_mother(node)
            if father and mother:
                num_nodes_with_parents += 1
                if (father, mother) in related_pairs or (mother, father) in related_pairs:
                    num_nodes_with_related_parents += 1
        return num_nodes_with_related_parents / num_nodes_with_parents if num_nodes_with_parents > 0 else 0

    def _calculate_has_children_proportion(self) -> float:
        num_nonleaf_parents = 0
        num_nonleaf_nodes = 0

        for node in self._ground_truth_pedigree.node_to_data:
            if self._node_to_generation[node] != max(self._node_to_generation.values()):
                num_nonleaf_nodes += 1
                if self._ground_truth_pedigree.get_children(node):
                    num_nonleaf_parents += 1
        return num_nonleaf_parents / num_nonleaf_nodes if num_nonleaf_nodes > 0 else 0

    def _calculate_mean_children_per_node(self) -> float:
        num_children = 0
        num_parents = 0

        for node in self._ground_truth_pedigree.node_to_data:
            if self._ground_truth_pedigree.get_children(node):
                num_children += len(self._ground_truth_pedigree.get_children(node))
                num_parents += 1
        return num_children / num_parents if num_parents > 0 else 0

    def get_metrics(self) -> dict[str, float]:
        metrics: dict[str, float] = dict()
        if self._algorithm_found_pedigree:
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
        else:
            metrics["Pairwise Relation Accuracy"] = 0
            metrics["Relation Precision"] = 0
            metrics["Relation Recall"] = 0
            metrics["Relation F1"] = 0
            metrics["Pairwise Degree Accuracy"] = 0
            metrics["Degree Precision"] = 0
            metrics["Degree Recall"] = 0
            metrics["Degree F1"] = 0
            metrics["Connectivity R-squared"] = 0
        return metrics

    @staticmethod
    def _calculate_tp_fp_fn(
        ground_truth_counts: defaultdict[str, int], algorithm_counts: defaultdict[str, int]
    ) -> tuple[int, int, int]:
        tp = 0  # True positives
        fp = 0  # False positives
        fn = 0  # False negatives
        relations = ground_truth_counts.keys() | algorithm_counts.keys()
        for relation in relations:
            true_count = ground_truth_counts[relation]
            algorithm_count = algorithm_counts[relation]

            if true_count == algorithm_count:
                tp += true_count
            elif true_count > algorithm_count:
                tp += algorithm_count
                fn += true_count - algorithm_count
            else:
                tp += true_count
                fp += algorithm_count - true_count
        return tp, fp, fn

    def _calculate_relation_metrics(self) -> tuple[float, float, float, float]:
        nodes: list[str] = self._final_nodes_df["id"].tolist()  # Use unmasked nodes
        correct_node_pairs: int = 0
        total_node_pairs: int = 0
        relation_tp: int = 0
        relation_fp: int = 0
        relation_fn: int = 0

        for node1, node2 in combinations(nodes, 2):
            ground_truth_relations: defaultdict[str, int] = self._ground_truth_pedigree.get_relations_between_nodes(
                node1, node2, include_maternal_paternal=True
            )
            algorithm_relations: defaultdict[str, int] = self._algorithm_pedigree.get_relations_between_nodes(
                node1, node2, include_maternal_paternal=True
            )

            if ground_truth_relations == algorithm_relations:
                correct_node_pairs += 1
            total_node_pairs += 1

            tp, fp, fn = self._calculate_tp_fp_fn(ground_truth_relations, algorithm_relations)
            relation_tp += tp
            relation_fp += fp
            relation_fn += fn

        pairwise_relation_accuracy = correct_node_pairs / total_node_pairs
        relation_precision = relation_tp / (relation_tp + relation_fp)
        relation_recall = relation_tp / (relation_tp + relation_fn)
        relation_f1 = (
            (2 * relation_precision * relation_recall) / (relation_precision + relation_recall)
            if relation_precision + relation_recall > 0
            else 0
        )
        return pairwise_relation_accuracy, relation_precision, relation_recall, relation_f1

    def _calculate_degree_metrics(self) -> tuple[float, float, float, float]:
        nodes: list[str] = self._final_nodes_df["id"].tolist()  # Use unmasked nodes
        correct_node_pairs: int = 0
        total_node_pairs: int = 0
        degree_tp: int = 0
        degree_fp: int = 0
        degree_fn: int = 0

        for node1, node2 in combinations(nodes, 2):
            ground_truth_relations: defaultdict[str, int] = self._ground_truth_pedigree.get_relations_between_nodes(
                node1, node2, include_maternal_paternal=True
            )
            algorithm_relations: defaultdict[str, int] = self._algorithm_pedigree.get_relations_between_nodes(
                node1, node2, include_maternal_paternal=True
            )

            ground_truth_degrees = defaultdict(int)
            algorithm_degrees = defaultdict(int)
            for relation in ["parent-child", "child-parent", "siblings"]:
                ground_truth_degrees["1"] += ground_truth_relations[relation]
                algorithm_degrees["1"] += algorithm_relations[relation]

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
                ground_truth_degrees["2"] += ground_truth_relations[relation]
                algorithm_degrees["2"] += algorithm_relations[relation]

            if ground_truth_degrees == algorithm_degrees:
                correct_node_pairs += 1
            total_node_pairs += 1

            tp, fp, fn = self._calculate_tp_fp_fn(ground_truth_degrees, algorithm_degrees)
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
        ground_truth_relation_counter: defaultdict[str, int] = defaultdict(int)
        algorithm_relation_counter: defaultdict[str, int] = defaultdict(int)

        nodes: list[str] = self._final_nodes_df["id"].tolist()  # Use unmasked nodes
        for node1, node2 in combinations(nodes, 2):
            ground_truth_relations: defaultdict[str, int] = self._ground_truth_pedigree.get_relations_between_nodes(
                node1, node2, include_maternal_paternal=True
            )
            algorithm_relations: defaultdict[str, int] = self._algorithm_pedigree.get_relations_between_nodes(
                node1, node2, include_maternal_paternal=True
            )

            for _relation, count in ground_truth_relations.items():
                ground_truth_relation_counter[node1] += count
                ground_truth_relation_counter[node2] += count

            for _relation, count in algorithm_relations.items():
                algorithm_relation_counter[node1] += count
                algorithm_relation_counter[node2] += count

        true_connectivities: list[int] = []
        algorithm_connectivities: list[int] = []
        for node in nodes:
            true_connectivities.append(ground_truth_relation_counter[node])
            algorithm_connectivities.append(algorithm_relation_counter[node])
        return r2_score(true_connectivities, algorithm_connectivities)
