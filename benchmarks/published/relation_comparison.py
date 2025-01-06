import tempfile
import pandas as pd
from collections import defaultdict, namedtuple
from itertools import combinations
from pedigree_package import Pedigree, PedigreeEnsemble

class RelationComparison:
    Relation = namedtuple("Relation", ["id1", "id2", "relation"])

    def __init__(self, published_relations_path: str, algorithm_nodes_path: str, algorithm_relations_path: str) -> None:
        self.published_relation_counts = self._load_published_relations(published_relations_path)
        self.algorithm_relation_counts = self._load_algorithm_relations(algorithm_nodes_path, algorithm_relations_path)

    def _load_published_relations(self, path: str) -> defaultdict[Relation, int]:
        published_relations_df = pd.read_csv(path, comment="#")
        relation_counts: defaultdict[Relation, int] = defaultdict(int)
        for id1, id2, relation in published_relations_df.itertuples(index=False):
            relation = self._sort_relation(RelationComparison.Relation(id1, id2, relation))
            relation_counts[relation] += 1
        return relation_counts

    def _load_algorithm_relations(self, nodes_path: str, relations_path: str) -> defaultdict[Relation, int]:
        self._algorithm_pedigree: Pedigree = self._run_algorithm(nodes_path, relations_path)
        algorithm_relations: defaultdict[Relation, int] = defaultdict(int)

        for id1, id2 in combinations(self._algorithm_pedigree.node_to_data, 2):
            if not id1.isnumeric() and not id2.isnumeric():  # Skip placeholder nodes
                shared_relations = self._algorithm_pedigree.get_relations_between_nodes(id1, id2, include_maternal_paternal=True)
                for shared_relation, count in shared_relations.items():
                    relation = self._sort_relation(RelationComparison.Relation(id1, id2, shared_relation))
                    algorithm_relations[relation] += count
        return algorithm_relations

    @staticmethod
    def _sort_relation(relation: Relation) -> Relation:
        FLIPPED_RELATIONS = {
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
            "2": "2"  # Symmetric
        }
        if relation.id2 < relation.id1:
            return RelationComparison.Relation(relation.id2, relation.id1, FLIPPED_RELATIONS[relation.relation])
        else:
            return relation

    @staticmethod
    def _run_algorithm(nodes_path: str, relations_path: str) -> Pedigree:
        with tempfile.TemporaryDirectory() as temp_dir:
            pedigree_ensemble = PedigreeEnsemble(relations_path,
                                                 nodes_path,
                                                 outputs_dir=temp_dir,
                                                 sample_count=100)
            return pedigree_ensemble.find_best_pedigree()

    @staticmethod
    def _fill_uncertain_relations(relation_counts: defaultdict[Relation, int]) -> defaultdict[Relation, int]:
        pass

    def get_metrics(self) -> dict[str, float]:
        metrics: dict[str, float] = dict()
        pairwise_relation_accuracy, relation_precision, relation_recall, relation_f1 = self._calculate_relation_metrics()
        # pairwise_degree_accuracy, degree_precision, degree_recall, degree_f1 = self._calculate_degree_metrics()

        metrics["Pairwise Relation Accuracy"] = pairwise_relation_accuracy
        metrics["Relation Precision"] = relation_precision
        metrics["Relation Recall"] = relation_recall
        metrics["Relation F1"] = relation_f1
        # metrics["Pairwise Degree Accuracy"] = pairwise_degree_accuracy
        # metrics["Degree Precision"] = degree_precision
        # metrics["Degree Recall"] = degree_recall
        # metrics["Degree F1"] = degree_f1
        # metrics["Connectivity R-squared"] = self._calculate_connectivity_r_squared()
        return metrics

    @staticmethod
    def _calculate_tp_fp_fn(ground_truth_counts: defaultdict(int), algorithm_counts: defaultdict(int)) -> tuple[int, int, int]:
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
        # Map node pairs to their shared relation counts to calculate pairwise relation accuracy
        published_node_pair_to_relation_counts: defaultdict[tuple[str, str], defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
        algorithm_node_pair_to_relation_counts: defaultdict[tuple[str, str], defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
        for (id1, id2, shared_relation), count in self.published_relation_counts.items():
            published_node_pair_to_relation_counts[(id1, id2)][shared_relation] += count
        for (id1, id2, shared_relation), count in self.algorithm_relation_counts.items():
            algorithm_node_pair_to_relation_counts[(id1, id2)][shared_relation] += count
        
        correct_node_pairs = 0
        total_node_pairs = 0
        for node1, node2 in combinations(self._algorithm_pedigree.node_to_data, 2):
            if not node1.isnumeric() and not node2.isnumeric():  # Skip placeholder nodes
                node1, node2 = sorted([node1, node2])
                if published_node_pair_to_relation_counts[(node1, node2)] == algorithm_node_pair_to_relation_counts[(node1, node2)]:
                    correct_node_pairs += 1
                total_node_pairs += 1
        pairwise_relation_accuracy = correct_node_pairs / total_node_pairs

        tp, fp, fn = self._calculate_tp_fp_fn(self.published_relation_counts, self.algorithm_relation_counts)
        relation_precision = tp / (tp + fp)
        relation_recall = tp / (tp + fn)
        relation_f1 = 2 * (relation_precision * relation_recall) / (relation_precision + relation_recall)
        return pairwise_relation_accuracy, relation_precision, relation_recall, relation_f1
    