import os
import copy
import random
import time
import logging
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from collections import defaultdict
from tqdm import tqdm
from typing import Any
from .pedigree import Pedigree


class PedigreeEnsemble:
    """
    Manages and builds up a collection of potential Pedigrees.
    """
    def __init__(self, relations_path: str, nodes_path: str, outputs_dir: str, sample_count: int = 1000, epsilon: float = 0.2, random_seed: Any = 42) -> None:
        self._start_time = time.time()
        self._process_node_data(nodes_path)
        self._process_relations_data(relations_path)
        self._outputs_dir = outputs_dir
        self._sample_count = sample_count  # Number of pedigrees to sample at each step of algorithm
        self._epsilon = epsilon  # Parameter for epsilon-greedy sampling when pruning pedigrees
        random.seed(random_seed)

        self._pedigrees: list[Pedigree] = [self._get_initial_pedigree()]
        self._pair_to_constraints: defaultdict[tuple[str, str], list[tuple[str, ...]]] = self._get_pair_to_constraints()

    def _process_node_data(self, nodes_path: str) -> None:
        """
        Read in node data and set types.
        """
        self._node_data = pd.read_csv(nodes_path, dtype=str, comment="#", keep_default_na=False)
        for column_name in ["id", "sex", "y_haplogroup", "mt_haplogroup", "can_have_children"]:
            if column_name not in self._node_data.columns:
                raise ValueError(f"Column \"{column_name}\" not found in input node data.")

        # Numeric IDs are used for placeholder nodes
        if self._node_data["id"].str.isnumeric().any():
            raise ValueError("Sample IDs cannot be completely numeric.")
        if not self._node_data["sex"].isin(["M", "F"]).all():
            raise ValueError("Node sex must be \"M\" or \"F\".")
        if not self._node_data["can_have_children"].isin(["True", "False", ""]).all():
            raise ValueError("can_have_children value must be \"True\", \"False\", or empty.")
        # Convert "can_have_children" column to booleans
        self._node_data["can_have_children"] = self._node_data["can_have_children"].map({"False": False, "": True})

    def _process_relations_data(self, relations_path: str) -> None:
        """
        Read in and separate relations data into first- and second-degree relations.
        """
        relations_data = pd.read_csv(relations_path, dtype=str, comment="#", keep_default_na=False)
        for column_name in ["id1", "id2", "degree", "constraints"]:
            if column_name not in relations_data.columns:
                raise ValueError(f"Column \"{column_name}\" not found in input relations data.")

        if not relations_data["id1"].isin(self._node_data["id"]).all() or not relations_data["id2"].isin(self._node_data["id"]).all():
            raise ValueError("All node IDs in relations data must be present in node data.")
        if not relations_data["degree"].isin(["1", "2", "3"]).all():
            raise ValueError("Degree must be 1, 2, or 3.")
            
        relations_data["pair_degree"] = relations_data.apply(lambda row: tuple(sorted([row["id1"], row["id2"], row["degree"]])), axis=1)
        grouped_relations = relations_data.groupby("pair_degree")
        # Check for groups with multiple non-empty constraints
        invalid_groups = grouped_relations.filter(
            lambda group: (group["constraints"] != "").sum() > 1
        )
        if not invalid_groups.empty:
            raise ValueError("Node pairs cannot have multiple non-empty constraints of the same degree.")
        relations_data.drop("pair_degree", axis=1, inplace=True)

        def split_and_validate_constraints(constraints: str) -> None:
            allowed_constraints = {
                "parent-child", "child-parent", "siblings", 
                "aunt/uncle-nephew/niece", "nephew/niece-aunt/uncle", 
                "grandparent-grandchild", "grandchild-grandparent", 
                "half-siblings"
            }
            if constraints:
                constraints_list = [c for c in constraints.split(";")]
                if any(c not in allowed_constraints for c in constraints_list):
                    raise ValueError(f"Invalid constraints found: {invalid_constraints}")
        relations_data["constraints"].apply(split_and_validate_constraints)

        def sort_nodes(row: pd.Series) -> pd.Series:
            """
            Ensure id1 and id2 are in a fixed (sorted) order and flip constraints as needed.
            """
            flipped_pair_to_constraints = {  # Mapping for flipping constraints
                "parent-child": "child-parent",
                "child-parent": "parent-child",
                "aunt/uncle-nephew/niece": "nephew/niece-aunt/uncle",
                "nephew/niece-aunt/uncle": "aunt/uncle-nephew/niece",
                "grandparent-grandchild": "grandchild-grandparent",
                "grandchild-grandparent": "grandparent-grandchild",
                "siblings": "siblings",  # Symmetric
                "half-siblings": "half-siblings"  # Symmetric
            }
            if row["id2"] < row["id1"]:
                constraints = row["constraints"]
                if constraints:  # Split constraints and map each to its flipped value
                    constraints_list = [c.strip() for c in constraints.split(";")]
                    flipped = [flipped_pair_to_constraints[c] for c in constraints_list]
                    flipped_constraints = ";".join(flipped)
                else:
                    flipped_constraints = ""
                # Swap id1 and id2, and flip constraints
                return pd.Series({"id1": row["id2"], "id2": row["id1"], "degree": row["degree"], "constraints": flipped_constraints})
            else:
                return row
        relations_data = relations_data.apply(sort_nodes, axis=1)

        self._DEFAULT_CONSTRAINTS = {
            "1": "parent-child;child-parent;siblings",
            "2": "aunt/uncle-nephew/niece;nephew/niece-aunt/uncle;grandparent-grandchild;grandchild-grandparent;half-siblings",
            "3": "half aunt/uncle-half nephew/niece;half nephew/niece-half aunt/uncle;greatgrandparent-greatgrandchild;greatgrandchild-greatgrandparent;grandaunt/granduncle-grandnephew/grandniece;grandnephew/grandniece-grandaunt/granduncle;first cousins"
        }

        def fill_constraints(row: pd.Series) -> pd.Series:
            if not row["constraints"]:
                constraints = self._DEFAULT_CONSTRAINTS[row["degree"]]
                return pd.Series({"id1": row["id1"], "id2": row["id2"], "degree": row["degree"], "constraints": constraints})
            return row
        relations_data = relations_data.apply(fill_constraints, axis=1)
        
        self._first_degree_relations = relations_data[relations_data["degree"] == "1"].reset_index(drop=True)
        self._second_degree_relations = relations_data[relations_data["degree"] == "2"].reset_index(drop=True)
        self._third_degree_relations = relations_data[relations_data["degree"] == "3"].reset_index(drop=True)
        self._first_and_second_degree_relations = pd.concat([self._first_degree_relations, self._second_degree_relations]).reset_index(drop=True)
        self._all_relations = pd.concat([self._first_degree_relations, self._second_degree_relations, self._third_degree_relations]).reset_index(drop=True)

    def _get_initial_pedigree(self):
        """
        Create the initial pedigree and add all nodes.
        """
        initial_pedigree = Pedigree()
        for node_id, sex, y_haplogroup, mt_haplogroup, can_have_children in self._node_data.itertuples(index=False):
            initial_pedigree.add_node(node_id, sex, y_haplogroup, mt_haplogroup, can_have_children)
        return initial_pedigree

    def find_best_pedigree(self) -> Pedigree:
        """
        Finds the configuration of relations that yields the "best" pedigree.
        Writes to output_dir the set of relations with the least changes from the original input data.
        """
        progress_bar = tqdm(
            self._first_and_second_degree_relations.iterrows(),
            total=self._first_and_second_degree_relations.shape[0],
            smoothing=0.5,
            bar_format="{l_bar}{bar} | {n_fmt}/{total_fmt} [{elapsed}<{remaining}]"
        )
        for idx, row in progress_bar:
            node1, node2, degree, _ = row
            logging.info(f"Current relation: {node1}, {node2}, {degree}")
            progress_bar.set_description(f"Processing relation {{{node1}, {node2}, {degree}}}")
            self._add_relation(node1, node2, degree=degree)
            self._clean_relation_dicts()

            relations_so_far = self._all_relations.iloc[:idx+1]
            if degree == "1" and len(relations_so_far) < len(self._first_and_second_degree_relations):
                self._prune_pedigrees(relations_so_far, check_half_siblings=False)
            else:
                self._prune_pedigrees(relations_so_far, check_half_siblings=True)
            logging.info(f"Remaining pedigrees after pruning: {len(self._pedigrees)}\t\tElapsed: {round(time.time() - self._start_time, 1)} s\n")

        self._final_pedigree.clean_up_relations()
        self._write_corrected_relations()
        return self._final_pedigree

    def _add_relation(self, node1: str, node2: str, degree: str) -> None:
        """
        Connects two nodes in every pedigree. 
        Does not care about input constraints; those will be handled by Pedigree.count_inconsistencies().
        """
        assert degree in ["1", "2"]

        new_pedigrees: list[Pedigree] = []
        for pedigree in self._pedigrees:
            if degree == "1":
                new_pedigrees.extend(PedigreeEnsemble._connect_first_degree_relation(pedigree, node1, node2, constraints=self._DEFAULT_CONSTRAINTS["1"]))
                new_pedigrees.extend(PedigreeEnsemble._connect_second_degree_relation(pedigree, node1, node2, constraints=self._DEFAULT_CONSTRAINTS["2"]))
            elif degree == "2":
                new_pedigrees.append(pedigree)  # No relation (i.e. false positive)
                new_pedigrees.extend(PedigreeEnsemble._connect_first_degree_relation(pedigree, node1, node2, constraints=self._DEFAULT_CONSTRAINTS["1"]))
                new_pedigrees.extend(PedigreeEnsemble._connect_second_degree_relation(pedigree, node1, node2, constraints=self._DEFAULT_CONSTRAINTS["2"]))
        self._pedigrees = new_pedigrees

    @staticmethod
    def _connect_first_degree_relation(pedigree: Pedigree, node1: str, node2: str, constraints: str) -> None:
        """
        Update pedigree with a first-degree relation.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data

        new_pedigrees: list[Pedigree] = []
        possible_relations: list[str] = constraints.split(";")

        for relation in possible_relations:
            if relation == "parent-child":
                new_pedigrees.extend(PedigreeEnsemble._connect_parent_relation(pedigree, node1, node2))
            if relation == "child-parent":
                new_pedigrees.extend(PedigreeEnsemble._connect_parent_relation(pedigree, node2, node1))
            if relation == "siblings":
                new_pedigrees.extend(PedigreeEnsemble._connect_sibling_relation(pedigree, node1, node2))
        return new_pedigrees

    @staticmethod
    def _connect_second_degree_relation(pedigree: Pedigree, node1: str, node2: str, constraints: str) -> None:
        """
        Update pedigree with a second-degree relation.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data

        new_pedigrees: list[Pedigree] = []
        possible_relations: list[str] = constraints.split(";")

        for relation in possible_relations:
            if relation == "aunt/uncle-nephew/niece":
                new_pedigrees.extend(PedigreeEnsemble._connect_aunt_uncle_relation(pedigree, node1, node2))
            if relation == "nephew/niece-aunt/uncle":
                new_pedigrees.extend(PedigreeEnsemble._connect_aunt_uncle_relation(pedigree, node2, node1))
            if relation == "grandparent-grandchild":
                new_pedigrees.extend(PedigreeEnsemble._connect_grandparent_relation(pedigree, node1, node2))
            if relation == "grandchild-grandparent":
                new_pedigrees.extend(PedigreeEnsemble._connect_grandparent_relation(pedigree, node2, node1))
            if relation == "half-siblings":
                new_pedigrees.extend(PedigreeEnsemble._connect_half_sibling_relation(pedigree, node1, node2))
        return new_pedigrees

    @staticmethod
    def _connect_parent_relation(pedigree: Pedigree, node1: str, node2: str) -> list[Pedigree]:
        """
        Adds a parent-child relation and merges nodes appropriately.
        Returns a list containing the resulting Pedigree, if successful.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node2)
        original_parent: str
        if new_pedigree.node_to_data[node1]["sex"] == "M":
            original_parent = new_pedigree.node_to_father[node2]
        else:
            original_parent = new_pedigree.node_to_mother[node2]

        if not new_pedigree.check_cycles_if_merged(node1, original_parent):
            new_pedigree.merge_nodes(node1, original_parent)
            ret.append(new_pedigree)
        return ret

    @staticmethod
    def _connect_sibling_relation(pedigree: Pedigree, node1: str, node2: str) -> list[Pedigree]:
        """
        Adds a sibling relation and merges nodes appropriately.
        Returns a list containing the resulting Pedigree, if successful.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node1)
        new_pedigree.fill_node_parents(node2)

        father1 = new_pedigree.node_to_father[node1]
        father2 = new_pedigree.node_to_father[node2]
        if not new_pedigree.check_cycles_if_merged(father1, father2):
            new_pedigree.merge_nodes(father1, father2)
            mother1 = new_pedigree.node_to_mother[node1]
            mother2 = new_pedigree.node_to_mother[node2]
            if not new_pedigree.check_cycles_if_merged(mother1, mother2):
                new_pedigree.merge_nodes(mother1, mother2)
                new_pedigree.add_sibling_relation(node1, node2)
                ret.append(new_pedigree)
        return ret

    @staticmethod
    def _connect_aunt_uncle_relation(pedigree: Pedigree, node1: str, node2: str, shared_relative_sex: str | None = None) -> list[Pedigree]:
        """
        Adds an aunt/uncle-nephew/niece relation and merges nodes appropriately.
        Returns a list containing the resulting Pedigree(s), if successful.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data
        assert shared_relative_sex in ["M", "F", None]

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node2)

        node2_parents: list[str]
        if shared_relative_sex == "M":
            node2_parents = [new_pedigree.node_to_father[node2]]
        elif shared_relative_sex == "F":
            node2_parents = [new_pedigree.node_to_mother[node2]]
        else:
            node2_parents = [new_pedigree.node_to_father[node2], new_pedigree.node_to_mother[node2]]

        for node2_parent in node2_parents:
            if node1 != node2_parent:
                ret.extend(PedigreeEnsemble._connect_sibling_relation(new_pedigree, node1, node2_parent))
        return ret

    @staticmethod
    def _connect_grandparent_relation(pedigree: Pedigree, node1: str, node2: str, shared_relative_sex: str | None = None) -> list[Pedigree]:
        """
        Adds a grandparent-grandchild relation and merges nodes appropriately.
        Returns a list containing the resulting Pedigree(s), if successful.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data
        assert shared_relative_sex in ["M", "F", None]

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node2)

        node2_parents: list[str]
        if shared_relative_sex == "M":
            node2_parents = [new_pedigree.node_to_father[node2]]
        elif shared_relative_sex == "F":
            node2_parents = [new_pedigree.node_to_mother[node2]]
        else:
            node2_parents = [new_pedigree.node_to_father[node2], new_pedigree.node_to_mother[node2]]

        for node2_parent in node2_parents:
            if node1 != node2_parent:
                ret.extend(PedigreeEnsemble._connect_parent_relation(new_pedigree, node1, node2_parent))
        return ret

    @staticmethod
    def _connect_half_sibling_relation(pedigree: Pedigree, node1: str, node2: str, shared_relative_sex: str | None = None) -> list[Pedigree]:
        """
        Adds a half-sibling relation and merges nodes appropriately.
        Returns a list containing the resulting Pedigree(s), if successful.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node1)
        new_pedigree.fill_node_parents(node2)

        node1_parents: list[str]
        node2_parents: list[str]
        if shared_relative_sex == "M":
            node1_parents = [new_pedigree.node_to_father[node1]]
            node2_parents = [new_pedigree.node_to_father[node2]]
        elif shared_relative_sex == "F":
            node1_parents = [new_pedigree.node_to_mother[node1]]
            node2_parents = [new_pedigree.node_to_mother[node2]]
        else:
            node1_parents = [new_pedigree.node_to_father[node1], new_pedigree.node_to_mother[node1]]
            node2_parents = [new_pedigree.node_to_father[node2], new_pedigree.node_to_mother[node2]]

        # Node 1 and Node 2 are half-siblings via one of Node 1's parents
        for node1_parent in node1_parents:
            if node1_parent != node2:
                ret.extend(PedigreeEnsemble._connect_parent_relation(new_pedigree, node1_parent, node2))
        # Node 1 and Node 2 are half-siblings via one of Node 2's parents
        for node2_parent in node2_parents:
            if node2_parent != node1:
                ret.extend(PedigreeEnsemble._connect_parent_relation(new_pedigree, node2_parent, node1))
        return ret

    def _clean_relation_dicts(self) -> None:
        """
        Remove unnecessary entries in Pedigree dicts.
        """
        for pedigree in self._pedigrees:
            pedigree.clean_up_relations()

    def _get_pair_to_constraints(self) -> defaultdict[tuple[str, str], list[tuple[str, ...]]]:
        """
        Turn DataFrame of relations/constraints into dict(s) of {node pairs: list of possible relations}.
        Dict values are lists of tuples (as opposed to just tuples) because a pair of nodes can share more than 1 relation.
        """
        pair_to_constraints: defaultdict[tuple[str, str], list[tuple[str, ...]]] = defaultdict(list)
        for node1, node2, _, constraints in self._all_relations.itertuples(index=False):
            pair_to_constraints[(node1, node2)].append(tuple(constraints.split(";")))
        for node_pair in pair_to_constraints:
            pair_to_constraints[node_pair].sort(key=lambda x: len(x))  # Sort by number of constraints so specific constraints are checked first when pruning
        return pair_to_constraints

    def _prune_pedigrees(self, relations_so_far: pd.DataFrame, check_half_siblings: bool) -> None:
        """
        Remove pedigrees with inconsistencies.
        """
        # Group all input relations (so far) by node pair
        relations_so_far_dict: defaultdict[tuple[str, str], tuple[str, str]] = defaultdict(list)
        for node1, node2, degree, constraints in relations_so_far.itertuples(index=False):
            relations_so_far_dict[(node1, node2)].append((degree, constraints))

        seen_topologies = set()
        new_potential_pedigrees = []
        for pedigree in self._pedigrees:
            if pedigree.validate_members(set(self._node_data["id"])) and pedigree.validate_can_have_children():
                pedigree.update_haplogroups()
                if pedigree.validate_haplogroups():
                    topology = pedigree.get_topo_sort()
                    if topology not in seen_topologies:
                        new_potential_pedigrees.append(pedigree)
                        seen_topologies.add(topology)
        random.shuffle(new_potential_pedigrees)  # Shuffle to avoid ordering bias in epsilon-greedy sampling

        strikes = []
        third_degree_strikes = []
        counts = defaultdict(int)
        for pedigree in new_potential_pedigrees:
            num_strikes, _ = pedigree.count_inconsistencies(self._pair_to_constraints, relations_so_far_dict, check_half_siblings)
            num_third_degree_strikes = pedigree.count_third_degree_inconcistencies(self._pair_to_constraints)
            strikes.append(num_strikes)
            third_degree_strikes.append(num_third_degree_strikes)
            counts[num_strikes] += 1
        logging.info(f"Strike counts before pruning: {str(dict(sorted(counts.items())))}")

        def epsilon_greedy_sample(pedigrees: list[Pedigree], strikes: list[int], third_degree_strikes: list[int], epsilon: float, sample_count: int) -> list[Pedigree]:
            assert len(pedigrees) == len(strikes)
            if len(pedigrees) <= sample_count:
                return pedigrees
            
            sorted_pedigrees = [pedigree for pedigree, _, _ in sorted(zip(pedigrees, strikes, third_degree_strikes), key=lambda x: (x[1], x[2]))]
            exploration_sample_count = int(epsilon * sample_count)
            exploitation_sample_count = sample_count - exploration_sample_count

            exploitation_pedigrees = sorted_pedigrees[:exploitation_sample_count]
            exploration_pedigrees = random.sample(sorted_pedigrees[exploitation_sample_count:], exploration_sample_count)
            return exploitation_pedigrees + exploration_pedigrees

        if len(relations_so_far) < len(self._first_and_second_degree_relations):
            self._pedigrees = epsilon_greedy_sample(new_potential_pedigrees, strikes, third_degree_strikes, epsilon=self._epsilon, sample_count=self._sample_count)
        else:  # Final iteration
            best_pedigrees = [pedigree for pedigree, num_strikes in zip(new_potential_pedigrees, strikes) if num_strikes == min(strikes)]
            third_degree_strikes = [pedigree.count_third_degree_inconcistencies(self._pair_to_constraints) for pedigree in best_pedigrees]  # Use 3rd-degree strikes as tiebreaker
            self._final_pedigree = best_pedigrees[third_degree_strikes.index(min(third_degree_strikes))]
            self._final_strike_count, self._final_strike_log = self._final_pedigree.count_inconsistencies(self._pair_to_constraints, relations_so_far_dict, check_half_siblings=True)

    def _write_corrected_relations(self) -> None:
        """
        Write corrected relations to file.
        """
        added_relations = []
        removed_relations = []
        for node1, node2, degree, constraints in self._final_strike_log:
            if degree[0] == "+":
                added_relations.append((node1, node2, degree[1], constraints))
            else:
                removed_relations.append((node1, node2, degree[1], constraints))
        removed_relations_set = set(removed_relations)
        
        # Separate out *changed* relations (added relation + removed relation pair, e.g., 1st-degree -> 2nd-degree)
        changed_node_pairs = set()
        for (add_node1, add_node2, _, _) in added_relations:
            for (remove_node1, remove_node2, _, _) in removed_relations:
                if (add_node1 == remove_node1 and add_node2 == remove_node2) or (add_node2 == remove_node1 and add_node1 == remove_node2):
                    changed_node_pairs.add((add_node1, add_node2))
        
        path = os.path.join(self._outputs_dir, "corrected_relations.csv")
        with open(path, "w") as file:
            file.write("id1,id2,degree,constraints\n")  # Header line
            file.write(f"# Final strike count: {self._final_strike_count}\n")

            def write_relations_line(node1, node2, degree, constraints, commented=False):
                if constraints == self._DEFAULT_CONSTRAINTS[degree]:
                    constraints = ""  # Don't write default constraints to file
                if commented:
                    file.write("# ")
                file.write(f"{node1},{node2},{degree},{constraints}\n")

            file.write("# Added relations\n")
            for node1, node2, degree, constraints in sorted(added_relations):  # Sort for consistency
                if (node1, node2) not in changed_node_pairs and (node2, node1) not in changed_node_pairs:
                    write_relations_line(node1, node2, degree, constraints)

            file.write("\n# Removed relations\n")
            for node1, node2, degree, constraints in sorted(removed_relations):
                if (node1, node2) not in changed_node_pairs and (node2, node1) not in changed_node_pairs:
                    write_relations_line(node1, node2, degree, constraints, commented=True)
            
            file.write("\n# Changed relations\n")
            for node1, node2 in sorted(changed_node_pairs):  # Pair up changed relations (add + remove)
                node1_to_write = None  # So we write the two nodes in the correct (original) order
                node2_to_write = None
                for node1_remove, node2_remove, degree_remove, constraints_remove in removed_relations:
                    if (node1_remove, node2_remove) == (node1, node2) or (node2_remove, node1_remove) == (node1, node2):
                        write_relations_line(node1_remove, node2_remove, degree_remove, constraints_remove, commented=True)
                        node1_to_write = node1_remove  # Because the removed nodes follow the original input order
                        node2_to_write = node2_remove
                for node1_add, node2_add, degree_add, constraints_add in added_relations:
                    if (node1_add, node2_add) == (node1, node2) or (node2_add, node1_add) == (node1, node2):
                        assert node1_to_write and node2_to_write
                        write_relations_line(node1_to_write, node2_to_write, degree_add, constraints_add)
            
            file.write("\n# Unchanged relations\n")
            for node1, node2, degree, constraints in self._all_relations.itertuples(index=False):
                if (node1, node2, degree, constraints) not in removed_relations_set:
                    assert (node2, node1, degree, constraints) not in removed_relations_set
                    write_relations_line(node1, node2, degree, constraints)
