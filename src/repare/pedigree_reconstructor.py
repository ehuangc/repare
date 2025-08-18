import copy
import logging
import os
import random
import time
from collections import defaultdict
from typing import Any

import pandas as pd
from tqdm import tqdm

from repare.pedigree import Pedigree

logger = logging.getLogger(__name__)


class PedigreeReconstructor:
    """
    Manages and builds up a collection of potential Pedigrees.
    """

    def __init__(
        self,
        relations_path: str,
        nodes_path: str,
        outputs_dir: str,
        max_candidate_pedigrees: int = 1000,
        epsilon: float = 0.2,
        plot: bool = True,
        write_alternate_pedigrees: bool = False,
        random_seed: Any = 42,
    ) -> None:
        self._start_time = time.time()
        self._validate_node_data(nodes_path)
        self._process_node_data()
        self._validate_relation_data(relations_path)
        self._process_relation_data()

        self._outputs_dir = outputs_dir
        # Number of pedigrees to downsample to after each iteration of algorithm
        self._max_candidate_pedigrees = max_candidate_pedigrees
        # Parameter for epsilon-greedy sampling when pruning pedigrees
        self._epsilon = epsilon
        # Whether to plot the reconstructed pedigree(s)
        self._plot = plot
        # Whether to write corrected relations and plots of alternate final pedigrees
        self._write_alternate_pedigrees = write_alternate_pedigrees
        self._random_seed = random_seed
        random.seed(self._random_seed)
        self._validate_arguments()

        # Maximum number of times to run the algorithm if no valid pedigree is found
        self._MAX_RUNS = 10
        self._candidate_pedigrees: list[Pedigree] = [self._get_initial_pedigree()]
        self._pair_to_constraints: defaultdict[tuple[str, str], list[tuple[str, ...]]] = self._get_pair_to_constraints()
        self._final_pedigrees: list[Pedigree] = []

    def _validate_node_data(self, nodes_path: str) -> None:
        """
        Validate node data input.
        """
        self._node_data = pd.read_csv(nodes_path, dtype=str, comment="#", keep_default_na=False)
        for mandatory_column in ["id", "sex", "y_haplogroup", "mt_haplogroup"]:
            if mandatory_column not in self._node_data.columns:
                raise ValueError(f'Column "{mandatory_column}" not found in input node data.')

        for optional_column in ["can_have_children", "can_be_inbred", "years_before_present"]:
            if optional_column not in self._node_data.columns:
                self._node_data[optional_column] = ""

        # Numeric IDs are used for placeholder nodes
        if self._node_data["id"].str.isnumeric().any():
            raise ValueError("Sample IDs cannot be completely numeric.")

        if self._node_data["id"].duplicated().any():
            raise ValueError("Sample IDs must be unique.")

        if self._node_data["id"].str.strip().eq("").any():
            raise ValueError("Sample IDs cannot be empty.")

        if not self._node_data["sex"].isin(["M", "F"]).all():
            raise ValueError('Node sex must be "M" or "F".')

        for haplogroup_column in ["y_haplogroup", "mt_haplogroup"]:
            for haplogroup in self._node_data[haplogroup_column]:
                if "*" in haplogroup[:-1]:
                    raise ValueError(
                        "Expandable haplogroups should contain one trailing asterisk. "
                        "No other asterisks are allowed in haplogroups."
                    )

        if not self._node_data["can_have_children"].isin(["True", "False", ""]).all():
            raise ValueError('can_have_children value must be "True", "False", or empty.')
        if not self._node_data["can_be_inbred"].isin(["True", "False", ""]).all():
            raise ValueError('can_be_inbred value must be "True", "False", or empty.')
        if not self._node_data["years_before_present"].apply(lambda x: x.isnumeric() or x == "").all():
            raise ValueError("years_before_present value must be integer or empty.")

    def _process_node_data(self) -> None:
        """
        Process node data input.
        """
        # Reorder node data columns and remove unnecessary columns
        self._node_data = self._node_data[
            ["id", "sex", "y_haplogroup", "mt_haplogroup", "can_have_children", "can_be_inbred", "years_before_present"]
        ]
        # Convert "can_have_children" and "can_be_inbred" columns to booleans
        self._node_data["can_have_children"] = self._node_data["can_have_children"].map(
            {"False": False, "True": True, "": True}
        )
        self._node_data["can_be_inbred"] = self._node_data["can_be_inbred"].map(
            {"False": False, "True": True, "": True}
        )
        # Convert "years_before_present" column to floats
        self._node_data["years_before_present"] = pd.to_numeric(
            self._node_data["years_before_present"], errors="coerce"
        )

    def _validate_relation_data(self, relations_path: str) -> None:
        """
        Validate relation data input.
        """
        self._relation_data = pd.read_csv(relations_path, dtype=str, comment="#", keep_default_na=False)
        for column_name in ["id1", "id2", "degree", "constraints"]:
            if column_name not in self._relation_data.columns:
                raise ValueError(f'Column "{column_name}" not found in input relation data.')

        for optional_column in ["force_constraints"]:
            if optional_column not in self._relation_data.columns:
                self._relation_data[optional_column] = ""

        excess_relation_nodes = set(self._relation_data["id1"]).union(set(self._relation_data["id2"])) - set(
            self._node_data["id"]
        )
        if excess_relation_nodes:
            raise ValueError(f"All node IDs in relation data must be present in node data: {excess_relation_nodes}.")

        if not self._relation_data["degree"].isin(["1", "2", "3"]).all():
            raise ValueError("Degree must be 1, 2, or 3.")
        if not self._relation_data["force_constraints"].isin(["True", "False", ""]).all():
            raise ValueError('can_have_children value must be "True", "False", or empty.')

        self._relation_data["pair_degree"] = self._relation_data.apply(
            lambda row: tuple(sorted([row["id1"], row["id2"], row["degree"]])), axis=1
        )
        grouped_relations = self._relation_data.groupby("pair_degree")
        # Check for groups with multiple non-empty constraints, which can lead to issues when counting inconsistencies
        invalid_groups = grouped_relations.filter(lambda group: (group["constraints"] != "").sum() > 1)
        if not invalid_groups.empty:
            raise ValueError("Node pairs cannot have multiple non-empty constraints of the same degree.")
        self._relation_data.drop("pair_degree", axis=1, inplace=True)

        def split_and_validate_constraints(constraints: str) -> None:
            allowed_constraints = {
                "parent-child",
                "child-parent",
                "siblings",
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
            }
            if constraints:
                constraints_list = [c for c in constraints.split(";")]
                if any(c not in allowed_constraints for c in constraints_list):
                    raise ValueError(
                        f"Invalid constraints found: {[c for c in constraints_list if c not in allowed_constraints]}."
                    )

        self._relation_data["constraints"].apply(split_and_validate_constraints)

    def _process_relation_data(self) -> None:
        """
        Process relation data input.
        """
        # Reorder relation data columns and remove unnecessary columns
        self._relation_data = self._relation_data[["id1", "id2", "degree", "constraints", "force_constraints"]]
        # Convert "force_constrains" column to booleans
        self._relation_data["force_constraints"] = self._relation_data["force_constraints"].map(
            {"False": False, "True": True, "": False}
        )

        def sort_nodes(row: pd.Series) -> pd.Series:
            """
            Ensure id1 and id2 are in a fixed (sorted) order and flip constraints as needed.
            """
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
            if row["id2"] < row["id1"]:
                constraints = row["constraints"]
                # Split constraints and map each to its flipped value
                if constraints:
                    constraints_list = [c.strip() for c in constraints.split(";")]
                    flipped = [flipped_constraints[c] for c in constraints_list]
                    relation_flipped_constraints = ";".join(flipped)
                else:
                    relation_flipped_constraints = ""
                # Swap id1 and id2, and flip constraints
                return pd.Series(
                    {
                        "id1": row["id2"],
                        "id2": row["id1"],
                        "degree": row["degree"],
                        "constraints": relation_flipped_constraints,
                        "force_constraints": row["force_constraints"],
                    }
                )
            else:
                return row

        self._relation_data = self._relation_data.apply(sort_nodes, axis=1)

        # Note: We don't use maternal/paternal 3rd-degree relations because those are not well-defined
        self._DEFAULT_CONSTRAINTS = {
            "1": ("parent-child;child-parent;siblings"),
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
                "paternal half-siblings"
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

        def fill_constraints(row: pd.Series) -> pd.Series:
            if not row["constraints"]:
                constraints = self._DEFAULT_CONSTRAINTS[row["degree"]]
                return pd.Series(
                    {
                        "id1": row["id1"],
                        "id2": row["id2"],
                        "degree": row["degree"],
                        "constraints": constraints,
                        "force_constraints": row["force_constraints"],
                    }
                )
            return row

        self._relation_data = self._relation_data.apply(fill_constraints, axis=1)

        self._first_degree_relations = self._relation_data[self._relation_data["degree"] == "1"].reset_index(drop=True)
        self._second_degree_relations = self._relation_data[self._relation_data["degree"] == "2"].reset_index(drop=True)
        self._third_degree_relations = self._relation_data[self._relation_data["degree"] == "3"].reset_index(drop=True)
        self._first_and_second_degree_relations = pd.concat(
            [self._first_degree_relations, self._second_degree_relations]
        ).reset_index(drop=True)
        self._all_relations = pd.concat(
            [self._first_degree_relations, self._second_degree_relations, self._third_degree_relations]
        ).reset_index(drop=True)

    def _validate_arguments(self) -> None:
        """
        Validate constructor arguments.
        """
        if not isinstance(self._max_candidate_pedigrees, int) or self._max_candidate_pedigrees <= 0:
            raise ValueError("max_candidate_pedigrees must be a positive integer.")
        if not (0 <= self._epsilon <= 1):
            raise ValueError("epsilon must be between 0 and 1.")

    def _shuffle_relations(self) -> None:
        """
        Shuffle relation DataFrames (when we want to restart the algorithm).
        """
        self._first_degree_relations = self._first_degree_relations.sample(
            frac=1, random_state=self._random_seed
        ).reset_index(drop=True)
        self._second_degree_relations = self._second_degree_relations.sample(
            frac=1, random_state=self._random_seed
        ).reset_index(drop=True)
        self._third_degree_relations = self._third_degree_relations.sample(
            frac=1, random_state=self._random_seed
        ).reset_index(drop=True)
        self._first_and_second_degree_relations = pd.concat(
            [self._first_degree_relations, self._second_degree_relations]
        ).reset_index(drop=True)
        self._all_relations = pd.concat(
            [self._first_degree_relations, self._second_degree_relations, self._third_degree_relations]
        ).reset_index(drop=True)

    def _get_initial_pedigree(self):
        """
        Create the initial pedigree and add all nodes.
        """
        initial_pedigree = Pedigree()
        for (
            node_id,
            sex,
            y_haplogroup,
            mt_haplogroup,
            can_have_children,
            can_be_inbred,
            years_before_present,
        ) in self._node_data.itertuples(index=False):
            initial_pedigree.add_node(
                node_id, sex, y_haplogroup, mt_haplogroup, can_have_children, can_be_inbred, years_before_present
            )
        return initial_pedigree

    def find_best_pedigree(self) -> Pedigree:
        """
        Finds the configuration of relations that yields the "best" pedigree.
        Writes to output_dir the set of relations with the least changes from the original input data.
        """
        for _ in range(self._MAX_RUNS):
            progress_bar = tqdm(
                self._first_and_second_degree_relations.iterrows(),
                total=self._first_and_second_degree_relations.shape[0],
                smoothing=0.5,
                bar_format="{l_bar}{bar} | {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
            )
            for idx, row in progress_bar:
                node1, node2, degree, constraints, force_constraints = row
                logger.info(f"Current relation: {node1}, {node2}, {degree}")
                progress_bar.set_description(f"Processing relation {{{node1}, {node2}, {degree}}}")
                self._add_relation(
                    node1, node2, degree=degree, constraints=constraints, force_constraints=force_constraints
                )
                self._clean_pedigree_data()
                self._validate_pedigree_structures()

                processed_relations = self._all_relations.iloc[: idx + 1]
                pair_to_relations_so_far = self._get_pair_to_relations_so_far(processed_relations)
                if degree == "1" and len(processed_relations) < len(self._first_and_second_degree_relations):
                    # Don't check for extraneous half-sibling relations because
                    # the 2 non-shared parents might be "merged" later
                    self._prune_pedigrees(pair_to_relations_so_far, check_half_siblings=False)
                else:
                    self._prune_pedigrees(pair_to_relations_so_far, check_half_siblings=True)
                logger.info(
                    f"Remaining pedigrees after pruning: {len(self._candidate_pedigrees)}"
                    "\t\tElapsed: "
                    f"{round(time.time() - self._start_time, 1)} s\n"
                )

            if not self._final_pedigrees:
                logger.warning("No valid pedigree found. Shuffling relations and restarting algorithm.\n")
                self._candidate_pedigrees = [self._get_initial_pedigree()]
                self._shuffle_relations()
            else:
                break

        if not self._final_pedigrees:
            logger.error(f"No valid pedigree found after {self._MAX_RUNS} runs. Exiting.")
            raise RuntimeError(f"No valid pedigree found after {self._MAX_RUNS} runs.")

        self._clean_pedigree_data()
        # Plot and write outputs of sample pedigree
        sample_idx = random.randint(0, len(self._final_pedigrees) - 1)
        self._sample_pedigree = self._final_pedigrees[sample_idx]
        self._sample_strike_count = self._final_strike_counts[sample_idx]
        self._sample_strike_log = self._final_strike_logs[sample_idx]
        self._write_corrected_input_relations(
            self._sample_strike_count,
            self._sample_strike_log,
            os.path.join(self._outputs_dir, "corrected_input_relations.csv"),
        )
        self._sample_pedigree.write_exact_relations(
            os.path.join(self._outputs_dir, "reconstructed_exact_relations.csv")
        )
        if self._plot:
            try:
                self._sample_pedigree.plot(os.path.join(self._outputs_dir, "reconstructed_pedigree.pdf"))
                pygraphviz_found = True
            except ImportError:
                logger.warning(
                    "PyGraphviz (https://pygraphviz.github.io/) must be installed to plot pedigrees. "
                    "Skipping plotting of reconstructed pedigree(s)."
                )
                pygraphviz_found = False

        # Plot and write outputs of alternate pedigrees
        if self._write_alternate_pedigrees:
            os.makedirs(os.path.join(self._outputs_dir, "alternate_pedigrees"), exist_ok=True)
            for idx, (pedigree, strike_count, strike_log) in enumerate(
                zip(self._final_pedigrees, self._final_strike_counts, self._final_strike_logs, strict=True)
            ):
                self._write_corrected_input_relations(
                    strike_count,
                    strike_log,
                    os.path.join(
                        self._outputs_dir, "alternate_pedigrees", f"pedigree_{idx}_corrected_input_relations.csv"
                    ),
                )
                pedigree.write_exact_relations(
                    os.path.join(self._outputs_dir, "alternate_pedigrees", f"pedigree_{idx}_exact_relations.csv")
                )
                if self._plot and pygraphviz_found:
                    pedigree.plot(os.path.join(self._outputs_dir, "alternate_pedigrees", f"pedigree_{idx}.png"))
        return self._sample_pedigree

    @staticmethod
    def _check_haplogroups(haplogroup1: str, haplogroup2: str) -> bool:
        """
        Checks if two haplogroups are compatible. Same semantics as pedigree.validate_haplogroups().
        "*" is wild card character.
        """
        if not haplogroup1 or not haplogroup2:  # empty OK
            return True
        haplogroup1_stripped, haplogroup2_stripped = haplogroup1.rstrip("*"), haplogroup2.rstrip("*")
        return haplogroup1_stripped.startswith(haplogroup2_stripped) or haplogroup2_stripped.startswith(
            haplogroup1_stripped
        )

    @staticmethod
    def _check_parent_child_haplogroups(pedigree: Pedigree, parent: str, child: str) -> bool:
        """
        Checks if the haplogroups of a parent and child are compatible.
        """
        if pedigree.get_data(parent)["sex"] == "M" and pedigree.get_data(child)["sex"] == "M":
            return PedigreeReconstructor._check_haplogroups(
                pedigree.get_data(parent)["y_haplogroup"], pedigree.get_data(child)["y_haplogroup"]
            )
        if pedigree.get_data(parent)["sex"] == "F":
            return PedigreeReconstructor._check_haplogroups(
                pedigree.get_data(parent)["mt_haplogroup"], pedigree.get_data(child)["mt_haplogroup"]
            )
        return True

    @staticmethod
    def _check_sibling_haplogroups(pedigree: Pedigree, sibling1: str, sibling2: str) -> bool:
        """
        Checks if the haplogroups of two full siblings are compatible.
        """
        if pedigree.get_data(sibling1)["sex"] == "M" and pedigree.get_data(sibling2)["sex"] == "M":
            # MT haplogroups still need to agree as well
            if not PedigreeReconstructor._check_haplogroups(
                pedigree.get_data(sibling1)["y_haplogroup"], pedigree.get_data(sibling2)["y_haplogroup"]
            ):
                return False
        # All full siblings should share MT haplogroups
        return PedigreeReconstructor._check_haplogroups(
            pedigree.get_data(sibling1)["mt_haplogroup"], pedigree.get_data(sibling2)["mt_haplogroup"]
        )

    @staticmethod
    def _check_aunt_uncle_nephew_niece_haplogroups(
        pedigree: Pedigree, aunt_uncle: str, nephew_niece: str, shared_relative_sex: str | None
    ) -> bool:
        """
        Checks if the haplogroups of an aunt/uncle and nephew/niece are compatible.
        """
        if not shared_relative_sex:
            return True

        if (
            shared_relative_sex == "M"
            and pedigree.get_data(aunt_uncle)["sex"] == "M"
            and pedigree.get_data(nephew_niece)["sex"] == "M"
        ):
            return PedigreeReconstructor._check_haplogroups(
                pedigree.get_data(aunt_uncle)["y_haplogroup"], pedigree.get_data(nephew_niece)["y_haplogroup"]
            )
        if shared_relative_sex == "F":
            return PedigreeReconstructor._check_haplogroups(
                pedigree.get_data(aunt_uncle)["mt_haplogroup"], pedigree.get_data(nephew_niece)["mt_haplogroup"]
            )
        return True

    @staticmethod
    def _check_grandparent_grandchild_haplogroups(
        pedigree: Pedigree, grandparent: str, grandchild: str, shared_relative_sex: str | None
    ) -> bool:
        """
        Checks if the haplogroups of a grandparent and grandchild are compatible.
        """
        if not shared_relative_sex:
            return True

        if (
            shared_relative_sex == "M"
            and pedigree.get_data(grandparent)["sex"] == "M"
            and pedigree.get_data(grandchild)["sex"] == "M"
        ):
            return PedigreeReconstructor._check_haplogroups(
                pedigree.get_data(grandparent)["y_haplogroup"], pedigree.get_data(grandchild)["y_haplogroup"]
            )
        if shared_relative_sex == "F" and pedigree.get_data(grandparent)["sex"] == "F":
            return PedigreeReconstructor._check_haplogroups(
                pedigree.get_data(grandparent)["mt_haplogroup"], pedigree.get_data(grandchild)["mt_haplogroup"]
            )
        return True

    @staticmethod
    def _check_half_sibling_haplogroups(
        pedigree: Pedigree, half_sibling1: str, half_sibling2: str, shared_relative_sex: str | None
    ) -> bool:
        """
        Checks if the haplogroups of two half-siblings are compatible.
        """
        if (
            shared_relative_sex == "M"
            and pedigree.get_data(half_sibling1)["sex"] == "M"
            and pedigree.get_data(half_sibling2)["sex"] == "M"
        ):
            return PedigreeReconstructor._check_haplogroups(
                pedigree.get_data(half_sibling1)["y_haplogroup"], pedigree.get_data(half_sibling2)["y_haplogroup"]
            )
        if shared_relative_sex == "F":
            return PedigreeReconstructor._check_haplogroups(
                pedigree.get_data(half_sibling1)["mt_haplogroup"], pedigree.get_data(half_sibling2)["mt_haplogroup"]
            )
        return True

    def _add_relation(self, node1: str, node2: str, degree: str, constraints: str, force_constraints: bool) -> None:
        """
        Connects two nodes in every pedigree.
        """
        assert degree in ["1", "2"]

        new_pedigrees: list[Pedigree] = []
        for pedigree in self._candidate_pedigrees:
            if degree == "1":
                if not force_constraints:
                    new_pedigrees.extend(
                        PedigreeReconstructor._connect_first_degree_relation(
                            pedigree, node1, node2, constraints=self._DEFAULT_CONSTRAINTS["1"]
                        )
                    )
                    new_pedigrees.extend(
                        PedigreeReconstructor._connect_second_degree_relation(
                            pedigree, node1, node2, constraints=self._DEFAULT_CONSTRAINTS["2"]
                        )
                    )
                else:
                    new_pedigrees.extend(
                        PedigreeReconstructor._connect_first_degree_relation(
                            pedigree, node1, node2, constraints=constraints
                        )
                    )

            elif degree == "2":
                if not force_constraints:
                    new_pedigrees.append(pedigree)  # No relation (i.e. false positive)
                    new_pedigrees.extend(
                        PedigreeReconstructor._connect_first_degree_relation(
                            pedigree, node1, node2, constraints=self._DEFAULT_CONSTRAINTS["1"]
                        )
                    )
                    new_pedigrees.extend(
                        PedigreeReconstructor._connect_second_degree_relation(
                            pedigree, node1, node2, constraints=self._DEFAULT_CONSTRAINTS["2"]
                        )
                    )
                else:
                    new_pedigrees.extend(
                        PedigreeReconstructor._connect_second_degree_relation(
                            pedigree, node1, node2, constraints=constraints
                        )
                    )
        self._candidate_pedigrees = new_pedigrees

    @staticmethod
    def _connect_first_degree_relation(pedigree: Pedigree, node1: str, node2: str, constraints: str) -> list[Pedigree]:
        """
        Update pedigree with a first-degree relation.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data

        new_pedigrees: list[Pedigree] = []
        possible_relations: list[str] = constraints.split(";")

        for relation in possible_relations:
            if relation == "parent-child":
                new_pedigrees.extend(PedigreeReconstructor._connect_parent_relation(pedigree, node1, node2))
            if relation == "child-parent":
                new_pedigrees.extend(PedigreeReconstructor._connect_parent_relation(pedigree, node2, node1))
            if relation == "siblings":
                new_pedigrees.extend(PedigreeReconstructor._connect_sibling_relation(pedigree, node1, node2))
        return new_pedigrees

    @staticmethod
    def _connect_second_degree_relation(pedigree: Pedigree, node1: str, node2: str, constraints: str) -> list[Pedigree]:
        """
        Update pedigree with a second-degree relation.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data

        new_pedigrees: list[Pedigree] = []
        possible_relations: list[str] = constraints.split(";")

        for relation in possible_relations:
            if relation == "maternal aunt/uncle-nephew/niece":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_aunt_uncle_relation(pedigree, node1, node2, shared_relative_sex="F")
                )
            if relation == "maternal nephew/niece-aunt/uncle":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_aunt_uncle_relation(pedigree, node2, node1, shared_relative_sex="F")
                )
            if relation == "paternal aunt/uncle-nephew/niece":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_aunt_uncle_relation(pedigree, node1, node2, shared_relative_sex="M")
                )
            if relation == "paternal nephew/niece-aunt/uncle":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_aunt_uncle_relation(pedigree, node2, node1, shared_relative_sex="M")
                )

            if relation == "maternal grandparent-grandchild":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_grandparent_relation(pedigree, node1, node2, shared_relative_sex="F")
                )
            if relation == "maternal grandchild-grandparent":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_grandparent_relation(pedigree, node2, node1, shared_relative_sex="F")
                )
            if relation == "paternal grandparent-grandchild":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_grandparent_relation(pedigree, node1, node2, shared_relative_sex="M")
                )
            if relation == "paternal grandchild-grandparent":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_grandparent_relation(pedigree, node2, node1, shared_relative_sex="M")
                )

            if relation == "maternal half-siblings":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_half_sibling_relation(
                        pedigree, node1, node2, shared_relative_sex="F"
                    )
                )
            if relation == "paternal half-siblings":
                new_pedigrees.extend(
                    PedigreeReconstructor._connect_half_sibling_relation(
                        pedigree, node1, node2, shared_relative_sex="M"
                    )
                )
        return new_pedigrees

    @staticmethod
    def _connect_parent_relation(pedigree: Pedigree, node1: str, node2: str) -> list[Pedigree]:
        """
        Adds a parent-child relation and merges nodes appropriately.
        Returns a list containing the resulting Pedigree, if successful.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data

        # Pre-check invalid relations to avoid unnecessary deep-copying
        if not PedigreeReconstructor._check_parent_child_haplogroups(pedigree, node1, node2):
            return []

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node2)
        original_parent: str
        if new_pedigree.get_data(node1)["sex"] == "M":
            original_parent = new_pedigree.get_father(node2)
        else:
            original_parent = new_pedigree.get_mother(node2)

        if new_pedigree.check_valid_merge(node1, original_parent):
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

        # Pre-check invalid relations to avoid unnecessary deep-copying
        if not PedigreeReconstructor._check_sibling_haplogroups(pedigree, node1, node2):
            return []

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node1)
        new_pedigree.fill_node_parents(node2)

        father1 = new_pedigree.get_father(node1)
        father2 = new_pedigree.get_father(node2)
        if new_pedigree.check_valid_merge(father1, father2):
            new_pedigree.merge_nodes(father1, father2)
            mother1 = new_pedigree.get_mother(node1)
            mother2 = new_pedigree.get_mother(node2)
            if new_pedigree.check_valid_merge(mother1, mother2):
                new_pedigree.merge_nodes(mother1, mother2)
                new_pedigree.add_sibling_relation(node1, node2)
                ret.append(new_pedigree)
        return ret

    @staticmethod
    def _connect_aunt_uncle_relation(
        pedigree: Pedigree, node1: str, node2: str, shared_relative_sex: str | None = None
    ) -> list[Pedigree]:
        """
        Adds an aunt/uncle-nephew/niece relation and merges nodes appropriately.
        Returns a list containing the resulting Pedigree(s), if successful.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data
        assert shared_relative_sex in ["M", "F", None]

        # Pre-check invalid relations to avoid unnecessary deep-copying
        if not PedigreeReconstructor._check_aunt_uncle_nephew_niece_haplogroups(
            pedigree, node1, node2, shared_relative_sex
        ):
            return []

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node2)

        node2_parents: list[str]
        if shared_relative_sex == "M":
            node2_parents = [new_pedigree.get_father(node2)]
        elif shared_relative_sex == "F":
            node2_parents = [new_pedigree.get_mother(node2)]
        else:
            node2_parents = [new_pedigree.get_father(node2), new_pedigree.get_mother(node2)]

        for node2_parent in node2_parents:
            if node1 != node2_parent:
                ret.extend(PedigreeReconstructor._connect_sibling_relation(new_pedigree, node1, node2_parent))
        return ret

    @staticmethod
    def _connect_grandparent_relation(
        pedigree: Pedigree, node1: str, node2: str, shared_relative_sex: str | None = None
    ) -> list[Pedigree]:
        """
        Adds a grandparent-grandchild relation and merges nodes appropriately.
        Returns a list containing the resulting Pedigree(s), if successful.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data
        assert shared_relative_sex in ["M", "F", None]

        # Pre-check invalid relations to avoid unnecessary deep-copying
        if not PedigreeReconstructor._check_grandparent_grandchild_haplogroups(
            pedigree, node1, node2, shared_relative_sex
        ):
            return []

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node2)

        node2_parents: list[str]
        if shared_relative_sex == "M":
            node2_parents = [new_pedigree.get_father(node2)]
        elif shared_relative_sex == "F":
            node2_parents = [new_pedigree.get_mother(node2)]
        else:
            node2_parents = [new_pedigree.get_father(node2), new_pedigree.get_mother(node2)]

        for node2_parent in node2_parents:
            if node1 != node2_parent:
                ret.extend(PedigreeReconstructor._connect_parent_relation(new_pedigree, node1, node2_parent))
        return ret

    @staticmethod
    def _connect_half_sibling_relation(
        pedigree: Pedigree, node1: str, node2: str, shared_relative_sex: str | None = None
    ) -> list[Pedigree]:
        """
        Adds a half-sibling relation and merges nodes appropriately.
        Returns a list containing the resulting Pedigree(s), if successful.
        """
        assert node1 in pedigree.node_to_data and node2 in pedigree.node_to_data

        # Pre-check invalid relations to avoid unnecessary deep-copying
        if not PedigreeReconstructor._check_half_sibling_haplogroups(pedigree, node1, node2, shared_relative_sex):
            return []

        ret: list[Pedigree] = []
        new_pedigree = copy.deepcopy(pedigree)
        new_pedigree.fill_node_parents(node1)
        new_pedigree.fill_node_parents(node2)

        node1_parents: list[str]
        node2_parents: list[str]
        if shared_relative_sex == "M":
            node1_parents = [new_pedigree.get_father(node1)]
            node2_parents = [new_pedigree.get_father(node2)]
        elif shared_relative_sex == "F":
            node1_parents = [new_pedigree.get_mother(node1)]
            node2_parents = [new_pedigree.get_mother(node2)]
        else:
            node1_parents = [new_pedigree.get_father(node1), new_pedigree.get_mother(node1)]
            node2_parents = [new_pedigree.get_father(node2), new_pedigree.get_mother(node2)]

        # Node 1 and Node 2 are half-siblings via one of Node 1's parents
        for node1_parent in node1_parents:
            if node1_parent != node2:
                ret.extend(PedigreeReconstructor._connect_parent_relation(new_pedigree, node1_parent, node2))
        # Node 1 and Node 2 are half-siblings via one of Node 2's parents
        for node2_parent in node2_parents:
            if node2_parent != node1:
                ret.extend(PedigreeReconstructor._connect_parent_relation(new_pedigree, node2_parent, node1))
        return ret

    def _clean_pedigree_data(self) -> None:
        """
        Remove unnecessary entries in Pedigree dicts.
        """
        for pedigree in self._candidate_pedigrees:
            pedigree.clean_data()

        for pedigree in self._final_pedigrees:
            pedigree.clean_data()

    def _validate_pedigree_structures(self) -> None:
        """
        Validate that all candidate pedigrees are consistent.
        """
        for pedigree in self._candidate_pedigrees:
            assert pedigree.validate_structure()

    def _get_pair_to_constraints(self) -> defaultdict[tuple[str, str], list[tuple[str, ...]]]:
        """
        Turn DataFrame of relations/constraints into dict(s) of {node pairs: list of possible relations}.
        Dict values are lists of tuples (as opposed to just tuples)
        because a pair of nodes can share more than 1 relation.
        """
        pair_to_constraints: defaultdict[tuple[str, str], list[tuple[str, ...]]] = defaultdict(list)
        for node1, node2, _, constraints, _ in self._all_relations.itertuples(index=False):
            pair_to_constraints[(node1, node2)].append(tuple(constraints.split(";")))
        for node_pair in pair_to_constraints:
            # Sort by number of constraints so specific constraints are checked first when pruning
            pair_to_constraints[node_pair].sort(key=lambda x: len(x))
        return pair_to_constraints

    def _get_pair_to_relations_so_far(
        self, processed_relations: pd.DataFrame
    ) -> defaultdict[tuple[str, str], list[tuple[str, str, bool]]]:
        """
        Turn DataFrame of relations/constraints processed so far
        into dict(s) of {node pairs: list of (degree, constraints) tuples}.
        """
        pair_to_relations_so_far: defaultdict[tuple[str, str], list[tuple[str, str, bool]]] = defaultdict(list)
        for node1, node2, degree, constraints, force_constraints in processed_relations.itertuples(index=False):
            pair_to_relations_so_far[(node1, node2)].append((degree, constraints, force_constraints))
        return pair_to_relations_so_far

    def _prune_pedigrees(
        self,
        pair_to_relations_so_far: defaultdict[tuple[str, str], list[tuple[str, str, bool]]],
        check_half_siblings: bool,
    ) -> None:
        """
        Remove pedigrees with inconsistencies.
        """
        seen_topologies = set()
        new_potential_pedigrees = []
        for pedigree in self._candidate_pedigrees:
            if (
                pedigree.validate_members(set(self._node_data["id"]))
                and pedigree.validate_can_have_children()
                and pedigree.validate_inbreeding()
                and pedigree.validate_years_before_present()
                and pedigree.validate_forced_constraints(pair_to_relations_so_far)
            ):
                pedigree.update_haplogroups()
                if pedigree.validate_haplogroups():
                    topology = pedigree.get_topo_sort()
                    if topology not in seen_topologies:
                        new_potential_pedigrees.append(pedigree)
                        seen_topologies.add(topology)
        # Shuffle to avoid ordering bias in epsilon-greedy sampling
        random.shuffle(new_potential_pedigrees)

        strikes = []
        third_degree_strikes = []
        counts: defaultdict[int, int] = defaultdict(int)
        for pedigree in new_potential_pedigrees:
            num_strikes, _ = pedigree.count_inconsistencies(
                self._pair_to_constraints, pair_to_relations_so_far, check_half_siblings
            )
            num_third_degree_strikes = pedigree.count_third_degree_inconsistencies(self._pair_to_constraints)
            strikes.append(num_strikes)
            third_degree_strikes.append(num_third_degree_strikes)
            counts[num_strikes] += 1
        logger.info(f"Strike counts before pruning: {str(dict(sorted(counts.items())))}")

        def epsilon_greedy_sample(
            pedigrees: list[Pedigree],
            strikes: list[int],
            third_degree_strikes: list[int],
            epsilon: float,
            max_candidate_pedigrees: int,
        ) -> list[Pedigree]:
            assert len(pedigrees) == len(strikes)
            if len(pedigrees) <= max_candidate_pedigrees:
                return pedigrees

            sorted_pedigrees = [
                pedigree
                for pedigree, _, _ in sorted(
                    zip(pedigrees, strikes, third_degree_strikes, strict=True), key=lambda x: (x[1], x[2])
                )
            ]
            exploitation_max_candidate_pedigrees = int((1 - epsilon) * max_candidate_pedigrees)
            exploration_max_candidate_pedigrees = max_candidate_pedigrees - exploitation_max_candidate_pedigrees

            exploitation_pedigrees = sorted_pedigrees[:exploitation_max_candidate_pedigrees]
            exploration_pedigrees = random.sample(
                sorted_pedigrees[exploitation_max_candidate_pedigrees:], exploration_max_candidate_pedigrees
            )
            return exploitation_pedigrees + exploration_pedigrees

        num_processed_relations = sum(len(relations) for relations in pair_to_relations_so_far.values())
        if num_processed_relations < len(self._first_and_second_degree_relations):
            self._candidate_pedigrees = epsilon_greedy_sample(
                new_potential_pedigrees,
                strikes,
                third_degree_strikes,
                epsilon=self._epsilon,
                max_candidate_pedigrees=self._max_candidate_pedigrees,
            )
        else:
            # Final iteration
            best_pedigrees = [
                pedigree
                for pedigree, num_strikes in zip(new_potential_pedigrees, strikes, strict=True)
                if num_strikes == min(strikes)
            ]
            # Use 3rd-degree strikes as tiebreaker
            third_degree_strikes = [
                pedigree.count_third_degree_inconsistencies(self._pair_to_constraints) for pedigree in best_pedigrees
            ]

            self._final_pedigrees.extend(
                [
                    pedigree
                    for pedigree, num_strikes in zip(best_pedigrees, third_degree_strikes, strict=True)
                    if num_strikes == min(third_degree_strikes)
                ]
            )
            self._final_strike_counts = []
            self._final_strike_logs = []
            for pedigree in self._final_pedigrees:
                strike_count, strike_log = pedigree.count_inconsistencies(
                    self._pair_to_constraints, pair_to_relations_so_far, check_half_siblings=True
                )
                self._final_strike_counts.append(strike_count)
                self._final_strike_logs.append(strike_log)

    def _write_corrected_input_relations(
        self, strike_count: int, strike_log: list[tuple[str, str, str, str]], path: str
    ) -> None:
        """
        Write corrected input relations to file. Includes information about added/removed/changed input relations.
        """
        added_relations = []
        removed_relations = []
        for node1, node2, degree, constraints in strike_log:
            if degree[0] == "+":
                added_relations.append((node1, node2, degree[1], constraints))
            else:
                removed_relations.append((node1, node2, degree[1], constraints))
        removed_relations_set = set(removed_relations)

        # Separate out *changed* relations (added relation + removed relation pair, e.g., 1st-degree -> 2nd-degree)
        changed_node_pairs = set()
        for add_node1, add_node2, _, _ in added_relations:
            for remove_node1, remove_node2, _, _ in removed_relations:
                if (add_node1 == remove_node1 and add_node2 == remove_node2) or (
                    add_node2 == remove_node1 and add_node1 == remove_node2
                ):
                    changed_node_pairs.add((add_node1, add_node2))

        with open(path, "w") as file:
            file.write("id1,id2,degree,constraints\n")  # Header line
            file.write(f"# Final inconsistency count: {strike_count}\n")

            def write_relations_line(node1, node2, degree, constraints, commented=False):
                if constraints == self._DEFAULT_CONSTRAINTS[degree]:
                    # Don't write default constraints to file
                    constraints = ""
                if commented:
                    file.write("# ")
                file.write(f"{node1},{node2},{degree},{constraints}\n")

            file.write("# Added relations\n")
            # Sort for consistency
            for node1, node2, degree, constraints in sorted(added_relations):
                if (node1, node2) not in changed_node_pairs and (node2, node1) not in changed_node_pairs:
                    write_relations_line(node1, node2, degree, constraints)

            file.write("\n# Removed relations\n")
            for node1, node2, degree, constraints in sorted(removed_relations):
                if (node1, node2) not in changed_node_pairs and (node2, node1) not in changed_node_pairs:
                    write_relations_line(node1, node2, degree, constraints, commented=True)

            file.write("\n# Changed relations\n")
            # Pair up changed relations (add + remove)
            for node1, node2 in sorted(changed_node_pairs):
                # We want to write the two nodes in the correct (original) order
                node1_to_write = None
                node2_to_write = None
                for node1_remove, node2_remove, degree_remove, constraints_remove in removed_relations:
                    if (node1_remove, node2_remove) == (node1, node2) or (node2_remove, node1_remove) == (node1, node2):
                        write_relations_line(
                            node1_remove, node2_remove, degree_remove, constraints_remove, commented=True
                        )
                        # The removed nodes follow the original input order
                        node1_to_write = node1_remove
                        node2_to_write = node2_remove
                for node1_add, node2_add, degree_add, constraints_add in added_relations:
                    if (node1_add, node2_add) == (node1, node2) or (node2_add, node1_add) == (node1, node2):
                        assert node1_to_write and node2_to_write
                        write_relations_line(node1_to_write, node2_to_write, degree_add, constraints_add)

            file.write("\n# Unchanged relations\n")
            for node1, node2, degree, constraints, _ in self._all_relations.itertuples(index=False):
                if (node1, node2, degree, constraints) not in removed_relations_set:
                    assert (node2, node1, degree, constraints) not in removed_relations_set
                    write_relations_line(node1, node2, degree, constraints)
