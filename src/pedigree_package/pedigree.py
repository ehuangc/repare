import math
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict, deque


class Pedigree:
    """
    Describes a pedigree configuration for a set of nodes.
    """
    def __init__(self) -> None:
        self.num_placeholders = 0
        self.node_to_data: dict[str, dict[str, str | bool | float]] = dict()
        self.node_to_father = defaultdict(str)
        self.node_to_mother = defaultdict(str)
        self.node_to_children = defaultdict(set)
        self.node_to_siblings = defaultdict(set)

    def __deepcopy__(self, memo: dict) -> "Pedigree":
        """
        Custom (faster) deepcopy implementation.
        """
        cls = self.__class__
        new_pedigree = cls.__new__(cls)
        memo[id(self)] = new_pedigree
        
        new_pedigree.num_placeholders = self.num_placeholders
        new_pedigree.node_to_data = dict()
        for k, v in self.node_to_data.items():
            new_pedigree.node_to_data[k] = v.copy()

        new_pedigree.node_to_father = self.node_to_father.copy()
        new_pedigree.node_to_mother = self.node_to_mother.copy()
        new_pedigree.node_to_children = defaultdict(set)
        for k, v in self.node_to_children.items():
            new_pedigree.node_to_children[k] = v.copy()
        new_pedigree.node_to_siblings = defaultdict(set)
        for k, v in self.node_to_siblings.items():
            new_pedigree.node_to_siblings[k] = v.copy()
        return new_pedigree

    def get_topo_sort(self) -> tuple[str, ...]:
        """
        Gets pedigree topological sort of the Pedigree. See https://doi.org/10.1089/cmb.2011.0254.
        """
        leaf_nodes: list[str] = sorted([node for node in self.node_to_data if not self.node_to_children[node]])
        result: list[str] = []
        # DFS
        def visit(node: str) -> None:
            if not self.node_to_father[node]:  # No father == no mother
                assert not self.node_to_mother[node]
                result.append(node)
            else:
                # Visit father first
                visit(self.node_to_father[node])
                visit(self.node_to_mother[node])
                result.append(node)

        for node in leaf_nodes:
            visit(node)
            result.append("")  # Break between each node's path

        # Re-label placeholder nodes
        placeholder_to_idx: dict[str, int] = {}
        for i, node in enumerate(result):
            if node.isnumeric():
                if node not in placeholder_to_idx:
                    placeholder_to_idx[node] = len(placeholder_to_idx)
                mt_haplogroup = self.node_to_data[node]["mt_haplogroup"]
                if self.node_to_data[node]["sex"] == "M":
                    y_haplogroup = self.node_to_data[node]["y_haplogroup"]
                    result[i] = f"M {placeholder_to_idx[node]} {mt_haplogroup} {y_haplogroup}" # Unique identifier for placeholder
                else:
                    result[i] = f"F {placeholder_to_idx[node]} {mt_haplogroup}"
        return tuple(result)

    def add_node(self, node_id: str, sex: str, y_haplogroup: str, mt_haplogroup: str, can_have_children: bool, can_be_inbred: bool, years_before_present: float) -> None:
        """
        Add a node to the pedigree. If haplogroup is unknown, set argument to empty string ("").
        """
        self.node_to_data[node_id]: dict[str, str | bool | float] = dict()
        self.node_to_data[node_id]["sex"] = sex
        if y_haplogroup and sex == "F":
            raise ValueError("Only males can have y_haplogroup values.")
        self.node_to_data[node_id]["y_haplogroup"] = y_haplogroup
        self.node_to_data[node_id]["mt_haplogroup"] = mt_haplogroup
        self.node_to_data[node_id]["can_have_children"] = can_have_children
        self.node_to_data[node_id]["can_be_inbred"] = can_be_inbred
        self.node_to_data[node_id]["years_before_present"] = years_before_present

    def add_parent_relation(self, node1: str, node2: str) -> None:
        """
        Adds a parent-child relationship to the tree.
        Node 1 is the parent and Node 2 is the child.
        Note: Overwrites existing parent, does not merge.
        """
        assert node1 != node2
        assert node1 in self.node_to_data and node2 in self.node_to_data

        if self.node_to_data[node1]["sex"] == "M":
            self.node_to_father[node2] = node1
        else:
            self.node_to_mother[node2] = node1
        self.node_to_children[node1].add(node2)

        # Add any sibling relations that are created
        node2_parents = [self.node_to_father[node2], self.node_to_mother[node2]]
        if node2_parents[0] and node2_parents[1]:  # Make sure Node 2 has 2 known parents
            for node1_child in self.node_to_children[node1]:
                if node1_child != node2 and [self.node_to_father[node1_child], self.node_to_mother[node1_child]] == node2_parents:
                    self.add_sibling_relation(node1_child, node2)

    def add_sibling_relation(self, node1: str, node2: str) -> None:
        """
        Adds a sibling relationship to the tree.
        Note: Does not merge parents.
        """
        assert node1 != node2
        assert node1 in self.node_to_data and node2 in self.node_to_data

        self.node_to_siblings[node1].add(node2)
        self.node_to_siblings[node2].add(node1)

        # Update siblings of siblings
        for node1_sibling in self.node_to_siblings[node1]:
            if node1_sibling != node2:
                self.node_to_siblings[node1_sibling].add(node2)
                self.node_to_siblings[node2].add(node1_sibling)
        for node2_sibling in self.node_to_siblings[node2]:
            if node1 != node2_sibling:
                self.node_to_siblings[node1].add(node2_sibling)
                self.node_to_siblings[node2_sibling].add(node1)

    def merge_nodes(self, node1: str, node2: str) -> None:
        """
        Merge the two nodes as if they were one person. Note this involves merging the nodes' ancestors.
        When merging nodes, if one is named and one is a placeholder, keep the named one's name.
        If both are placeholders (doesn't matter) or both are named (will be pruned), just keep any name.
        """
        assert node1 in self.node_to_data and node2 in self.node_to_data

        pair_queue: deque[tuple[str, str]] = deque([(node1, node2)])  # Pairs of nodes to merge
        while pair_queue:
            node1, node2 = pair_queue.popleft()
            if node1 != node2:
                name_to_keep = node1 if not node1.isnumeric() else node2
                name_to_discard = node2 if name_to_keep == node1 else node1

                # Update relations for relatives of the discarded node
                name_to_discard_father = self.node_to_father[name_to_discard]
                name_to_discard_mother = self.node_to_mother[name_to_discard]
                if name_to_discard_father:
                    self.node_to_children[name_to_discard_father].remove(name_to_discard)
                if name_to_discard_mother:
                    self.node_to_children[name_to_discard_mother].remove(name_to_discard)

                for child in self.node_to_children[name_to_discard]:
                    if name_to_keep == child:  # Merging a parent and child - we will see this when there is inbreeding
                        if self.node_to_data[name_to_keep]["sex"] == "M":
                            del self.node_to_father[name_to_keep]
                        else:
                            del self.node_to_mother[name_to_keep]
                    else:
                        self.add_parent_relation(name_to_keep, child)  # This also handles having the correct name of the merged parents from last loop iteration
                
                # Remove all occurrences of name_to_discard in its sibling's sibling sets first so that self.add_sibling_relation() does not add it back in
                for sibling in self.node_to_siblings[name_to_discard]:
                    self.node_to_siblings[sibling].remove(name_to_discard)
                for sibling in self.node_to_siblings[name_to_discard]:
                    if sibling != name_to_keep:
                        self.add_sibling_relation(sibling, name_to_keep)

                # Recursively merge parent relations of Node 1 and Node 2
                father1 = self.node_to_father[name_to_keep]
                father2 = self.node_to_father[name_to_discard]
                mother1 = self.node_to_mother[name_to_keep]
                mother2 = self.node_to_mother[name_to_discard]
                if father1 and father2:
                    pair_queue.append((father1, father2))
                elif father2:
                    if father2 == name_to_keep:
                        del self.node_to_father[name_to_keep]
                    else:
                        self.add_parent_relation(father2, name_to_keep)  # Switch name_to_keep's father to name_to_discard's father
                if mother1 and mother2:
                    pair_queue.append((mother1, mother2))
                elif mother2:
                    if mother2 == name_to_keep:
                        del self.node_to_mother[name_to_keep]
                    else:
                        # Switch name_to_keep's mother to name_to_discard's mother
                        self.add_parent_relation(mother2, name_to_keep)

                # Update any nodes in the queue whose names might have been changed
                for idx, (node1, node2) in enumerate(pair_queue):
                    if node1 == name_to_discard and node2 == name_to_discard:
                        pair_queue[idx] = (name_to_keep, name_to_keep)
                    elif node1 == name_to_discard:
                        pair_queue[idx] = (name_to_keep, node2)
                    elif node2 == name_to_discard:
                        pair_queue[idx] = (node1, name_to_keep)

                del self.node_to_data[name_to_discard]
                del self.node_to_father[name_to_discard]
                del self.node_to_mother[name_to_discard]
                del self.node_to_children[name_to_discard]
                del self.node_to_siblings[name_to_discard]

    def check_cycles_if_merged(self, node1: str, node2: str) -> bool:
        """
        Returns True if merging Node 1 and Node 2 (and their ancestors) would result in a cycle.
        """
        # Get sets of nodes that would be merged if Node 1 and Node 2 were merged (i.e. ancestors of Node 1 and Node 2)
        # Note that we get sets and not just pairs because of potential inbreeding (e.g. one node is both a parent and grandparent of another node)
        merge_sets: list[set[str]] = []
        merge_queue = deque([(node1, node2)])
        while merge_queue:
            curr_node1, curr_node2 = merge_queue.popleft()
            if curr_node1 != curr_node2:  # Update merge sets
                updated = False
                for merge_set in merge_sets:
                    if curr_node1 in merge_set or curr_node2 in merge_set:
                        merge_set.update([curr_node1, curr_node2])
                        updated = True
                        break
                if not updated:
                    merge_sets.append({curr_node1, curr_node2})
                # Add parents to the queue
                curr_father1 = self.node_to_father[curr_node1]
                curr_father2 = self.node_to_father[curr_node2]
                curr_mother1 = self.node_to_mother[curr_node1]
                curr_mother2 = self.node_to_mother[curr_node2]
                if curr_father1 and curr_father2:
                    merge_queue.append((curr_father1, curr_father2))
                if curr_mother1 and curr_mother2:
                    merge_queue.append((curr_mother1, curr_mother2))

        # DFS cycle detection
        def dfs(node):
            merged_nodes: set[str] | None = None
            for merge_set in merge_sets:
                if node in merge_set:
                    merged_nodes = merge_set
                    break

            if merged_nodes:
                if node in in_progress:
                    return True
                if node in explored:
                    return False
                in_progress.update(merged_nodes)
                for child in [child for node in merged_nodes for child in self.node_to_children[node]]:
                    if dfs(child):
                        return True
                in_progress.difference_update(merged_nodes)
                explored.update(merged_nodes)
            else:
                if node in in_progress:
                    return True
                if node in explored:
                    return False
                in_progress.add(node)
                for child in self.node_to_children[node]:
                    if dfs(child):
                        return True
                in_progress.remove(node)
                explored.add(node)
            return False

        explored: set[str] = set()
        in_progress: set[str] = set()
        # Check for cycles starting from each node
        for node in self.node_to_data:
            if dfs(node):
                return True  # Cycle detected
        return False

    def fill_node_parents(self, node: str) -> None:
        """
        If the given node doesn't have parents, add placeholder parents.
        If it does, do nothing.
        """
        assert node in self.node_to_data

        father = self.node_to_father[node]
        mother = self.node_to_mother[node]

        if not father:
            father_id = str(self.num_placeholders)
            self.add_node(node_id=father_id, sex="M", y_haplogroup="*", mt_haplogroup="*", can_have_children=True, can_be_inbred=True, years_before_present=math.nan)

            self.add_parent_relation(father_id, node)
            for sibling in self.node_to_siblings[node]:
                self.add_parent_relation(father_id, sibling)
            self.num_placeholders += 1

        if not mother:
            mother_id = str(self.num_placeholders)
            self.add_node(node_id=mother_id, sex="F", y_haplogroup="", mt_haplogroup="*", can_have_children=True, can_be_inbred=True, years_before_present=math.nan)
            
            self.add_parent_relation(mother_id, node)
            for sibling in self.node_to_siblings[node]:
                self.add_parent_relation(mother_id, sibling)
            self.num_placeholders += 1

    def update_haplogroups(self) -> None:
        """
        Update haplogroups of placeholder nodes.
        """
        for node in self.node_to_data:
            y_haplogroup: str = self.node_to_data[node]["y_haplogroup"]
            y_lineage: deque[str] = deque([self.node_to_father[node]] + [child for child in self.node_to_children[node] if self.node_to_data[child]["sex"] == "M"])
            
            while y_lineage:
                curr_node = y_lineage.popleft()
                if not curr_node or "*" not in self.node_to_data[curr_node]["y_haplogroup"] or self.node_to_data[curr_node]["y_haplogroup"].rstrip("*") == y_haplogroup.rstrip("*"):
                    continue
                # Overwrite/extend Y haplogroup if it contains a "*" and is a strict subset of the "leaf" haplogroup
                if y_haplogroup.startswith(self.node_to_data[curr_node]["y_haplogroup"].rstrip("*")):
                    self.node_to_data[curr_node]["y_haplogroup"] = y_haplogroup if y_haplogroup[-1] == "*" else y_haplogroup + "*"
                    y_lineage.append(self.node_to_father[curr_node])
                    for curr_node_child in self.node_to_children[curr_node]:
                        # Only males have Y chromosome
                        if self.node_to_data[curr_node_child]["sex"] == "M":
                            y_lineage.append(curr_node_child)

            mt_haplogroup: str = self.node_to_data[node]["mt_haplogroup"]
            mt_lineage: deque[str] = deque([self.node_to_mother[node]])
            # Only females pass on mitochondrial DNA to children
            if self.node_to_data[node]["sex"] == "F":
                mt_lineage.extend(self.node_to_children[node])

            while mt_lineage:
                curr_node = mt_lineage.popleft()
                if not curr_node or "*" not in self.node_to_data[curr_node]["mt_haplogroup"] or self.node_to_data[curr_node]["mt_haplogroup"].rstrip("*") == mt_haplogroup.rstrip("*"):
                    continue
                # Overwrite/extend mitochondrial haplogroup if it contains a "*" and is a strict subset of the "leaf" haplogroup
                if mt_haplogroup.startswith(self.node_to_data[curr_node]["mt_haplogroup"].rstrip("*")):
                    self.node_to_data[curr_node]["mt_haplogroup"] = mt_haplogroup if mt_haplogroup[-1] == "*" else mt_haplogroup + "*"
                    mt_lineage.append(self.node_to_mother[curr_node])
                    if self.node_to_data[curr_node]["sex"] == "F":
                        mt_lineage.extend(self.node_to_children[curr_node])

    def validate_members(self, members: set) -> bool:
        """
        Validates this tree based on the member nodes it should contain.
        """
        non_placeholder_nodes = self.get_non_placeholder_nodes()
        # Return False if pedigree doesn't have all the nodes it's supposed to (because of invalid merging)
        return non_placeholder_nodes == members

    def validate_haplogroups(self) -> bool:
        """
        Validates that all haplogroups are consistent.
        """
        def haplogroups_agree(haplogroup1: str, haplogroup2: str) -> bool:
            if "*" not in haplogroup1 and "*" not in haplogroup2:
                return haplogroup1 == haplogroup2
            elif "*" in haplogroup1 and "*" in haplogroup2:
                return haplogroup1.startswith(haplogroup2.rstrip("*")) or haplogroup2.startswith(haplogroup1.rstrip("*"))
            elif "*" in haplogroup1:
                return haplogroup2.startswith(haplogroup1.rstrip("*"))
            else:
                return haplogroup1.startswith(haplogroup2.rstrip("*"))
            
        for parent, child in self.get_parent_child_pairs():
            if self.node_to_data[parent]["sex"] == "F":
                if not haplogroups_agree(self.node_to_data[parent]["mt_haplogroup"], self.node_to_data[child]["mt_haplogroup"]):
                    return False
            elif self.node_to_data[parent]["sex"] == "M" and self.node_to_data[child]["sex"] == "M":
                if not haplogroups_agree(self.node_to_data[parent]["y_haplogroup"], self.node_to_data[child]["y_haplogroup"]):
                    return False
        return True

    def validate_can_have_children(self) -> bool:
        """
        Validates that nodes that cannot have children do not have children.
        """
        for node in self.get_non_placeholder_nodes():
            if len(self.node_to_children[node]) > 0 and not self.node_to_data[node]["can_have_children"]:
                return False
        return True

    def validate_inbreeding(self) -> bool:
        """
        Validates that nodes that are known to be not inbred are not inbred.
        """
        related_pairs = self.get_related_pairs()
        for node in self.get_non_placeholder_nodes():
            if not self.node_to_data[node]["can_be_inbred"]:
                father = self.node_to_father[node]
                mother = self.node_to_mother[node]
                if (father, mother) in related_pairs or (mother, father) in related_pairs:
                    return False
        return True

    def validate_years_before_present(self) -> bool:
        """
        Validates that nodes do not postdate their descendants.
        """
        leaf_nodes: list[str] = [node for node in self.node_to_data if not self.node_to_children[node]]
        # DFS
        def visit(node: str, curr_years_before_present: float) -> None:
            years_before_present = self.node_to_data[node]["years_before_present"]
            if not math.isnan(years_before_present):
                if years_before_present < curr_years_before_present:  # Node postdates its descendants
                    return False
                else:
                    curr_years_before_present = years_before_present

            if self.node_to_father[node]:
                assert self.node_to_mother[node]
                if not visit(self.node_to_father[node], curr_years_before_present):
                    return False
                if not visit(self.node_to_mother[node], curr_years_before_present):
                    return False
            return True

        for node in leaf_nodes:
            if not visit(node, float("-inf")):
                return False
        return True

    def validate_forced_constraints(self, pair_to_relations_so_far: defaultdict[tuple[str, str], list[tuple[str, str, bool]]]) -> bool:
        """
        Validates that forced constraints so far are present in the pedigree.
        Note: Additional relations between two nodes are allowed as long as the forced constraints are present.
        """
        for (node1, node2), degree_constraints in pair_to_relations_so_far.items():
            for degree, constraints, force_constraints in degree_constraints:
                if force_constraints and not self.is_relation_in_pedigree(node1, node2, constraints.split(";")):
                    return False
        return True

    def count_inconsistencies(self, pair_to_constraints: defaultdict[tuple[str, str], list[tuple[str, ...]]], pair_to_relations_so_far: defaultdict[tuple[str, str], list[tuple[str, str, bool]]], check_half_siblings: bool) -> tuple[int, list[tuple[str, str, str]]]:
        """
        Validates this tree based on the input relation data.
        If check_half_siblings is False, don't check for extraneous half-sibling relations because the 2 non-shared parents might be merged later.
        Returns count of inconsistencies with the input data as well as a log of the inconsistencies.
        Note: pair_to_constraints values must be sorted by increasing length so that specific constraints are checked first.
        """
        for (node1, node2) in pair_to_constraints:
            assert (node2, node1) not in pair_to_constraints  # Ensure no duplicate/symmetric entries
        pair_to_constraints_seen_entries: defaultdict[tuple[str, str], set[int]] = defaultdict(set)  # Marks which entries in pair_to_constraints have been seen/used

        def is_relation_in_input_data(node1: str, node2: str, relation: str) -> bool:
            if (node1, node2) in pair_to_constraints:
                for idx, constraints in enumerate(pair_to_constraints[(node1, node2)]):
                    if idx not in pair_to_constraints_seen_entries[(node1, node2)] and relation in constraints:
                        return True
            return False

        def remove_relation_from_input_data(node1: str, node2: str, relation: str) -> None:
            if (node1, node2) in pair_to_constraints:
                for idx, constraints in enumerate(pair_to_constraints[(node1, node2)]):
                    if relation in constraints:
                        pair_to_constraints_seen_entries[(node1, node2)].add(idx)
                        break

        def validate_relation(node1: str, node2: str, relation: str, strike_log: list[str, str, str, str]) -> None:
            relation_to_degree = {
                "parent-child": "1", "child-parent": "1", "siblings": "1",
                "maternal aunt/uncle-nephew/niece": "2", "maternal nephew/niece-aunt/uncle": "2",
                "paternal aunt/uncle-nephew/niece": "2", "paternal nephew/niece-aunt/uncle": "2",
                "maternal grandparent-grandchild": "2", "maternal grandchild-grandparent": "2", 
                "paternal grandparent-grandchild": "2", "paternal grandchild-grandparent": "2",
                "maternal half-siblings": "2", "paternal half-siblings": "2"
            }
            flipped_relations = {
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
                "paternal half-siblings": "paternal half-siblings"  # Symmetric
            }
            if not is_relation_in_input_data(node1, node2, relation) and not is_relation_in_input_data(node2, node1, flipped_relations[relation]):
                strike_log.append((node1, node2, f"+{relation_to_degree[relation]}", ""))
            remove_relation_from_input_data(node1, node2, relation)
            remove_relation_from_input_data(node2, node1, flipped_relations[relation])
        
        strike_log: list[str, str, str, str] = []  # (node1, node2, +/- relation degree, constraints)
        # Check that relations in the pedigree are present in the input data
        for parent, child in self.get_parent_child_pairs(include_placeholders=False):
            validate_relation(parent, child, "parent-child", strike_log)
        for sibling1, sibling2 in self.get_sibling_pairs(include_placeholders=False):
            validate_relation(sibling1, sibling2, "siblings", strike_log)

        for aunt_uncle, nephew_niece in self.get_aunt_uncle_nephew_niece_pairs(include_placeholders=False, shared_relative_sex="F"):
            validate_relation(aunt_uncle, nephew_niece, "maternal aunt/uncle-nephew/niece", strike_log)
        for aunt_uncle, nephew_niece in self.get_aunt_uncle_nephew_niece_pairs(include_placeholders=False, shared_relative_sex="M"):
            validate_relation(aunt_uncle, nephew_niece, "paternal aunt/uncle-nephew/niece", strike_log)

        for grandparent, grandchild in self.get_grandparent_grandchild_pairs(include_placeholders=False, shared_relative_sex="F"):
            validate_relation(grandparent, grandchild, "maternal grandparent-grandchild", strike_log)
        for grandparent, grandchild in self.get_grandparent_grandchild_pairs(include_placeholders=False, shared_relative_sex="M"):
            validate_relation(grandparent, grandchild, "paternal grandparent-grandchild", strike_log)

        if check_half_siblings:
            for half_sibling1, half_sibling2 in self.get_half_sibling_pairs(include_placeholders=False, shared_relative_sex="F"):
                validate_relation(half_sibling1, half_sibling2, "maternal half-siblings", strike_log)
            for half_sibling1, half_sibling2 in self.get_half_sibling_pairs(include_placeholders=False, shared_relative_sex="M"):
                validate_relation(half_sibling1, half_sibling2, "paternal half-siblings", strike_log)

        # Check for "dropped" input relations
        # Note: We use constrained relations instead of all relations because we want to catch half-siblings that explicitly should be some other relation even when check_half_siblings is False
        # The purpose of check_half_siblings is to avoid marking *incidental* half-siblings, not half-siblings that should be something else
        for (node1, node2), degrees_constraints in pair_to_relations_so_far.items():
            if len(degrees_constraints) == 1:  # If only one input relation between these two nodes, simple check is much faster
                degree, constraints, _ = degrees_constraints[0]
                if not self.is_relation_in_pedigree(node1, node2, constraints.split(";")):
                    strike_log.append((node1, node2, f"-{degree}", constraints))
            else:
                pedigree_shared_relations: defaultdict(int) = self.get_relations_between_nodes(node1, node2, include_maternal_paternal=True)
                for degree, constraints, _ in degrees_constraints:
                    present_flag = False
                    for constraint in constraints.split(";"):
                        if constraint in pedigree_shared_relations:
                            present_flag = True
                            pedigree_shared_relations[constraint] -= 1
                            if pedigree_shared_relations[constraint] == 0:
                                del pedigree_shared_relations[constraint]
                            break
                    if not present_flag:
                        strike_log.append((node1, node2, f"-{degree}", constraints))

        # Count # of strikes (will not equal len(strike_log) because we don't want to double-count *changed* relations)
        strike_count: int = 0
        node_pair_strike_balances: defaultdict[tuple[str, str], int] = defaultdict(int)
        for node1, node2, strike, _ in strike_log:
            if strike[0] == "+":
                if node_pair_strike_balances[(node1, node2)] >= 0:
                    strike_count += 1
                node_pair_strike_balances[(node1, node2)] += 1
                node_pair_strike_balances[(node2, node1)] += 1

            elif strike[0] == "-":
                if node_pair_strike_balances[(node1, node2)] <= 0:
                    strike_count += 1
                node_pair_strike_balances[(node1, node2)] -= 1
                node_pair_strike_balances[(node2, node1)] -= 1
        return strike_count, strike_log

    def count_third_degree_inconcistencies(self, pair_to_constraints: defaultdict[tuple[str, str], list[tuple[str, ...]]]) -> int:
        """
        Counts only one-sided inconsistencies in third-degree relations. 
        Used as a "tie-breaker" for 1st- and 2nd-degree inconsistences.
        """
        for (node1, node2) in pair_to_constraints:
            assert (node2, node1) not in pair_to_constraints  # Ensure no duplicate/symmetric entries
        pair_to_constraints_seen_entries: defaultdict[tuple[str, str], set[int]] = defaultdict(set)  # Marks which entries in pair_to_constraints have been seen/used

        def is_relation_in_input_data(node1: str, node2: str, relation: str) -> bool:
            if (node1, node2) in pair_to_constraints:
                for idx, constraints in enumerate(pair_to_constraints[(node1, node2)]):
                    if idx not in pair_to_constraints_seen_entries[(node1, node2)] and relation in constraints:
                        return True
            return False

        def remove_relation_from_input_data(node1: str, node2: str, relation: str) -> None:
            if (node1, node2) in pair_to_constraints:
                for idx, constraints in enumerate(pair_to_constraints[(node1, node2)]):
                    if relation in constraints:
                        pair_to_constraints_seen_entries[(node1, node2)].add(idx)
                        break
        
        def validate_relation(node1: str, node2: str, relation: str) -> bool:
            flipped_relations = {
                "half aunt/uncle-half nephew/niece": "half nephew/niece-half aunt/uncle",
                "half nephew/niece-half aunt/uncle": "half aunt/uncle-half nephew/niece",
                "greatgrandparent-greatgrandchild": "greatgrandchild-greatgrandparent",
                "greatgrandchild-greatgrandparent": "greatgrandparent-greatgrandchild",
                "grandaunt/granduncle-grandnephew/grandniece": "grandnephew/grandniece-grandaunt/granduncle",
                "grandnephew/grandniece-grandaunt/granduncle": "grandaunt/granduncle-grandnephew/grandniece",
                "first cousins": "first cousins"  # Symmetric
            }
            ret = False
            if not is_relation_in_input_data(node1, node2, relation) and not is_relation_in_input_data(node2, node1, flipped_relations[relation]):
                ret = True
            remove_relation_from_input_data(node1, node2, relation)
            remove_relation_from_input_data(node2, node1, flipped_relations[relation])
            return ret
            
        strike_count: int = 0
        for half_aunt_uncle, half_nephew_niece in self.get_half_aunt_uncle_nephew_niece_pairs(include_placeholders=False):
            strike_count += validate_relation(half_aunt_uncle, half_nephew_niece, "half aunt/uncle-half nephew/niece")
        for greatgrandparent, greatgrandchild in self.get_greatgrandparent_greatgrandchild_pairs(include_placeholders=False):
            strike_count += validate_relation(greatgrandparent, greatgrandchild, "greatgrandparent-greatgrandchild")
        for grandaunt_granduncle, grandnephew_grandniece in self.get_grandaunt_granduncle_grandnephew_grandniece_pairs(include_placeholders=False):
            strike_count += validate_relation(grandaunt_granduncle, grandnephew_grandniece, "grandaunt/granduncle-grandnephew/grandniece")
        for first_cousin1, first_cousin2 in self.get_first_cousin_pairs(include_placeholders=False):
            strike_count += validate_relation(first_cousin1, first_cousin2, "first cousins")
        return strike_count

    def is_relation_in_pedigree(self, node1: str, node2: str, relations_list: list[str]) -> bool:
        """
        Returns True if *any* of the relations in relations_list are present between node1 and node2 in the pedigree.
        """
        assert node1 in self.node_to_data and node2 in self.node_to_data

        for relation in relations_list:
            if relation == "parent-child":
                if node2 in self.node_to_children[node1]:
                    return True
            if relation == "child-parent":
                if node1 in self.node_to_children[node2]:
                    return True
            if relation == "siblings":
                if node2 in self.node_to_siblings[node1]:
                    assert node1 in self.node_to_siblings[node2]
                    return True

            if relation == "aunt/uncle-nephew/niece":
                for sibling in self.node_to_siblings[node1]:
                    if node2 in self.node_to_children[sibling]:
                        return True
            if relation == "nephew/niece-aunt/uncle":
                for sibling in self.node_to_siblings[node2]:
                    if node1 in self.node_to_children[sibling]:
                        return True
            if relation == "grandparent-grandchild":
                for child in self.node_to_children[node1]:
                    if node2 in self.node_to_children[child]:
                        return True
            if relation == "grandchild-grandparent":
                for child in self.node_to_children[node2]:
                    if node1 in self.node_to_children[child]:
                        return True
            if relation == "half-siblings":
                if self.node_to_father[node2]:
                    if node1 in self.node_to_children[self.node_to_father[node2]] and self.node_to_mother[node1] != self.node_to_mother[node2]:
                        return True
                if self.node_to_mother[node2]:
                    if node1 in self.node_to_children[self.node_to_mother[node2]] and self.node_to_father[node1] != self.node_to_father[node2]:
                        return True

            if relation == "maternal aunt/uncle-nephew/niece":
                for sibling in self.node_to_siblings[node1]:
                    if self.node_to_data[sibling]["sex"] == "F" and node2 in self.node_to_children[sibling]:
                        return True
            if relation == "paternal aunt/uncle-nephew/niece":
                for sibling in self.node_to_siblings[node1]:
                    if self.node_to_data[sibling]["sex"] == "M" and node2 in self.node_to_children[sibling]:
                        return True
            if relation == "maternal nephew/niece-aunt/uncle":
                for sibling in self.node_to_siblings[node2]:
                    if self.node_to_data[sibling]["sex"] == "F" and node1 in self.node_to_children[sibling]:
                        return True
            if relation == "paternal nephew/niece-aunt/uncle":
                for sibling in self.node_to_siblings[node2]:
                    if self.node_to_data[sibling]["sex"] == "M" and node1 in self.node_to_children[sibling]:
                        return True

            if relation == "maternal grandparent-grandchild":
                for child in self.node_to_children[node1]:
                    if self.node_to_data[child]["sex"] == "F" and node2 in self.node_to_children[child]:
                        return True
            if relation == "paternal grandparent-grandchild":
                for child in self.node_to_children[node1]:
                    if self.node_to_data[child]["sex"] == "M" and node2 in self.node_to_children[child]:
                        return True
            if relation == "maternal grandchild-grandparent":
                for child in self.node_to_children[node2]:
                    if self.node_to_data[child]["sex"] == "F" and node1 in self.node_to_children[child]:
                        return True
            if relation == "paternal grandchild-grandparent":
                for child in self.node_to_children[node2]:
                    if self.node_to_data[child]["sex"] == "M" and node1 in self.node_to_children[child]:
                        return True

            if relation == "paternal half-siblings":
                if self.node_to_father[node2]:
                    if node1 in self.node_to_children[self.node_to_father[node2]] and self.node_to_mother[node1] != self.node_to_mother[node2]:
                        return True
            if relation == "maternal half-siblings":
                if self.node_to_mother[node2]:
                    if node1 in self.node_to_children[self.node_to_mother[node2]] and self.node_to_father[node1] != self.node_to_father[node2]:
                        return True
        return False

    def get_relations_between_nodes(self, node1: str, node2: str, include_maternal_paternal: bool = False) -> defaultdict[str, int]:
        """
        Returns a dictionary of the *1st- and 2nd-degree* relations between node1 and node2.
        """
        relations = defaultdict(int)
        if self.is_relation_in_pedigree(node1, node2, ["parent-child"]):
            relations["parent-child"] += 1
        if self.is_relation_in_pedigree(node1, node2, ["child-parent"]):
            relations["child-parent"] += 1
        if self.is_relation_in_pedigree(node1, node2, ["siblings"]):
            relations["siblings"] += 1
        
        if self.is_relation_in_pedigree(node1, node2, ["maternal aunt/uncle-nephew/niece"]):
            relations["maternal aunt/uncle-nephew/niece"] += 1
        if self.is_relation_in_pedigree(node1, node2, ["paternal aunt/uncle-nephew/niece"]):
            relations["paternal aunt/uncle-nephew/niece"] += 1
        if self.is_relation_in_pedigree(node1, node2, ["maternal nephew/niece-aunt/uncle"]):
            relations["maternal nephew/niece-aunt/uncle"] += 1
        if self.is_relation_in_pedigree(node1, node2, ["paternal nephew/niece-aunt/uncle"]):
            relations["paternal nephew/niece-aunt/uncle"] += 1
        
        if self.is_relation_in_pedigree(node1, node2, ["maternal grandparent-grandchild"]):
            relations["maternal grandparent-grandchild"] += 1
        if self.is_relation_in_pedigree(node1, node2, ["paternal grandparent-grandchild"]):
            relations["paternal grandparent-grandchild"] += 1
        if self.is_relation_in_pedigree(node1, node2, ["maternal grandchild-grandparent"]):
            relations["maternal grandchild-grandparent"] += 1
        if self.is_relation_in_pedigree(node1, node2, ["paternal grandchild-grandparent"]):
            relations["paternal grandchild-grandparent"] += 1
        
        if self.is_relation_in_pedigree(node1, node2, ["maternal half-siblings"]):
            relations["maternal half-siblings"] += 1
        if self.is_relation_in_pedigree(node1, node2, ["paternal half-siblings"]):
            relations["paternal half-siblings"] += 1
        
        if not include_maternal_paternal:
            relations["aunt/uncle-nephew/niece"] = relations["maternal aunt/uncle-nephew/niece"] + relations["paternal aunt/uncle-nephew/niece"]
            relations["nephew/niece-aunt/uncle"] = relations["maternal nephew/niece-aunt/uncle"] + relations["paternal nephew/niece-aunt/uncle"]
            relations["grandparent-grandchild"] = relations["maternal grandparent-grandchild"] + relations["paternal grandparent-grandchild"]
            relations["grandchild-grandparent"] = relations["maternal grandchild-grandparent"] + relations["paternal grandchild-grandparent"]
            relations["half-siblings"] = relations["maternal half-siblings"] + relations["paternal half-siblings"]
            for relation in ["maternal aunt/uncle-nephew/niece", "paternal aunt/uncle-nephew/niece", 
                             "maternal nephew/niece-aunt/uncle", "paternal nephew/niece-aunt/uncle", 
                             "maternal grandparent-grandchild", "paternal grandparent-grandchild", 
                             "maternal grandchild-grandparent", "paternal grandchild-grandparent", 
                             "maternal half-siblings", "paternal half-siblings"]:
                del relations[relation]
        
        relations_to_remove = set()
        for relation, count in relations.items():
            if count == 0:
                relations_to_remove.add(relation)
        for relation in relations_to_remove:
            del relations[relation]
        return relations

    def get_parent_child_pairs(self, include_placeholders: bool = True) -> list[tuple[str, str]]:
        """
        Gets all (parent, child) pairs in the tree.
        """
        parent_child_pairs: list[tuple[str, str]] = []
        for parent in self.node_to_children:
            for child in self.node_to_children[parent]:
                if include_placeholders or (not parent.isnumeric() and not child.isnumeric()):
                    parent_child_pairs.append((parent, child))
        return parent_child_pairs

    def get_sibling_pairs(self, include_placeholders: bool = True) -> list[tuple[str, str]]:
        """
        Gets all (sibling, sibling) pairs in the tree.
        Note: Only gets *full* siblings. See self.get_half_sibling_pairs().
        """
        sibling_pairs: list[tuple[str, str]] = []
        for sibling1 in self.node_to_siblings:
            for sibling2 in self.node_to_siblings[sibling1]:
                if include_placeholders or (not sibling1.isnumeric() and not sibling2.isnumeric()):
                    if (sibling2, sibling1) not in sibling_pairs:  # Don't add symmetric duplicates
                        sibling_pairs.append((sibling1, sibling2))
        return sibling_pairs

    def get_aunt_uncle_nephew_niece_pairs(self, include_placeholders: bool = True, shared_relative_sex: str | None = None) -> list[tuple[str, str]]:
        """
        Gets all (aunt/uncle, nephew/niece) pairs in the tree.
        Includes duplicates if, for example, shared_relative_sex=None and an aunt is both a maternal and paternal aunt to a nephew (i.e., full-sib mating).
        """
        aunt_uncle_nephew_niece_pairs: list[tuple[str, str]] = []
        for parent, child in self.get_parent_child_pairs():
            for parent_sibling in self.node_to_siblings[parent]:
                if not shared_relative_sex or self.node_to_data[parent]["sex"] == shared_relative_sex:
                    if include_placeholders or (not parent_sibling.isnumeric() and not child.isnumeric()):
                        aunt_uncle_nephew_niece_pairs.append((parent_sibling, child))
        return aunt_uncle_nephew_niece_pairs

    def get_grandparent_grandchild_pairs(self, include_placeholders: bool = True, shared_relative_sex: str | None = None) -> list[tuple[str, str]]:
        """
        Gets all (grandparent, grandchild) pairs in the tree.
        Includes duplicates if, for example, a grandparent is both a maternal and paternal grandparent to a grandchild.
        """
        grandparent_grandchild_pairs: list[tuple[str, str]] = []
        for parent, child in self.get_parent_child_pairs():
            for child_child in self.node_to_children[child]:
                if not shared_relative_sex or self.node_to_data[child]["sex"] == shared_relative_sex:
                    if include_placeholders or (not parent.isnumeric() and not child_child.isnumeric()):
                        grandparent_grandchild_pairs.append((parent, child_child))
        return grandparent_grandchild_pairs

    def get_half_sibling_pairs(self, include_placeholders: bool = True, shared_relative_sex: str | None = None) -> list[tuple[str, str]]:
        """
        Gets all (half-sibling, half-sibling) pairs in the tree.
        """
        half_sibling_pairs: list[tuple[str, str]] = []
        for parent, child in self.get_parent_child_pairs():
            for other_child in self.node_to_children[parent]:
                if child != other_child and other_child not in self.node_to_siblings[child]:
                    if not shared_relative_sex or self.node_to_data[parent]["sex"] == shared_relative_sex:
                        if include_placeholders or (not child.isnumeric() and not other_child.isnumeric()):
                            if (other_child, child) not in half_sibling_pairs:  # Don't add symmetric duplicates
                                half_sibling_pairs.append((child, other_child))
        return half_sibling_pairs

    def get_half_aunt_uncle_nephew_niece_pairs(self, include_placeholders: bool = True) -> list[tuple[str, str]]:
        """
        Gets all (half-aunt/half-uncle, half-nephew/half-niece) pairs in the tree.
        """
        half_aunt_uncle_nephew_niece_pairs: list[tuple[str, str]] = []
        for half_sibling1, half_sibling2 in self.get_half_sibling_pairs():
            for half_sibling1_child in self.node_to_children[half_sibling1]:
                if half_sibling1_child != half_sibling2:
                    if include_placeholders or (not half_sibling2.isnumeric() and not half_sibling1_child.isnumeric()):
                        half_aunt_uncle_nephew_niece_pairs.append((half_sibling2, half_sibling1_child))

            for half_sibling2_child in self.node_to_children[half_sibling2]:
                if half_sibling2_child != half_sibling1:
                    if include_placeholders or (not half_sibling1.isnumeric() and not half_sibling2_child.isnumeric()):
                        half_aunt_uncle_nephew_niece_pairs.append((half_sibling1, half_sibling2_child))
        return half_aunt_uncle_nephew_niece_pairs

    def get_greatgrandparent_greatgrandchild_pairs(self, include_placeholders: bool = True) -> list[tuple[str, str]]:
        """
        Gets all (greatgrandparent, greatgrandchild) pairs in the tree.
        """
        greatgrandparent_greatgrandchild_pairs: list[tuple[str, str]] = []
        for grandparent, grandchild in self.get_grandparent_grandchild_pairs():
            for grandchild_child in self.node_to_children[grandchild]:
                if include_placeholders or (not grandparent.isnumeric() and not grandchild_child.isnumeric()):
                    greatgrandparent_greatgrandchild_pairs.append((grandparent, grandchild_child))
        return greatgrandparent_greatgrandchild_pairs

    def get_grandaunt_granduncle_grandnephew_grandniece_pairs(self, include_placeholders: bool = True) -> list[tuple[str, str]]:
        """
        Gets all (grandaunt/uncle, grandnephew/niece) pairs in the tree.
        """
        grandaunt_granduncle_grandnephew_grandniece_pairs: list[tuple[str, str]] = []
        for aunt_uncle, nephew_niece in self.get_aunt_uncle_nephew_niece_pairs():
            for nephew_niece_child in self.node_to_children[nephew_niece]:
                if include_placeholders or (not aunt_uncle.isnumeric() and not nephew_niece_child.isnumeric()):
                    grandaunt_granduncle_grandnephew_grandniece_pairs.append((aunt_uncle, nephew_niece_child))
        return grandaunt_granduncle_grandnephew_grandniece_pairs

    def get_first_cousin_pairs(self, include_placeholders: bool = True) -> list[tuple[str, str]]:
        """
        Gets all (first cousin, first cousin) pairs in the tree.
        """
        cousin_pairs: list[tuple[str, str]] = []
        for aunt_uncle, child in self.get_aunt_uncle_nephew_niece_pairs():
            for aunt_uncle_child in self.node_to_children[aunt_uncle]:
                if include_placeholders or (not child.isnumeric() and not aunt_uncle_child.isnumeric()):
                    if aunt_uncle_child != child and (aunt_uncle_child, child) not in cousin_pairs:  # Don't add symmetric duplicates
                        cousin_pairs.append((child, aunt_uncle_child))
        return cousin_pairs

    def get_related_pairs(self, include_placeholders: bool = True) -> set[tuple[str, str]]:
        """
        Gets all related pairs (up to and including 3rd-degree relations) in the pedigree.
        """
        related_pairs: set[tuple[str, str]] = set()
        related_pairs.update(self.get_parent_child_pairs(include_placeholders=include_placeholders))
        related_pairs.update(self.get_sibling_pairs(include_placeholders=include_placeholders))
        related_pairs.update(self.get_aunt_uncle_nephew_niece_pairs(include_placeholders=include_placeholders))
        related_pairs.update(self.get_grandparent_grandchild_pairs(include_placeholders=include_placeholders))
        related_pairs.update(self.get_half_sibling_pairs(include_placeholders=include_placeholders))
        related_pairs.update(self.get_half_aunt_uncle_nephew_niece_pairs(include_placeholders=include_placeholders))
        related_pairs.update(self.get_greatgrandparent_greatgrandchild_pairs(include_placeholders=include_placeholders))
        related_pairs.update(self.get_grandaunt_granduncle_grandnephew_grandniece_pairs(include_placeholders=include_placeholders))
        related_pairs.update(self.get_first_cousin_pairs(include_placeholders=include_placeholders))
        return related_pairs

    def get_non_placeholder_nodes(self) -> set[str]:
        """
        Gets all non-placeholder nodes in the tree.
        """
        return set([node for node in self.node_to_data if not node.isnumeric()])

    def clean_up_relations(self) -> None:
        """
        Remove any empty entries in the relation dictionaries.
        Also remove unnecessary placeholder nodes to standardize topological sort output.
        """
        for relation_dict in [self.node_to_father, self.node_to_mother, self.node_to_children, self.node_to_siblings]:
            keys_to_remove = set()
            for k, v in relation_dict.items():
                if not v:
                    keys_to_remove.add(k)
            for key_to_remove in keys_to_remove:
                del relation_dict[key_to_remove]

        placeholder_nodes_to_remove: set[str] = set()
        for node in self.node_to_data:
            mother = self.node_to_mother[node]
            father = self.node_to_father[node]
            if mother.isnumeric() and father.isnumeric():
                if len(self.node_to_children[mother]) == 1 and len(self.node_to_children[father]) == 1:
                    if (not self.node_to_mother[mother] and not self.node_to_father[mother]
                        and not self.node_to_mother[father] and not self.node_to_father[father]):
                        placeholder_nodes_to_remove.add(mother)
                        placeholder_nodes_to_remove.add(father)
        
        for node in placeholder_nodes_to_remove:
            for data_dict in [self.node_to_data, self.node_to_father, self.node_to_mother, self.node_to_children, self.node_to_siblings]:
                if node in data_dict:
                    del data_dict[node]
            
            for node in self.node_to_data:
                assert node not in self.node_to_siblings[node] and node not in self.node_to_children[node]
                if self.node_to_father[node] in placeholder_nodes_to_remove:
                    del self.node_to_father[node]
                if self.node_to_mother[node] in placeholder_nodes_to_remove:
                    del self.node_to_mother[node]

    def plot(self, path: str) -> None:
        """
        Plot the pedigree to the given path.
        """
        tree = nx.from_dict_of_lists(self.node_to_children, create_using=nx.DiGraph)
        for node in self.node_to_data:  # Add childless nodes
            if node not in tree.nodes:
                tree.add_node(node)
        male_named_nodes = [node for node in self.node_to_data if self.node_to_data[node]["sex"] == "M" and not node.isnumeric()]
        male_placeholder_nodes = [node for node in self.node_to_data if self.node_to_data[node]["sex"] == "M" and node.isnumeric()]
        female_named_nodes = [node for node in self.node_to_data if self.node_to_data[node]["sex"] == "F" and not node.isnumeric()]
        female_placeholder_nodes = [node for node in self.node_to_data if self.node_to_data[node]["sex"] == "F" and node.isnumeric()]
        
        non_placeholder_labels = dict()
        for node in tree.nodes:
            mt_haplogroup = self.node_to_data[node]["mt_haplogroup"][:3]
            y_haplogroup = self.node_to_data[node]["y_haplogroup"][:3]
            if node.isnumeric():
                if y_haplogroup:
                    non_placeholder_labels[node] = f"{mt_haplogroup}\n{y_haplogroup}"
                else:
                    non_placeholder_labels[node] = f"{mt_haplogroup}"
            else:
                if y_haplogroup:
                    non_placeholder_labels[node] = f"{node}\n{mt_haplogroup}\n{y_haplogroup}"
                else:
                    non_placeholder_labels[node] = f"{node}\n{mt_haplogroup}"
        
        # Create colormap for MT haplogroups
        cmap = plt.get_cmap("tab20")
        mt_haplogroups = sorted(set([self.node_to_data[node]["mt_haplogroup"] for node in self.node_to_data if not node.isnumeric()]))
        mt_haplogroup_to_color = {haplogroup: cmap(i / len(mt_haplogroups)) for i, haplogroup in enumerate(mt_haplogroups)}
        male_named_node_colors = [mt_haplogroup_to_color[self.node_to_data[node]["mt_haplogroup"]] for node in male_named_nodes]
        female_named_node_colors = [mt_haplogroup_to_color[self.node_to_data[node]["mt_haplogroup"]] for node in female_named_nodes]

        plt.figure(figsize=(12, 4.8), dpi=1200)
        node_size = min(1000, 10000 / len(tree.nodes))  # Scale sizes based on pedigree node count
        font_size = min(4.5, 150 / len(tree.nodes))

        pos = nx.nx_agraph.graphviz_layout(tree, prog="dot")
        nx.draw_networkx_nodes(tree, pos=pos, nodelist=male_named_nodes, node_shape="s", node_size=node_size, node_color=male_named_node_colors, edgecolors="black", linewidths=0.2)
        nx.draw_networkx_nodes(tree, pos=pos, nodelist=female_named_nodes, node_shape="o", node_size=node_size, node_color=female_named_node_colors, edgecolors="black", linewidths=0.2)
        nx.draw_networkx_nodes(tree, pos=pos, nodelist=male_placeholder_nodes, node_shape="s", node_size=node_size, node_color="#e5e5e5", edgecolors="black", linewidths=0.2)
        nx.draw_networkx_nodes(tree, pos=pos, nodelist=female_placeholder_nodes, node_shape="o", node_size=node_size, node_color="#e5e5e5", edgecolors="black", linewidths=0.2)
        nx.draw_networkx_labels(tree, pos=pos, labels=non_placeholder_labels, font_size=font_size)
        nx.draw_networkx_edges(tree, pos=pos, node_shape="s", node_size=node_size, width=0.2, arrowsize=font_size*1.2)

        plt.axis("off")
        plt.savefig(path, bbox_inches="tight")
        plt.close()
