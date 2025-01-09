import os
from relation_comparison import RelationComparison

def main():
    data_dir = os.path.join(os.path.dirname(__file__), "rivollat")
    published_relations_path = os.path.join(data_dir, "relations_published.csv")
    algorithm_nodes_path = os.path.join(data_dir, "nodes.csv")
    algorithm_relations_path = os.path.join(data_dir, "relations_READv2.csv")

    relation_comparison = RelationComparison(published_relations_path=published_relations_path,
                                             algorithm_nodes_path=algorithm_nodes_path,
                                             algorithm_relations_path=algorithm_relations_path)
    for relation, count in relation_comparison.get_metrics().items():
        print(f"{relation}: {round(count, 2)}")

if __name__ == "__main__":
    main()
