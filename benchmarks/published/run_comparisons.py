import os
from relation_comparison import RelationComparison

def main():
    data_dir = os.path.join(os.path.dirname(__file__), "fowler")
    published_relations_path = os.path.join(data_dir, "published_relations.csv")
    algorithm_nodes_path = os.path.join(data_dir, "algorithm_nodes.csv")
    algorithm_relations_path = os.path.join(data_dir, "algorithm_relations_custom.csv")

    relation_comparison = RelationComparison(published_relations_path=published_relations_path,
                                             algorithm_nodes_path=algorithm_nodes_path,
                                             algorithm_relations_path=algorithm_relations_path)
    for relation, count in relation_comparison.get_metrics().items():
        print(f"{relation}: {round(count, 2)}")

if __name__ == "__main__":
    main()
