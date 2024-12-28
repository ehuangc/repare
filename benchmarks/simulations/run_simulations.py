from simulated_pedigree import SimulatedPedigree


def simulate():
    simulated_pedigree = SimulatedPedigree()
    simulated_pedigree.create_pedigree()
    simulated_pedigree.mask_data()
    simulated_pedigree.run_algorithm()
    metrics = simulated_pedigree.get_metrics()
    print(metrics)

def run_experiment():
    simulate()

def main():
    run_experiment()

if __name__ == "__main__":
    main()
