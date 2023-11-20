import random

import matplotlib
import matplotlib.pyplot as plt


def load_vcf(filename):
    individuals = {}

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('#'):
                # This is the header line
                parts = line.strip().split('\t')
                individual_ids = parts[6:]
                for individual_id in individual_ids:
                    individuals[individual_id] = {"Maternal chromosome": "", "Paternal chromosome": ""}
            else:
                # These are the data lines with genotype information
                parts = line.strip().split('\t')
                genotype_info = parts[6:]  # Extract genotype information
                for individual_id, genotype in zip(individual_ids, genotype_info):
                    if genotype == "0|0":
                        individuals[individual_id]["Maternal chromosome"] += "0"
                        individuals[individual_id]["Paternal chromosome"] += "0"
                    elif genotype == "1|0":
                        individuals[individual_id]["Maternal chromosome"] += "1"
                        individuals[individual_id]["Paternal chromosome"] += "0"
                    elif genotype == "0|1":
                        individuals[individual_id]["Maternal chromosome"] += "0"
                        individuals[individual_id]["Paternal chromosome"] += "1"
                    else:
                        individuals[individual_id]["Maternal chromosome"] += "1"
                        individuals[individual_id]["Paternal chromosome"] += "1"

    return individuals


def select_parent(population):
    """Select a parent randomly from the population."""
    return random.choice(list(population.values()))


def reproduce(parent1, parent2, L):
    """Produce offspring from two parents."""
    crossover_point = random.randint(0, L - 1)

    child_maternal_chromosome = parent1["Maternal chromosome"][:crossover_point] + parent2["Maternal chromosome"][
                                                                                   crossover_point:]
    child_paternal_chromosome = parent1["Paternal chromosome"][:crossover_point] + parent2["Paternal chromosome"][
                                                                                   crossover_point:]

    return {"Maternal chromosome": child_maternal_chromosome, "Paternal chromosome": child_paternal_chromosome}


def simulate_generation(population, L):
    """Simulate a generation."""
    new_population = {}
    for i in range(len(population)):
        parent1 = select_parent(population)
        parent2 = select_parent(population)
        while parent1 == parent2:
            parent2 = select_parent(population)

        child = reproduce(parent1, parent2, L)
        new_population[f"Individual_{i}"] = child

    return new_population


def count_extinct_snps(population, L):
    """Count the number of SNPs that went extinct."""
    extinct_count = 0
    for i in range(L):
        if all(individual["Maternal chromosome"][i] == '0' and individual["Paternal chromosome"][i] == '0' for
               individual in population.values()):
            extinct_count += 1
    return extinct_count


def update_allele_frequencies(population, L):
    """
    Update the allele frequencies for the population.
    """
    allele_frequencies = [0] * L
    for individual in population.values():
        for i, (maternal, paternal) in enumerate(
                zip(individual["Maternal chromosome"][:L], individual["Paternal chromosome"][:L])):
            allele_frequencies[i] += int(maternal) + int(paternal)
    # Divide by the total number of alleles for each SNP to get the frequency
    allele_frequencies = [freq / (2 * len(population)) for freq in allele_frequencies]
    return allele_frequencies


def plot_allele_frequencies():
    """
    Plot the allele frequencies.
    """
    # Initialize variables
    vcf_file_path = 'initial_population.vcf'
    initial_population = load_vcf(vcf_file_path)
    L = 100  # We are only interested in the first 100 SNPs
    num_generations = 20
    allele_frequency_over_time = [[] for _ in range(L)]

    # Simulate the evolution
    current_population = initial_population
    for generation in range(num_generations):
        # Update population
        current_population = simulate_generation(current_population, L)
        # Update allele frequencies for the first 100 SNPs
        current_frequencies = update_allele_frequencies(current_population, L)
        # Record the frequencies
        for i in range(L):
            allele_frequency_over_time[i].append(current_frequencies[i])

    # Plot the results
    plt.figure(figsize=(20, 12))
    for i in range(L):
        plt.plot(range(num_generations), allele_frequency_over_time[i], label=f'SNP {i + 1}')
    # plt.rc('axes', labelsize=100)  # fontsize of the x and y labels
    plt.xlabel('Generation', fontsize=30)
    plt.ylabel('Alternate allele frequency', fontsize=30)
    plt.title('Alternate Allele Frequency Over 20 Generations', fontsize=30)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))  # Move the legend out of the plot
    plt.show()


def calculate_simulation_prob():
    # Load the initial population data
    vcf_file_path = 'initial_population.vcf'
    initial_population = load_vcf(vcf_file_path)
    print(len(initial_population))
    print(initial_population["IND0"])

    # Define the number of SNPs (L)
    L = len(next(iter(initial_population.values()))["Maternal chromosome"])

    # Simulate evolution for one generation (since we are interested in Gen 0 to Gen 1 transition)
    next_generation_population = simulate_generation(initial_population, L)
    print(len(next_generation_population))
    print(next_generation_population["Individual_0"])

    # Count the number of extinct SNPs in the next generation
    extinct_snp_count = count_extinct_snps(next_generation_population, L)

    # Estimate the probability
    estimated_probability = extinct_snp_count / L
    print(f"Estimated probability of allele extinction in generation 1: {estimated_probability}")

    # Compare with analytical result
    analytical_probability = ((2 * len(initial_population) - 1) / (2 * len(initial_population))) ** (
            2 * len(initial_population))
    print(f"Analytical probability of allele extinction in generation 1: {analytical_probability}")


def simulate_evolution(initial_population, L, num_generations):
    extinction_probabilities = []
    fixation_probabilities = []
    current_population = initial_population

    for generation in range(num_generations):
        current_population = simulate_generation(current_population, L)
        frequencies = update_allele_frequencies(current_population, L)

        # Calculate extinction and fixation probabilities
        extinct = sum(1 for freq in frequencies if freq == 0)
        fixed = sum(1 for freq in frequencies if freq == 1)

        extinction_probabilities.append(extinct / L)
        fixation_probabilities.append(fixed / L)

    return extinction_probabilities, fixation_probabilities


def plot_1e():
    # Load the initial population and define parameters
    vcf_file_path = 'initial_population.vcf'
    initial_population = load_vcf(vcf_file_path)
    L = 10000  # Total number of SNPs to track
    num_generations = 1000  # Number of generations to simulate

    # Simulate evolution and calculate probabilities
    extinction_probabilities, fixation_probabilities = simulate_evolution(initial_population, L, num_generations)

    # Plot the results
    plt.figure(figsize=(14, 7))

    # Plot extinction probabilities
    plt.subplot(1, 2, 1)
    plt.plot(range(num_generations), extinction_probabilities, label='Extinction Probability', color='red')
    plt.xlabel('Generation')
    plt.ylabel('Probability of Extinction')
    plt.title('Probability of Allele Extinction Over Generations')

    # Plot fixation probabilities
    plt.subplot(1, 2, 2)
    plt.plot(range(num_generations), fixation_probabilities, label='Fixation Probability', color='blue')
    plt.xlabel('Generation')
    plt.ylabel('Probability of Fixation')
    plt.title('Probability of Allele Fixation Over Generations')

    plt.tight_layout()
    plt.show()


def calculate_fitness_f(individual, snp_index=42):
    genotype = (individual["Maternal chromosome"][snp_index], individual["Paternal chromosome"][snp_index])
    if genotype in [('0', '1'), ('1', '0')]:
        return 1.5  # Heterozygous at SNP42
    elif genotype == ('1', '1'):
        return 2  # Homozygous alternate at SNP42
    return 1  # Default fitness for other genotypes


def select_parent_f(population, fitness_scores):
    total_fitness = sum(fitness_scores.values())
    selection_probs = [fitness / total_fitness for fitness in fitness_scores.values()]
    selected_individual = random.choices(list(population.keys()), weights=selection_probs, k=1)[0]
    return selected_individual


def reproduce(parent1, parent2, L):
    crossover_point = random.randint(0, L - 1)
    child_maternal_chromosome = parent1["Maternal chromosome"][:crossover_point] + parent2["Maternal chromosome"][
                                                                                   crossover_point:]
    child_paternal_chromosome = parent1["Paternal chromosome"][:crossover_point] + parent2["Paternal chromosome"][
                                                                                   crossover_point:]
    return {"Maternal chromosome": child_maternal_chromosome, "Paternal chromosome": child_paternal_chromosome}


def simulate_generation_f(population, L, snp_index=42):
    new_population = {}
    fitness_scores = {individual_id: calculate_fitness_f(individual, snp_index) for individual_id, individual in
                      population.items()}

    for i in range(len(population)):
        parent1_id = select_parent_f(population, fitness_scores)
        parent2_id = select_parent_f(population, fitness_scores)
        while parent1_id == parent2_id:
            parent2_id = select_parent_f(population, fitness_scores)

        parent1 = population[parent1_id]
        parent2 = population[parent2_id]

        child = reproduce(parent1, parent2, L)
        new_population[f"Individual_{i}"] = child

    return new_population


def estimate_extinction_probability_f(num_runs=1000):
    # Load initial population
    vcf_file_path = 'initial_population.vcf'
    initial_population = load_vcf(vcf_file_path)
    L = 10000  # Total number of SNPs
    snp_index = 42  # SNP of interest

    extinction_count = 0
    for _ in range(num_runs):
        # Simulate one generation
        next_generation_population = simulate_generation_f(initial_population, L, snp_index)

        # Check if the alternate allele at SNP42 is extinct in the next generation
        if all(individual["Maternal chromosome"][snp_index] == '0' and individual["Paternal chromosome"][
            snp_index] == '0' for individual in next_generation_population.values()):
            extinction_count += 1

    # Calculate the probability of extinction
    probability_of_extinction = extinction_count / num_runs
    return probability_of_extinction


def calculate_1f():
    # Call the function to perform the estimation
    probability_of_extinction = estimate_extinction_probability_f()
    print(f"Probability of allele extinction at SNP42 after one generation: {probability_of_extinction}")


def simulate_multiple_generations(population, L, snp_index, num_generations):
    current_population = population
    for _ in range(num_generations):
        current_population = simulate_generation_f(current_population, L, snp_index)
    return current_population


def estimate_probabilities(L, snp_index, num_generations, num_runs):
    extinction_count = 0
    fixation_count = 0

    initial_population = load_vcf('initial_population.vcf')

    for _ in range(num_runs):
        final_population = simulate_multiple_generations(initial_population, L, snp_index, num_generations)

        # Check the status of SNP42 in the final population
        allele_frequencies = update_allele_frequencies(final_population, L)
        if allele_frequencies[snp_index] == 0:
            extinction_count += 1
        elif allele_frequencies[snp_index] == 1:
            fixation_count += 1

    extinction_probability = extinction_count / num_runs
    fixation_probability = fixation_count / num_runs

    return extinction_probability, fixation_probability


def calculate_1g():
    # Parameters
    L = 100  # Total number of SNPs
    snp_index = 42  # SNP of interest
    num_generations = 100  # Number of generations to simulate
    num_runs = 1000  # Number of simulation runs

    # Estimate probabilities
    extinction_probability, fixation_probability = estimate_probabilities(L, snp_index, num_generations, num_runs)
    print(f"Probability of extinction: {extinction_probability}")
    print(f"Probability of fixation: {fixation_probability}")


def calculate_fitness_h(individual, snp_index=42):
    genotype = (individual["Maternal chromosome"][snp_index], individual["Paternal chromosome"][snp_index])
    if genotype in [('0', '1'), ('1', '0')]:
        return 0.9  # Heterozygous at SNP42 (deleterious)
    elif genotype == ('1', '1'):
        return 0.8  # Homozygous alternate at SNP42 (deleterious)
    return 1  # Default fitness for other genotypes


def simulate_generation_h(population, L, snp_index):
    new_population = {}
    fitness_scores = {individual_id: calculate_fitness_h(individual, snp_index) for individual_id, individual in
                      population.items()}

    for i in range(len(population)):
        parent1_id = select_parent_f(population, fitness_scores)
        parent2_id = select_parent_f(population, fitness_scores)
        while parent1_id == parent2_id:
            parent2_id = select_parent_f(population, fitness_scores)

        parent1 = population[parent1_id]
        parent2 = population[parent2_id]

        child = reproduce(parent1, parent2, L)
        new_population[f"Individual_{i}"] = child

    return new_population


def estimate_probabilities_h(L, snp_index, num_generations, num_runs):
    extinction_count = 0
    fixation_count = 0

    initial_population = load_vcf('initial_population.vcf')

    for _ in range(num_runs):
        final_population = simulate_multiple_generations(initial_population, L, snp_index, num_generations)

        # Check the status of SNP42 in the final population
        allele_frequencies = update_allele_frequencies(final_population, L)
        if allele_frequencies[snp_index] == 0:
            extinction_count += 1
        elif allele_frequencies[snp_index] == 1:
            fixation_count += 1

    extinction_probability = extinction_count / num_runs
    fixation_probability = fixation_count / num_runs

    return extinction_probability, fixation_probability


def calculate_1h():
    # Parameters for the simulation
    L = 100  # Reduced number of SNPs for simplicity
    snp_index = 42  # SNP of interest
    num_generations = 100  # Number of generations to simulate
    num_runs = 1000  # Number of simulation runs for averaging

    # Estimate probabilities under deleterious allele conditions
    extinction_probability, fixation_probability = estimate_probabilities_h(L, snp_index, num_generations, num_runs)
    print(f"Probability of extinction (deleterious): {extinction_probability}")
    print(f"Probability of fixation (deleterious): {fixation_probability}")


if __name__ == "__main__":
    # random.seed(123)
    calculate_1h()
    # calculate_1g()
    # calculate_1f()
    # plot_1e()
    # plot_allele_frequencies()
    # calculate_simulation_prob()
