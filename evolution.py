import pandas as pd
import random


def load_vcf(filename):
    individuals = []
    genotypes = []

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('##'):
                continue  # Skip header lines with meta-information
            if line.startswith('#'):
                # This is the header line with individual names
                parts = line.strip().split('\t')
                individuals = parts[9:]  # Assuming the genotype data starts at the 10th column
            else:
                # These are the data lines with genotype information
                parts = line.strip().split('\t')
                genotype_info = parts[9:]  # Extract genotype information
                genotype_info = [[int(allele) for allele in geno.replace('|', '')] for geno in
                                 genotype_info]  # Convert to integers
                genotypes.append(genotype_info)

    # Convert genotypes to a list of tuples representing each individual's genotype (maternal and paternal)
    population_genotypes = [tuple(zip(*individual)) for individual in zip(*genotypes)]

    return population_genotypes


def calculate_fitness(individual):
    # Ensure individual is a tuple containing two lists: (maternal_chromosome, paternal_chromosome)
    maternal_chromosome, paternal_chromosome = individual
    return sum(maternal_chromosome) + sum(paternal_chromosome)


def select_parent(population, fitness_scores):
    total_fitness = sum(fitness_scores.values())
    selection_probs = [fitness / total_fitness for fitness in fitness_scores.values()]
    selected_index = random.choices(range(len(population)), weights=selection_probs, k=1)[0]
    return selected_index


def reproduce(parent1, parent2, L):
    # Ensure parent1 and parent2 are tuples containing two lists each
    maternal_chromosome1, paternal_chromosome1 = parent1
    maternal_chromosome2, paternal_chromosome2 = parent2

    # Choose random crossover points for maternal and paternal chromosomes
    crossover_point = random.randint(0, L-1)
    # Create the child's chromosomes by crossing over at the chosen points
    child_maternal = maternal_chromosome1[:crossover_point] + maternal_chromosome2[crossover_point:]
    child_paternal = paternal_chromosome1[:crossover_point] + paternal_chromosome2[crossover_point:]

    return child_maternal, child_paternal


def simulate_generation(population, L):
    fitness_scores = {i: calculate_fitness(individual) for i, individual in enumerate(population)}
    new_population = []

    for _ in range(len(population)):
        parent1_index = select_parent(population, fitness_scores)
        parent2_index = select_parent(population, fitness_scores)
        while parent1_index == parent2_index:
            parent2_index = select_parent(population, fitness_scores)

        parent1 = population[parent1_index]
        parent2 = population[parent2_index]

        child = reproduce(parent1, parent2, L)
        new_population.append(child)

    return new_population


if __name__ == "__main__":
    # Load the initial population data
    vcf_file_path = 'initial_population.vcf'
    initial_population = load_vcf(vcf_file_path)

    # Define the number of SNPs (L) and generations to simulate
    number_of_snps = 5  # This should match the actual number of SNPs in the VCF file
    number_of_generations = 10

    # Simulate evolution for a given number of generations
    for _ in range(number_of_generations):
        initial_population = simulate_generation(initial_population, number_of_snps)
        # Analyze population here if necessary

    # Display the final population
    for individual_id, genotype in enumerate(initial_population):
        print(f"Individual {individual_id}: Maternal: {genotype[0]}, Paternal: {genotype[1]}")