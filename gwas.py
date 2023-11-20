import scipy.stats as stats


def parse_vcf(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Skip the header
    lines = lines[1:]

    # Initialize a dictionary to hold our parsed data
    snp_data = {}

    for line in lines:
        # Split the line by whitespace
        split_line = line.strip().split('\t')

        # The first two columns are not SNPs
        snp_id = split_line[2]  # The SNP ID is in the third column
        genotypes = split_line[5:]  # The genotypes start from the 6th column

        # Convert genotype strings to integers for analysis
        # Assuming '|' or '/' separated alleles and taking only the first allele for simplicity
        # This will need to be adjusted if the format is different
        int_genotypes = []
        for gt in genotypes:
            if gt == "0|0":
                int_genotypes.append(0)
            elif gt == "1|0" or gt == "0|1":
                int_genotypes.append(1)
            else:
                int_genotypes.append(2)

        # Store the genotypes in the dictionary
        snp_data[snp_id] = int_genotypes

    return snp_data


# Function to read phenotype data
def read_phenotypes(filename):
    with open(filename) as f:
        phenotypes = [int(line.split('\t')[1]) for line in f]
        print(f"Count 0: {phenotypes.count(0)}")
        print(f"Count 1: {phenotypes.count(1)}")
    return phenotypes


# Function to read genotype data (this will depend on the actual format of your VCF file)
def read_genotypes(vcf_filename):
    # This is a placeholder function. You'll need to implement the actual parsing.
    genotypes = parse_vcf(vcf_filename)  # A dictionary where the key is the SNP ID and the value is a list of genotypes
    return genotypes


# Perform the GWAS
def gwas(vcf_filename, phenotype_filename):
    phenotypes = read_phenotypes(phenotype_filename)
    genotypes = read_genotypes(vcf_filename)

    # Dictionary to store p-values for each SNP
    p_values = {}

    for snp_id, snp_genotypes in genotypes.items():
        # Construct the contingency table for the current SNP
        contingency_table = [[0, 0, 0],
                             [0, 0, 0]]  # [[disease_yes_ref, disease_no_ref], [disease_yes_alt, disease_no_alt]]

        for genotype, disease in zip(snp_genotypes, phenotypes):
            if genotype == 0:  # Assuming '0' is homozygous reference allel
                if disease == 1:
                    contingency_table[0][0] += 1
                else:
                    contingency_table[1][0] += 1
            elif genotype == 1:  # Assuming '1' is heterozygous allel
                if disease == 1:
                    contingency_table[0][1] += 1
                else:
                    contingency_table[1][1] += 1
            else:  # Assuming '2' is homozygous alternate allel
                if disease == 1:
                    contingency_table[0][2] += 1
                else:
                    contingency_table[1][2] += 1

        while 0 in contingency_table[0]:
            index = contingency_table[0].index(0)
            contingency_table[0].remove(contingency_table[0][index])
            contingency_table[1].remove(contingency_table[1][index])

        # Perform the chi-squared test
        chi2, p, dof, expected = stats.chi2_contingency(contingency_table, correction=False)
        p_values[snp_id] = p

    return p_values


# Perform the GWAS and Bonferroni correction
def gwas_with_bonferroni_correction(vcf_filename, phenotype_filename, alpha=0.05):
    phenotypes = read_phenotypes(phenotype_filename)
    genotypes = read_genotypes(vcf_filename)

    # Total number of tests
    num_tests = len(genotypes)

    # Bonferroni corrected alpha
    corrected_alpha = alpha / 1000

    # Dictionary to store results
    results = {}

    for snp_id, snp_genotypes in genotypes.items():
        # Construct the contingency table for the current SNP
        contingency_table = [[1, 1, 1],
                             [1, 1, 1]]  # [[disease_yes_ref, disease_no_ref], [disease_yes_alt, disease_no_alt]]

        for genotype, disease in zip(snp_genotypes, phenotypes):
            if genotype == 0:  # Assuming '0' is homozygous reference allel
                if disease == 1:
                    contingency_table[0][0] += 1
                else:
                    contingency_table[1][0] += 1
            elif genotype == 1:  # Assuming '1' is heterozygous allel
                if disease == 1:
                    contingency_table[0][1] += 1
                else:
                    contingency_table[1][1] += 1
            else:  # Assuming '2' is homozygous alternate allel
                if disease == 1:
                    contingency_table[0][2] += 1
                else:
                    contingency_table[1][2] += 1

        # while 0 in contingency_table[0]:
        #     index = contingency_table[0].index(0)
        #     contingency_table[0].remove(contingency_table[0][index])
        #     contingency_table[1].remove(contingency_table[1][index])

        # Perform the chi-squared test
        chi2, p, dof, expected = stats.chi2_contingency(contingency_table, correction=False)

        # Store results if p-value is less than the Bonferroni corrected alpha
        if p < corrected_alpha:
            # Calculate odds ratios
            odds_ratio_het = (contingency_table[0][1] / contingency_table[1][1]) / (
                    contingency_table[0][0] / contingency_table[1][0])
            odds_ratio_homo_alt = (contingency_table[0][2] / contingency_table[1][2]) / (
                    contingency_table[0][0] / contingency_table[1][0])
            results[snp_id] = {
                'uncorrected_p': p,
                'corrected_p': p * num_tests,  # Applying Bonferroni correction
                'odds_ratio_het': odds_ratio_het,
                'odds_ratio_homo_alt': odds_ratio_homo_alt
            }

    return results


if __name__ == '__main__':
    # # Use the function and print the results
    # p_values = gwas('gwas_population.vcf', 'gwas_phenotypes.txt')
    # significant_snps = [snp_id for snp_id, p_value in p_values.items() if p_value < 0.05]
    # print(len(significant_snps))
    # Use the function and print the results
    results = gwas_with_bonferroni_correction('gwas_population.vcf', 'gwas_phenotypes.txt')

    # Now print the table as specified
    print("SNP ID\tUncorrected p-value\tCorrected p-value\tOdds Ratio (Het)\tOdds Ratio (Homo-Alt)")
    for snp_id, data in results.items():
        print(
            f"{snp_id}\t{data['uncorrected_p']:.4g}\t{data['corrected_p']:.4g}\t{data['odds_ratio_het']:.4g}\t{data['odds_ratio_homo_alt']:.4g}")
