import pandas as pd
from scipy.spatial.distance import hamming, pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

# Load the reference and input data
ref_df = pd.read_excel('/Users/savaglisic/Desktop/reference.xlsx')
input_df = pd.read_excel('/Users/savaglisic/Desktop/input.xlsx')

# Remove the non-genotype columns (SSR and 'Allele')
ref_genotypes = ref_df.drop(columns=['SSR', 'Allele'])
input_genotype = input_df.drop(columns=['SSR', 'Allele']).iloc[:, 0]  # Convert the input dataframe to a Series

def calculate_distances(ref_genotypes, input_genotype):
    # Initialize a dictionary to store the Hamming distances
    distances = {}

    # Calculate the Hamming distance for each genotype in the reference data
    for col in ref_genotypes.columns:
        distances[col] = hamming(ref_genotypes[col], input_genotype)

    return distances

# Calculate the distances
distances = calculate_distances(ref_genotypes, input_genotype)

# If an exact match is found
if 0 in distances.values():
    match = [genotype for genotype, dist in distances.items() if dist == 0]
    print(f'Exact match found: {match[0]}')

# Sort distances in ascending order
sorted_distances = sorted(distances.items(), key=lambda x: x[1])

# Select top N most similar genotypes
N = 7  # adjust as needed
top_N_genotypes = sorted_distances[:N]

# Print the most similar genotypes with percentages
for genotype, dist in top_N_genotypes:
    print(f'{genotype}: {100 * (1 - dist):.2f}%')

# Convert the distances to a format suitable for scipy's clustering functions
dist_matrix = pd.DataFrame.from_dict(dict(top_N_genotypes), orient='index')
dist_array = pdist(dist_matrix)

# Perform the hierarchical clustering
clusters = linkage(dist_array, method='average')

# Create a dendrogram
plt.figure(figsize=(10, 7))
dendrogram(clusters, labels=dist_matrix.index.tolist())
plt.show()
