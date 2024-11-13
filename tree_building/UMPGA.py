from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def save_to_phylip(distances, labels, file_path):
    """
    Saves the distance matrix to a Phylip-formatted file.
    """
    n = len(labels)
    with open(file_path, 'w') as f:
        f.write(f"{n}\n")
        for i, label in enumerate(labels):
            # Phylip requires labels to be truncated or padded to 10 characters
            formatted_label = label[:10].ljust(10)
            f.write(f"{formatted_label}")
            f.write(" ".join(f"{distances[i, j]:.5f}" for j in range(n)))
            f.write("\n")

def load_distance_matrix(file_path):
    # Load the CSV file with "-" as missing values in the lower triangle and diagonal
    df = pd.read_csv(file_path, index_col=0)
    # Replace "-" with 0 for handling
    df.replace("-", 0, inplace=True)
    # Convert DataFrame to a NumPy array, with NaN in the lower triangle and diagonal
    distance_matrix = df.to_numpy(dtype=float)
    # Fill the lower triangle by mirroring the upper triangular part
    i_upper = np.triu_indices_from(distance_matrix, 1)
    distance_matrix[(i_upper[1], i_upper[0])] = distance_matrix[i_upper]
    return distance_matrix, df.columns.tolist()

def upgma_tree(distance_matrix, labels, title, ax, save=False):
    # Perform UPGMA clustering using the average linkage method
    linkage_matrix = linkage(distance_matrix[np.triu_indices_from(distance_matrix, 1)], method='weighted')
    # Plot the dendrogram
    dendrogram(linkage_matrix, labels=labels, orientation='left', leaf_rotation=0, leaf_font_size=10, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Distance")
    ax.set_ylabel("Clusters")
    if save:
        # Save as a Phylip distance matrix
        save_to_phylip(distance_matrix, labels, "UMPGA_Dashing.phy")


# Load the three distance matrices and labels
file_paths = ['distance_matrix_whole_genome.csv', 'distance_matrix_gene_level.csv', 'distance_matrix_cds.csv']
distance_matrices = []
labels = None

for file_path in file_paths:
    matrix, labels = load_distance_matrix(file_path)
    distance_matrices.append(matrix)

# Calculate the consensus distance matrix by averaging the three matrices
consensus_matrix = np.nanmean(distance_matrices, axis=0)
print(consensus_matrix)
# Create a 2x2 subplot layout
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Generate individual UPGMA trees for each matrix
upgma_tree(distance_matrices[0], labels, "Whole Genome UPGMA Tree", axes[0, 0])
upgma_tree(distance_matrices[1], labels, "Gene Level UPGMA Tree", axes[0, 1])
upgma_tree(distance_matrices[2], labels, "CDS Level UPGMA Tree", axes[1, 0])

# Generate the consensus UPGMA tree
upgma_tree(consensus_matrix, labels, "Consensus UPGMA Tree", axes[1, 1], save=True)

# Adjust layout and display the plot
plt.tight_layout()
plt.show()
