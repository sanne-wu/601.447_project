from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from Bio import Phylo
from Bio.Phylo.Newick import Tree, Clade

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
            f.write(f"{formatted_label} ")
            f.write(" ".join(f"{distances[i, j]:.5f}" for j in range(n)))
            f.write("\n")

def load_distance_matrix(file_path):
    """
    Loads a distance matrix from a CSV file, handling missing values.
    """
    # Load the CSV file with "-" as missing values in the lower triangle and diagonal
    df = pd.read_csv(file_path, index_col=0)
    # Replace "-" with NaN to handle missing values
    df.replace("-", np.nan, inplace=True)
    # Convert DataFrame to a NumPy array of floats
    distance_matrix = df.to_numpy(dtype=float)
    # Ensure the distance matrix is symmetric by mirroring the upper triangle
    i_upper = np.triu_indices_from(distance_matrix, 1)
    distance_matrix[(i_upper[1], i_upper[0])] = distance_matrix[i_upper]
    return distance_matrix, df.columns.tolist()

def linkage_to_phylo(linkage_matrix, labels):
    """
    Converts a SciPy linkage matrix to a Biopython Phylo Tree object.
    """
    tree_root, nodelist = to_tree(linkage_matrix, rd=True)
    
    def build_clade(node):
        if node.is_leaf():
            return Clade(name=labels[node.id])
        else:
            clade = Clade()
            clade.clades.append(build_clade(node.left))
            clade.clades.append(build_clade(node.right))
            clade.branch_length = node.dist
            return clade
    
    phylo_tree = Tree(root=build_clade(tree_root))
    return phylo_tree

def upgma_tree(distance_matrix, labels, title, ax, save=False, consensus=False, output_dir=None):
    """
    Performs UPGMA clustering, plots the dendrogram, and optionally saves the tree in Nexus format.
    """
    # Perform UPGMA clustering using the average linkage method
    # Convert the distance matrix to a condensed distance matrix required by linkage
    condensed_dist = distance_matrix[np.triu_indices_from(distance_matrix, 1)]
    linkage_matrix = linkage(condensed_dist, method='average')
    
    # Plot the dendrogram
    dendrogram(linkage_matrix, labels=labels, orientation='left', leaf_rotation=0, leaf_font_size=10, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Distance")
    ax.set_ylabel("Clusters")
    
    if save and consensus and output_dir:
        # Convert linkage matrix to a Phylo tree and save in Nexus format
        tree = linkage_to_phylo(linkage_matrix, labels)
        nexus_path = os.path.join(output_dir, "Consensus_UPGMA_Tree.nex")
        Phylo.write(tree, nexus_path, "nexus")
        print(f"Consensus tree saved in Nexus format at: {nexus_path}")

def main():
    # Define Input and Output Paths
    
    # Specify the paths to your input distance matrix CSV files
    input_file_paths = [
        '/users/harry/desktop/Computational Genetics/Final project/code/distance_matrix_whole_genome.csv',
        '/users/harry/desktop/Computational Genetics/Final project/code/distance_matrix_gene_level.csv',
        '/users/harry/desktop/Computational Genetics/Final project/code/distance_matrix_cds.csv'
    ]
    
    # Specify the output directory where all outputs will be saved
    output_dir = '/users/harry/desktop/Computational Genetics/Final project/code/result/UMPGA'  
    
    # Ensure the Output Directory Exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Load Distance Matrices
    distance_matrices = []
    labels = None
    
    for file_path in input_file_paths:
        if not os.path.isfile(file_path):
            print(f"Error: File '{file_path}' does not exist. Please check the path.")
            return
        matrix, labels = load_distance_matrix(file_path)
        distance_matrices.append(matrix)
        print(f"Loaded distance matrix from '{file_path}'.")
    
    # Calculate Consensus Matrix
    consensus_matrix = np.nanmean(distance_matrices, axis=0)
    print("Consensus Distance Matrix:")
    print(consensus_matrix)
    
    # Plotting Setup
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Generate Individual UPGMA Trees
    tree_titles = [
        "Whole Genome UPGMA Tree",
        "Gene Level UPGMA Tree",
        "CDS Level UPGMA Tree"
    ]
    
    for i in range(3):
        upgma_tree(
            distance_matrix=distance_matrices[i],
            labels=labels,
            title=tree_titles[i],
            ax=axes[i//2, i%2],
            save=False
        )
    
    # Generate Consensus UPGMA Tree
    upgma_tree(
        distance_matrix=consensus_matrix,
        labels=labels,
        title="Consensus UPGMA Tree",
        ax=axes[1, 1],
        save=True,
        consensus=True,
        output_dir=output_dir
    )
    
    # Save Consensus Distance Matrix in Phylip Format
    phylip_path = os.path.join(output_dir, "Consensus_UPGMA_Distance_Matrix.phy")
    save_to_phylip(consensus_matrix, labels, phylip_path)
    print(f"Consensus distance matrix saved in Phylip format at: {phylip_path}")
    
    # Save and Show Plot
    plot_path = os.path.join(output_dir, "UPGMA_Trees.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.show()
    print(f"Dendrogram plots saved at: {plot_path}")
    
    print("UPGMA clustering and tree generation completed successfully.")

if __name__ == "__main__":
    main()
