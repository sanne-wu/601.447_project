import pandas as pd
import numpy as np
from skbio import DistanceMatrix
from skbio.tree import nj
from dendropy import Tree as DendroPyTree
from dendropy import TaxonNamespace
import glob
import matplotlib.pyplot as plt
from io import StringIO
from Bio import Phylo
import os
import dendropy
import io

def read_similarity_csv(file_path):
    """
    Reads a CSV file containing the similarity matrix.
    Replaces non-numeric entries with NaN and converts all data to float.
    """
    similarity_df = pd.read_csv(file_path, index_col=0, na_values=['-'])
    similarity_df = similarity_df.astype(float)
    return similarity_df

def similarity_to_distance(similarity_df):
    """
    Converts a similarity matrix to a distance matrix using the transformation:
    distance = -log(similarity + epsilon)
    Ensures that similarity values are clipped between 0 and 1 and adds a small epsilon to avoid log(0).
    """
    similarity_df = similarity_df.clip(lower=0, upper=1)
    epsilon = 1e-10
    similarity_df += epsilon
    distance_df = -np.log(similarity_df)
    np.fill_diagonal(distance_df.values, 0)
    return distance_df

def construct_nj_tree(distance_df):
    """
    Constructs a Neighbor-Joining (NJ) phylogenetic tree from a distance matrix using skbio.
    """
    ids = distance_df.index.tolist()
    distance_array = distance_df.values
    dm = DistanceMatrix(distance_array, ids)
    tree = nj(dm)
    return tree

def skbio_to_dendropy_tree(skbio_tree):
    """
    Converts a skbio NJ tree to a DendroPy Tree object.
    """
    # Use StringIO to capture the Newick output as a string
    newick_io = io.StringIO()
    skbio_tree.write(newick_io, format='newick')
    newick_str = newick_io.getvalue()
    newick_io.close()
    
    # Convert Newick string to a DendroPy Tree
    dendropy_tree = DendroPyTree.get(data=newick_str, schema="newick", taxon_namespace=TaxonNamespace())
    return dendropy_tree

def visualize_tree(dendropy_tree, output_image_path):
    """
    Visualizes the DendroPy tree and saves it as a PNG image.
    """
    newick_str = dendropy_tree.as_string(schema='newick')

    # Convert Newick string to Bio.Phylo Tree
    handle = StringIO(newick_str)
    phylo_tree = Phylo.read(handle, 'newick')

    # Draw the tree using matplotlib
    fig = plt.figure(figsize=(12, 8))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(phylo_tree, do_show=False, axes=axes)

    # Save the figure to a file
    plt.savefig(output_image_path, format='png', dpi=300)
    plt.close(fig)
    print(f"Annotated phylogenetic tree visualization saved to {output_image_path}.")

def save_tree_nexus(dendropy_tree, output_tree_path):
    """
    Saves the DendroPy tree in Nexus format.
    """
    dendropy_tree.write(path=output_tree_path, schema='nexus', suppress_annotations=True)
    print(f"Annotated phylogenetic tree saved to {output_tree_path} in Nexus format.")

def main(similarity_csv_files, output_dir):
    """
    Main function to process similarity matrices, generate a consensus tree,
    and save the results in the specified output directory.
    
    Parameters:
    - similarity_csv_files: List of paths to the three similarity matrix CSV files.
    - output_dir: Path to the output directory where results will be saved.
    """
    # Ensure the output directory exists
    nn_output_dir = os.path.join(output_dir, "NN")
    os.makedirs(nn_output_dir, exist_ok=True)
    
    # Read all similarity matrices first to determine common labels
    similarity_dfs = []
    all_labels = []
    
    for file_path in similarity_csv_files:
        similarity_df = read_similarity_csv(file_path)
        similarity_dfs.append(similarity_df)
        all_labels.extend(similarity_df.index.tolist())
        print(f"Loaded similarity matrix from '{file_path}'.")
    
    # Determine the intersection of all labels
    common_labels = set(similarity_dfs[0].index)
    for df in similarity_dfs[1:]:
        common_labels = common_labels.intersection(set(df.index))
    
    if not common_labels:
        raise ValueError("No common labels found across all similarity matrices.")
    
    common_labels = sorted(common_labels)  # Sort for consistent ordering
    print(f"Common labels across all matrices: {common_labels}")
    
    # Reindex all similarity matrices to include only common labels and align their order
    distance_matrices = []
    
    for idx, df in enumerate(similarity_dfs):
        # Reindex to common labels
        df = df.reindex(index=common_labels, columns=common_labels)
        
        # Handle missing values by filling with a default similarity score, e.g., 0
        # Alternatively, you can choose to exclude such labels or handle differently
        df = df.fillna(0)
        print(f"Reindexed similarity matrix {idx+1} to include only common labels.")
        
        # Convert similarity to distance
        distance_df = similarity_to_distance(df)
        distance_matrices.append(distance_df)
        print(f"Converted similarity matrix {idx+1} to distance matrix.")
    
    # Calculate consensus distance matrix by averaging the distance matrices
    distance_arrays = [df.values for df in distance_matrices]
    consensus_distance_matrix = np.nanmean(distance_arrays, axis=0)
    consensus_distance_df = pd.DataFrame(consensus_distance_matrix, index=common_labels, columns=common_labels)
    print("Calculated consensus distance matrix by averaging the distance matrices.")
    
    # Construct consensus NJ tree
    consensus_nj_tree_skbio = construct_nj_tree(consensus_distance_df)
    print("Constructed consensus phylogenetic tree using Neighbor-Joining method.")
    
    consensus_dendropy_tree = skbio_to_dendropy_tree(consensus_nj_tree_skbio)
    print("Converted consensus tree to DendroPy format.")
    
    # Save consensus tree in Nexus format
    output_tree_path = os.path.join(nn_output_dir, "Consensus_Phylogenetic_Tree.nex")
    save_tree_nexus(consensus_dendropy_tree, output_tree_path)
    
    # Visualize and save the consensus tree
    output_image_path = os.path.join(nn_output_dir, "Consensus_Phylogenetic_Tree.png")
    visualize_tree(consensus_dendropy_tree, output_image_path)
    
    print("Consensus tree generation and saving completed successfully.")

if __name__ == "__main__":
    # Define the paths to your three similarity matrix CSV files
    similarity_csv_files = [
        '/users/harry/desktop/Computational Genetics/Final project/code/distance_matrix_whole_genome.csv',
        '/users/harry/desktop/Computational Genetics/Final project/code/distance_matrix_gene_level.csv',
        '/users/harry/desktop/Computational Genetics/Final project/code/distance_matrix_cds.csv'
    ]
    
    # Define the output directory
    output_directory = '/users/harry/desktop/Computational Genetics/Final project/code/result'
    
    # Run the main function
    main(similarity_csv_files, output_directory)
