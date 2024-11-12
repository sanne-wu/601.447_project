import pandas as pd
import numpy as np
from skbio import TreeNode, DistanceMatrix
from skbio.tree import nj

def read_similarity_csv(file_path):
    # Reads a CSV file containing the similarity matrix.
    similarity_df = pd.read_csv(file_path, index_col=0)
    return similarity_df

def similarity_to_distance(similarity_df):
    # Ensure all similarities are between 0 and 1
    similarity_df = similarity_df.clip(lower=0, upper=1)
    epsilon = 1e-10
    similarity_df += epsilon
    distance_df = -np.log(similarity_df)
    np.fill_diagonal(distance_df.values, 0)
    return distance_df

def construct_tree(distance_df):
    ids = distance_df.index.tolist()
    distance_array = distance_df.values
    dm = DistanceMatrix(distance_array, ids)
    tree = nj(dm)
    return tree

def main(similarity_csv_path, output_tree_path):
    similarity_df = read_similarity_csv(similarity_csv_path)
    print("Similarity matrix loaded.")
    distance_df = similarity_to_distance(similarity_df)
    print("Similarity matrix converted to distance matrix.")
    tree = construct_tree(distance_df)
    print("Phylogenetic tree constructed using Neighbor-Joining method.")
    tree.write(output_tree_path)
    print(f"Phylogenetic tree saved to {output_tree_path}.")
    print("\nConstructed Phylogenetic Tree:")
    print(tree.ascii_art())

if __name__ == "__main__":
    similarity_csv_path = "/users/harry/desktop/Computational Genetics/Final project/code/cds_similarity_matrix.csv"  # Replace with your CSV file path
    output_tree_path = "/users/harry/desktop/Computational Genetics/Final project/phylogenetic_tree.nwk"     # Output Newick file path
    main(similarity_csv_path, output_tree_path)