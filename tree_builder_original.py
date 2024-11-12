import pandas as pd
import numpy as np
from skbio import TreeNode, DistanceMatrix
from skbio.tree import nj

# Step 1: Read the Similarity Matrix from CSV
def read_distance_csv(file_path):
    # Reads a CSV file containing the similarity matrix.
    similarity_df = pd.read_csv(file_path, index_col=0)
    return similarity_df

# Step 3: Construct the Phylogenetic Tree using Neighbor-Joining
def construct_tree(distance_df):
    # Create a list of sample names
    ids = distance_df.index.tolist()
    # Convert the distance DataFrame to a numpy array
    distance_array = distance_df.values
    # Create a scikit-bio DistanceMatrix
    dm = DistanceMatrix(distance_array, ids)
    # Construct the tree using Neighbor-Joining
    tree = nj(dm)
    return tree

# Step 4: Main Function to Tie Everything Together
def main(distance_path, output_tree_path):
    # Read the similarity matrix
    distance_df = read_distance_csv(distance_path)
    print("Similarity matrix loaded.")

    # Construct the phylogenetic tree
    tree = construct_tree(distance_df)
    print("Phylogenetic tree constructed using Neighbor-Joining method.")

    # Optional: Write the tree to a file in Newick format
    tree.write(output_tree_path)
    print(f"Phylogenetic tree saved to {output_tree_path}.")

    # Optional: Print the tree to the console
    print("\nConstructed Phylogenetic Tree:")
    print(tree.ascii_art())

# Example usage
if __name__ == "__main__":
    distance_path = "/users/harry/desktop/Computational Genetics/Final project/distance_matrix.csv"  # Replace with your CSV file path
    output_tree_path = "/users/harry/desktop/Computational Genetics/Final project/phylogenetic_tree.nwk"     # Output Newick file path
    main(distance_path, output_tree_path)