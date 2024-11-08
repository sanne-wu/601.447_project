import pandas as pd
import numpy as np
from skbio import TreeNode, DistanceMatrix
from skbio.tree import nj

# Step 1: Read the Similarity Matrix from CSV
def read_similarity_csv(file_path):
    """
    Reads a CSV file containing the similarity matrix.

    :param file_path: Path to the CSV file.
    :return: pandas DataFrame of the similarity matrix.
    """
    similarity_df = pd.read_csv(file_path, index_col=0)
    return similarity_df

# Step 2: Convert Similarity Matrix to Distance Matrix
def similarity_to_distance(similarity_df):
    """
    Converts a similarity matrix to a distance matrix, accounting for HLL approximations.

    :param similarity_df: pandas DataFrame of the similarity matrix.
    :return: pandas DataFrame of the distance matrix.
    """
    # Ensure all similarities are between 0 and 1
    similarity_df = similarity_df.clip(lower=0, upper=1)
    
    # Add a small epsilon to avoid log(0)
    epsilon = 1e-10
    similarity_df += epsilon
    
    # Convert similarity to distance using -log(similarity)
    distance_df = -np.log(similarity_df)
    
    # Since HLL estimates are approximate, consider scaling or smoothing if necessary
    return distance_df

# Step 3: Construct the Phylogenetic Tree using Neighbor-Joining
def construct_tree(distance_df):
    """
    Constructs a phylogenetic tree using the Neighbor-Joining method.

    :param distance_df: pandas DataFrame of the distance matrix.
    :return: skbio TreeNode representing the phylogenetic tree.
    """
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
def main(similarity_csv_path, output_tree_path):
    # Read the similarity matrix
    similarity_df = read_similarity_csv(similarity_csv_path)
    print("Similarity matrix loaded.")

    # Convert similarity to distance
    distance_df = similarity_to_distance(similarity_df)
    print("Similarity matrix converted to distance matrix.")

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
    similarity_csv_path = "similarity_matrix.csv"  # Replace with your CSV file path
    output_tree_path = "phylogenetic_tree.nwk"     # Output Newick file path
    main(similarity_csv_path, output_tree_path)