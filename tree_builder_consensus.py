import pandas as pd
import numpy as np
from skbio import DistanceMatrix
from skbio.tree import nj

# Step 1: Read Multiple Similarity Matrices from CSV Files
def read_similarity_csvs(file_paths):
    """
    Reads multiple CSV files containing similarity matrices.

    :param file_paths: List of paths to the CSV files.
    :return: List of pandas DataFrames of the similarity matrices.
    """
    similarity_dfs = []
    for file_path in file_paths:
        similarity_df = pd.read_csv(file_path, index_col=0)
        similarity_dfs.append(similarity_df)
    return similarity_dfs

# Step 2: Convert Similarity Matrices to Distance Matrices
def similarity_to_distance(similarity_dfs):
    """
    Converts a list of similarity matrices to distance matrices.

    :param similarity_dfs: List of pandas DataFrames of similarity matrices.
    :return: List of pandas DataFrames of distance matrices.
    """
    distance_dfs = []
    for similarity_df in similarity_dfs:
        # Ensure all similarities are between 0 and 1
        similarity_df = similarity_df.clip(lower=0, upper=1)

        # Add a small epsilon to avoid log(0)
        epsilon = 1e-10
        similarity_df += epsilon

        # Convert similarity to distance using -log(similarity)
        distance_df = -np.log(similarity_df)

        # Set the diagonal to zero for a valid distance matrix
        np.fill_diagonal(distance_df.values, 0)

        distance_dfs.append(distance_df)
    return distance_dfs

# Step 3: Compute the Consensus Distance Matrix
def compute_consensus_distance(distance_dfs):
    """
    Computes the consensus distance matrix by averaging the distances.

    :param distance_dfs: List of pandas DataFrames of distance matrices.
    :return: pandas DataFrame of the consensus distance matrix.
    """
    # Stack the distance matrices along a new axis
    stacked_distances = np.stack([df.values for df in distance_dfs], axis=2)

    # Compute the mean across the stacked distance matrices
    consensus_distances = np.mean(stacked_distances, axis=2)

    # Create a DataFrame with the consensus distances
    ids = distance_dfs[0].index.tolist()
    consensus_df = pd.DataFrame(consensus_distances, index=ids, columns=ids)

    return consensus_df

# Step 4: Construct the Phylogenetic Tree using Neighbor-Joining
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

# Step 5: Main Function to Tie Everything Together
def main(similarity_csv_paths, output_tree_path):
    # Read the similarity matrices
    similarity_dfs = read_similarity_csvs(similarity_csv_paths)
    print(f"{len(similarity_dfs)} similarity matrices loaded.")

    # Convert similarities to distances
    distance_dfs = similarity_to_distance(similarity_dfs)
    print("Similarity matrices converted to distance matrices.")

    # Compute the consensus distance matrix
    consensus_distance_df = compute_consensus_distance(distance_dfs)
    print("Consensus distance matrix computed by averaging.")

    # Construct the phylogenetic tree
    tree = construct_tree(consensus_distance_df)
    print("Phylogenetic tree constructed using Neighbor-Joining method.")

    # Write the tree to a file in Newick format
    with open(output_tree_path, 'w') as f:
        f.write(tree.write("/users/harry/desktop/Computational Genetics/Final project/phylogenetic_tree.nwk") + ";")
    print(f"Phylogenetic tree saved to {output_tree_path}.")

    # Optional: Print the tree to the console
    print("\nConstructed Phylogenetic Tree:")
    print(tree.ascii_art())

# Example usage
if __name__ == "__main__":
    # List of similarity matrix CSV file paths obtained through resampling
    similarity_csv_paths = [
        "/users/harry/desktop/Computational Genetics/Final project/similarity_matrix1.csv",
        "/users/harry/desktop/Computational Genetics/Final project/similarity_matrix2.csv",
        # Add more file paths as needed
    ]

    output_tree_path = "/users/harry/desktop/Computational Genetics/Final project/phylogenetic_tree.nwk"
    main(similarity_csv_paths, output_tree_path)
