import pandas as pd
import numpy as np
from skbio import DistanceMatrix
from skbio.tree import nj

def read_similarity_csvs(file_paths):
    similarity_dfs = []
    for file_path in file_paths:
        similarity_df = pd.read_csv(file_path, index_col=0)
        similarity_dfs.append(similarity_df)
    return similarity_dfs

def similarity_to_distance(similarity_dfs):
    distance_dfs = []
    for similarity_df in similarity_dfs:
        # Ensure all similarities are between 0 and 1
        similarity_df = similarity_df.clip(lower=0, upper=1)
        epsilon = 1e-10
        similarity_df += epsilon
        distance_df = -np.log(similarity_df)

        # Set the diagonal to zero for a valid distance matrix
        np.fill_diagonal(distance_df.values, 0)
        distance_dfs.append(distance_df)
    return distance_dfs

def compute_consensus_distance(distance_dfs):
    # Stack the distance matrices along a new axis
    stacked_distances = np.stack([df.values for df in distance_dfs], axis=2)
    consensus_distances = np.mean(stacked_distances, axis=2)
    ids = distance_dfs[0].index.tolist()
    consensus_df = pd.DataFrame(consensus_distances, index=ids, columns=ids)
    return consensus_df

def construct_tree(distance_df):
    ids = distance_df.index.tolist()
    distance_array = distance_df.values
    dm = DistanceMatrix(distance_array, ids)
    tree = nj(dm)
    return tree

def main(similarity_csv_paths, output_tree_path):

    similarity_dfs = read_similarity_csvs(similarity_csv_paths)
    print(f"{len(similarity_dfs)} similarity matrices loaded.")

    distance_dfs = similarity_to_distance(similarity_dfs)
    print("Similarity matrices converted to distance matrices.")

    consensus_distance_df = compute_consensus_distance(distance_dfs)
    print("Consensus distance matrix computed by averaging.")

    tree = construct_tree(consensus_distance_df)
    print("Phylogenetic tree constructed using Neighbor-Joining method.")

    with open(output_tree_path, 'w') as f:
        f.write(tree.write("/users/harry/desktop/Computational Genetics/Final project/phylogenetic_tree.nwk") + ";")
    print(f"Phylogenetic tree saved to {output_tree_path}.")

    print("\nConstructed Phylogenetic Tree:")
    print(tree.ascii_art())

if __name__ == "__main__":
    similarity_csv_paths = [
        "/users/harry/desktop/Computational Genetics/Final project/distance_matrix_genome.csv",
        "/users/harry/desktop/Computational Genetics/Final project/cds_similarity_matrix.csv",
    ]
    output_tree_path = "/users/harry/desktop/Computational Genetics/Final project/phylogenetic_tree.nwk"
    main(similarity_csv_paths, output_tree_path)
