import pandas as pd
import numpy as np
from skbio import DistanceMatrix
from skbio.tree import nj
from dendropy import Tree, TaxonNamespace
from dendropy.calculate import treecompare
import glob
import matplotlib.pyplot as plt
from io import StringIO
from Bio import Phylo
import io
import dendropy

# Step 1: Read the Similarity Matrix from CSV
def read_similarity_csv(file_path):
    # Reads a CSV file containing the similarity matrix.
    similarity_df = pd.read_csv(file_path, index_col=0)
    return similarity_df

# Step 2: Convert Similarity Matrix to Distance Matrix
def similarity_to_distance(similarity_df):
    # Ensure all similarities are between 0 and 1
    similarity_df = similarity_df.clip(lower=0, upper=1)
    # Add a small epsilon to avoid log(0)
    epsilon = 1e-10
    similarity_df += epsilon
    # Convert similarity to distance using -log(similarity)
    distance_df = -np.log(similarity_df)
    # Set the diagonal to zero for a valid distance matrix
    np.fill_diagonal(distance_df.values, 0)
    return distance_df

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

# Step 4: Convert scikit-bio TreeNode to DendroPy Tree
def skbio_to_dendropy_tree(skbio_tree):
    # Use StringIO to capture the Newick output as a string
    newick_io = io.StringIO()
    skbio_tree.write(newick_io, format='newick')
    newick_str = newick_io.getvalue()
    newick_io.close()
    
    # Convert Newick string to a DendroPy Tree
    dendropy_tree = dendropy.Tree.get(data=newick_str, schema="newick")
    return dendropy_tree

# Step 5: Calculate Bootstrap Support
def calculate_bootstrap_support(original_tree, bootstrap_trees):

    if not bootstrap_trees:
        raise ValueError("The list of bootstrap trees is empty. Please provide bootstrap trees for support calculation.")
    
    # Initialize support values for each node
    for node in original_tree.postorder_node_iter():
        node.bootstrap_support = 0

    # Encode bipartitions in the original tree and initialize split counts
    original_tree.encode_bipartitions()
    split_counts = {bipartition: 0 for bipartition in original_tree.bipartition_edge_map}

    # Count how many times each bipartition appears in the bootstrap trees
    for bt in bootstrap_trees:
        bt.encode_bipartitions()
        for bipartition in bt.bipartition_edge_map:
            if bipartition in split_counts:
                split_counts[bipartition] += 1

    # Calculate support values and assign to nodes
    num_bootstrap_trees = len(bootstrap_trees)
    for edge in original_tree.postorder_edge_iter():
        if edge.head_node is not None and not edge.head_node.is_leaf():
            bipartition = edge.bipartition
            support = (split_counts.get(bipartition, 0) / num_bootstrap_trees) * 100
            edge.head_node.label = f"{support:.1f}%"

    return original_tree

# Step 6: Visualize and Save the Annotated Tree
def visualize_tree(dendropy_tree, output_image_path):
    # Convert DendroPy tree to Newick string
    newick_str = dendropy_tree.as_string(schema='newick')

    # Read the Newick string into a Biopython Phylo tree
    from Bio import Phylo
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

# Step 7: Main Function to Tie Everything Together
def main(original_similarity_csv, bootstrap_similarity_dir, output_tree_path, output_image_path):
    # Read the original similarity matrix
    similarity_df = read_similarity_csv(original_similarity_csv)
    print("Original similarity matrix loaded.")

    # Convert similarity to distance
    distance_df = similarity_to_distance(similarity_df)
    print("Original similarity matrix converted to distance matrix.")

    # Construct the original phylogenetic tree
    original_skbio_tree = construct_tree(distance_df)
    print("Original phylogenetic tree constructed using Neighbor-Joining method.")

    # Convert to DendroPy tree
    original_tree = skbio_to_dendropy_tree(original_skbio_tree)

    # Read bootstrap similarity matrices and construct trees
    bootstrap_trees = []
    bootstrap_files = glob.glob(f"{bootstrap_similarity_dir}/*.csv")
    print(f"Found {len(bootstrap_files)} bootstrap similarity matrices.")

    for i, file_path in enumerate(bootstrap_files):
        # Read the similarity matrix
        sim_df = read_similarity_csv(file_path)
        # Convert similarity to distance
        dist_df = similarity_to_distance(sim_df)
        # Construct the tree
        skbio_tree = construct_tree(dist_df)
        # Convert to DendroPy tree
        dendropy_tree = skbio_to_dendropy_tree(skbio_tree)
        bootstrap_trees.append(dendropy_tree)
        print(f"Processed bootstrap tree {i+1}/{len(bootstrap_files)}")

    # Calculate bootstrap support values
    annotated_tree = calculate_bootstrap_support(original_tree, bootstrap_trees)
    print("Bootstrap support values calculated and annotated on the original tree.")

    # Write the annotated tree to a file in Newick format
    annotated_tree.write(path=output_tree_path, schema='newick')
    print(f"Annotated phylogenetic tree saved to {output_tree_path}.")

    # Visualize and save the annotated tree
    visualize_tree(annotated_tree, output_image_path)

# Example usage
if __name__ == "__main__":
    original_similarity_csv = "/users/harry/desktop/Computational Gentics/Final project/similarity_matrix.csv"  # Path to the original similarity matrix
    bootstrap_similarity_dir = "/users/harry/desktop/Computational Gentics/Final project/bootstrap_matrices"    # Directory containing bootstrap similarity matrices
    output_tree_path = "/users/harry/desktop/Computational Gentics/Final project/annotated_phylogenetic_tree.nwk"  # Output Newick file path
    output_image_path = "/users/harry/desktop/Computational Gentics/Final project/annotated_phylogenetic_tree.png" # Output image file path
    main(original_similarity_csv, bootstrap_similarity_dir, output_tree_path, output_image_path)