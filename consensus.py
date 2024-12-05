# Import necessary modules
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from collections import Counter
import os
import matplotlib.pyplot as plt

# File paths
tree_files = [
    '/users/harry/desktop/Computational Genetics/Final project/code/result/NN/Consensus_Phylogenetic_Tree.nex',
    '/users/harry/desktop/Computational Genetics/Final project/code/result/UMPGA/Consensus_UPGMA_Tree.nex'
]

# Function to extract splits from a tree
def get_splits(tree):
    splits = []
    clades = list(tree.find_clades(order='level'))

    # For each internal node (excluding the root), get the set of taxa under it
    for clade in clades:
        if clade.is_terminal():
            continue  # Skip leaf nodes
        taxa = clade.get_terminals()
        taxa_names = frozenset(taxon.name for taxon in taxa)
        splits.append(taxa_names)
    return splits

# Read trees and collect splits
split_counts = Counter()
all_taxa = set()
num_trees = len(tree_files)

for file_path in tree_files:
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        continue

    # Read the tree from the Nexus file
    trees = list(Phylo.parse(file_path, 'nexus'))
    if not trees:
        print(f"No trees found in file: {file_path}")
        continue

    # Assuming one tree per file
    tree = trees[0]
    # Collect all taxa
    all_taxa.update(taxon.name for taxon in tree.get_terminals())
    # Extract splits
    splits = get_splits(tree)
    # Update split counts
    split_counts.update(splits)

# Determine majority splits (appear in more than 50% of trees)
majority_threshold = num_trees / 2
majority_splits = {split for split, count in split_counts.items() if count > majority_threshold}

# Function to build consensus tree from majority splits
def build_consensus_tree(all_taxa, majority_splits):
    # Start with all taxa in an unresolved polytomy
    consensus_tree = Tree()
    consensus_tree.root.clades = [Phylo.Newick.Clade(name=taxon) for taxon in all_taxa]

    # Iteratively resolve nodes based on majority splits
    changed = True
    while changed:
        changed = False
        for split in majority_splits.copy():
            # Find the minimal clade that contains all taxa in the split
            for clade in consensus_tree.get_nonterminals(order='postorder'):
                clade_taxa = {taxon.name for taxon in clade.get_terminals()}
                if split < clade_taxa and split != clade_taxa:
                    # Split the clade into two: taxa in the split and taxa not in the split
                    taxa_in_split = []
                    taxa_not_in_split = []
                    for subclade in clade.clades:
                        subclade_taxa = {taxon.name for taxon in subclade.get_terminals()}
                        if subclade_taxa & split:
                            taxa_in_split.append(subclade)
                        else:
                            taxa_not_in_split.append(subclade)
                    if taxa_in_split and taxa_not_in_split:
                        # Create new clade for taxa in split
                        new_clade = Phylo.Newick.Clade()
                        new_clade.clades = taxa_in_split
                        # Update the current clade
                        clade.clades = [new_clade] + taxa_not_in_split
                        changed = True
                    majority_splits.remove(split)
                    break
    return consensus_tree

# Build the consensus tree
consensus_tree = build_consensus_tree(all_taxa, majority_splits)

# Save the consensus tree to a Nexus file
output_file = '/users/harry/desktop/Computational Genetics/Final project/code/result/Consensus_Majority_Rule_Tree.nex'
Phylo.write(consensus_tree, output_file, 'nexus')

print(f"Consensus tree saved to {output_file}")

# Visualization
plt.figure(figsize=(10, 8))
Phylo.draw(consensus_tree, do_show=False)
plt.title('Majority-Rule Consensus Tree')
plt.show()

