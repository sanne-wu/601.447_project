import dendropy
from dendropy.calculate import consensus

tree_files = [
    '/users/harry/desktop/Computational Genetics/Final project/code/result/NN/Consensus_Phylogenetic_Tree.nex',
    'users/harry/desktop/Computational Genetics/Final project/code/result/UMPGA/Consensus_UPGMA_Tree.nex',
]

# Define a common taxon namespace
taxon_namespace = dendropy.TaxonNamespace()

# Read all trees into a TreeList
trees = dendropy.TreeList()
for file in tree_files:
    tree = dendropy.Tree.get(
        path=file,
        schema="newick",
        taxon_namespace=taxon_namespace
    )
    trees.append(tree)

# Compute the majority-rule consensus tree
# You can adjust the consensus method as needed
consensus_tree = consensus.majority_rule_consensus(trees, majority=0.5)

# Alternatively, for strict consensus:
# consensus_tree = consensus.strict_consensus(trees)

# Print the consensus tree in Newick format
print(consensus_tree.as_string(schema="newick"))

# Optionally, save the consensus tree to a file
consensus_tree.write(
    path="consensus_tree.nwk",
    schema="newick",
    suppress_rooting=True
)
