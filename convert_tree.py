import re
import os

def convert_to_unweighted_newick(weighted_newick):
    unweighted_newick = re.sub(r':\d+(\.\d+)?', '', weighted_newick)
    return unweighted_newick

result_path = os.path.join("test","result", "Consensus")
paths = ["Consensus_Majority_Rule_Tree_cds.nwk", "Consensus_Majority_Rule_Tree_gene.nwk", "Consensus_Majority_Rule_Tree_whole_genome.nwk"]
for p in paths:
    with open(os.path.join(result_path, p), 'r') as file:
        weighted_newick_tree = file.readline()
    print(p)
    unweighted_newick_tree = convert_to_unweighted_newick(weighted_newick_tree)
    print(unweighted_newick_tree)
    print()
