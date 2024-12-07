import re
import os

def convert_to_unweighted_newick(weighted_newick):
    unweighted_newick = re.sub(r':\d+(\.\d+)?', '', weighted_newick)
    return unweighted_newick

result_path = os.path.join("test","result", "Consensus")
# paths = ["Consensus_Majority_Rule_Tree_cds.nwk", "Consensus_Majority_Rule_Tree_gene.nwk", "Consensus_Majority_Rule_Tree_whole_genome.nwk"]
# for p in paths:
#     with open(os.path.join(result_path, p), 'r') as file:
#         weighted_newick_tree = file.readline()
#     print(p)
#     unweighted_newick_tree = convert_to_unweighted_newick(weighted_newick_tree)
#     print(unweighted_newick_tree)
#     print()
weighted_newick_tree = "((((CB9615:0.1,(Sakai:0.1,EDL933:0.1):0.1):0.1,((E24377A:0.1,(IAI1:0.1,SE11:0.1):0.1):0.2,(HS:0.2,((BW2952:0.1,DH10B:0.1):0.1,(MG1655:0.1,W3110:0.1):0.1):0.1):0.1):0.1):0.1,UMN026:0.1):0.1,((IAI39:0.1,SMS35:0.1):0.1,(E234869:0.1,(536:0.1,((S88:0.1,(APEC01:0.1,UTI89:0.1):0.1):0.1,(ED1a:0.1,CFT073:0.1):0.1):0.1):0.1):0.1):0.1):0.1;"
print(convert_to_unweighted_newick(weighted_newick_tree))
