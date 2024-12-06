import pandas as pd
import numpy as np
from skbio import DistanceMatrix
from dendropy import Tree as DendroPyTree
from dendropy import TaxonNamespace
import glob
import matplotlib.pyplot as plt
from io import StringIO
from Bio import Phylo
import os
import dendropy
import io
import subprocess
import tempfile

def read_similarity_csv(filePath):
    similarityDf = pd.read_csv(filePath, index_col=0, na_values=['-'])
    similarityDf = similarityDf.astype(float)
    return similarityDf

def similarity_to_distance(similarityDf):
    similarityDf = similarityDf.clip(lower=0, upper=1)
    epsilon = 1e-10
    similarityDf += epsilon
    distanceDf = -np.log(similarityDf)
    np.fill_diagonal(distanceDf.values, 0)
    return distanceDf

def run_minimal_evolution(distanceDf):
    # Convert the distance matrix DataFrame into PHYLIP-like format for RapidNJ
    ids = distanceDf.index.tolist()
    n = len(ids)
    
    # RapidNJ expects a PHYLIP-formatted distance file, where:
    # First line: number of taxa
    # Following lines: TaxonName Distances...
    # 
    # Example:
    #  4
    #  Taxon1 0.0  0.1  0.2  0.3
    #  Taxon2 0.1  0.0  0.2  0.2
    #  ...
    
    with tempfile.NamedTemporaryFile('w', delete=False) as dist_file:
        dist_file_name = dist_file.name
        dist_file.write(f"{n}\n")
        for i, taxon in enumerate(ids):
            dist_line = " ".join(f"{d:.6f}" for d in distanceDf.iloc[i])
            dist_file.write(f"{taxon} {dist_line}\n")

    # Run RapidNJ with minimal evolution method
    # -i pd : input is a phylip distance file
    # -o t  : output is a tree in newick format
    # -m ME : minimal evolution method
    # 
    # RapidNJ usage example:
    # rapidnj input_dist_file -i pd -o t -m ME > output_tree.txt
    
    with tempfile.NamedTemporaryFile('w', delete=False) as tree_file:
        tree_file_name = tree_file.name

    cmd = ["rapidnj", dist_file_name, "-i", "pd", "-o", "t", "-m", "ME"]
    with open(tree_file_name, 'w') as outfile:
        subprocess.run(cmd, stdout=outfile, check=True)

    # Read the resulting Newick tree
    with open(tree_file_name, 'r') as f:
        newickStr = f.read().strip()

    # Cleanup temporary files
    os.remove(dist_file_name)
    os.remove(tree_file_name)
    
    # Convert to DendroPy tree
    dendropyTree = DendroPyTree.get(data=newickStr, schema="newick", taxon_namespace=TaxonNamespace())
    return dendropyTree

def visualize_tree(dendropyTree, outputImagePath, treeLabel):
    newickStr = dendropyTree.as_string(schema='newick')
    handle = StringIO(newickStr)
    phyloTree = Phylo.read(handle, 'newick')
    fig = plt.figure(figsize=(12, 8))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(phyloTree, do_show=False, axes=axes)
    plt.title(f"Phylogenetic Tree (ME): {treeLabel}", fontsize=16)
    plt.savefig(outputImagePath, format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Phylogenetic tree visualization saved to '{outputImagePath}'.")

def save_trees_nexus(dendropyTrees, outputTreePath):
    with open(outputTreePath, 'w') as nexusFile:
        nexusFile.write("#NEXUS\n")
        nexusFile.write("BEGIN TREES;\n")
        for idx, (tree, label) in enumerate(dendropyTrees, 1):
            nexusFile.write(f"    TREE {label} = {tree.as_string(schema='newick').strip()}\n")
        nexusFile.write("END;\n")
    print(f"All phylogenetic trees saved to '{outputTreePath}' in Nexus format.")

def process_single_csv(filePath, outputDir):
    baseName = os.path.basename(filePath)
    treeLabel = os.path.splitext(baseName)[0]
    print(f"\nProcessing file: '{filePath}' with label '{treeLabel}'.")
    similarityDf = read_similarity_csv(filePath)
    print(" - Loaded similarity matrix.")
    distanceDf = similarity_to_distance(similarityDf)
    print(" - Converted similarity matrix to distance matrix.")
    dendropyTree = run_minimal_evolution(distanceDf)
    print(" - Constructed Minimal Evolution tree using RapidNJ.")
    outputImagePath = os.path.join(outputDir, f"{treeLabel}_ME_Phylogenetic_Tree.png")
    visualize_tree(dendropyTree, outputImagePath, treeLabel)
    return dendropyTree, treeLabel

def ME_pipeline(generalFolderPath):
    inputDir = os.path.join(generalFolderPath, 'input')
    outputDir = os.path.join(generalFolderPath, 'result', 'ME')
    os.makedirs(outputDir, exist_ok=True)
    similarityCsvFiles = glob.glob(os.path.join(inputDir, '*.csv'))
    if not similarityCsvFiles:
        print("No CSV files found in the specified directory. Please check the file paths.")
    else:
        dendropyTrees = []
        for filePath in similarityCsvFiles:
            dendropyTree, treeLabel = process_single_csv(filePath, outputDir)
            dendropyTrees.append((dendropyTree, treeLabel))
            print(f"Completed processing for '{treeLabel}'.")
        outputTreePath = os.path.join(outputDir, "All_Phylogenetic_Trees_ME.nex")
        save_trees_nexus(dendropyTrees, outputTreePath)
        print("\nAll phylogenetic tree generation and saving completed successfully (Minimal Evolution).")

if __name__ == "__main__":
    generalFolderPath = '/users/harry/desktop/Computational_Genetics/Final_project/code/'
    ME_pipeline(generalFolderPath)
