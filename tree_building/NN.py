import pandas as pd
import numpy as np
from skbio import DistanceMatrix
from skbio.tree import nj
from dendropy import Tree as DendroPyTree
from dendropy import TaxonNamespace
import glob
import matplotlib.pyplot as plt
from io import StringIO
from Bio import Phylo
import os
import dendropy
import io

def read_distance_csv(filePath):
    """Reads a CSV file containing the distance matrix."""
    distanceDf = pd.read_csv(filePath, index_col=0, na_values=['-'])
    distanceDf = distanceDf.astype(float)
    return distanceDf

def construct_nj_tree(distanceDf):
    """Constructs a Neighbor-Joining tree from a distance matrix."""
    ids = distanceDf.index.tolist()
    distanceArray = distanceDf.values
    dm = DistanceMatrix(distanceArray, ids)
    tree = nj(dm)
    return tree

def skbio_to_dendropy_tree(skbioTree):
    """Converts a skbio tree to a dendropy tree."""
    newickIo = io.StringIO()
    skbioTree.write(newickIo, format='newick')
    newickStr = newickIo.getvalue()
    newickIo.close()
    dendropyTree = DendroPyTree.get(data=newickStr, schema="newick", taxon_namespace=TaxonNamespace())
    return dendropyTree

def visualize_tree(dendropyTree, outputImagePath, treeLabel):
    """Visualizes the tree and saves it as an image."""
    newickStr = dendropyTree.as_string(schema='newick')
    handle = StringIO(newickStr)
    phyloTree = Phylo.read(handle, 'newick')
    fig = plt.figure(figsize=(12, 8))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(phyloTree, do_show=False, axes=axes)
    plt.title(f"Phylogenetic Tree: {treeLabel}", fontsize=16)
    plt.savefig(outputImagePath, format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Phylogenetic tree visualization saved to '{outputImagePath}'.")

def save_trees_nexus(dendropyTrees, outputTreePath):
    """Saves all trees in Nexus format."""
    with open(outputTreePath, 'w') as nexusFile:
        nexusFile.write("#NEXUS\n")
        nexusFile.write("BEGIN TREES;\n")
        for idx, (tree, label) in enumerate(dendropyTrees, 1):
            nexusFile.write(f"    TREE {label} = {tree.as_string(schema='newick').strip()}\n")
        nexusFile.write("END;\n")
    print(f"All phylogenetic trees saved to '{outputTreePath}' in Nexus format.")

def process_single_csv(filePath, outputDir):
    """Processes a single CSV file containing a distance matrix."""
    baseName = os.path.basename(filePath)
    treeLabel = os.path.splitext(baseName)[0]
    print(f"\nProcessing file: '{filePath}' with label '{treeLabel}'.")
    distanceDf = read_distance_csv(filePath)
    print(" - Loaded distance matrix.")
    njTreeSkbio = construct_nj_tree(distanceDf)
    print(" - Constructed Neighbor-Joining tree using skbio.")
    dendropyTree = skbio_to_dendropy_tree(njTreeSkbio)
    print(" - Converted skbio tree to DendroPy format.")
    outputImagePath = os.path.join(outputDir, f"{treeLabel}_Phylogenetic_Tree.png")
    visualize_tree(dendropyTree, outputImagePath, treeLabel)
    return dendropyTree, treeLabel

def NN(generalFolderPath):
    """Main function to process all CSV files and generate trees."""
    inputDir = os.path.join(generalFolderPath, 'input')
    outputDir = os.path.join(generalFolderPath, 'result', 'NN')
    os.makedirs(outputDir, exist_ok=True)
    distanceCsvFiles = glob.glob(os.path.join(inputDir, '*.csv'))
    if not distanceCsvFiles:
        print("No CSV files found in the specified directory. Please check the file paths.")
    else:
        dendropyTrees = []
        for filePath in distanceCsvFiles:
            dendropyTree, treeLabel = process_single_csv(filePath, outputDir)
            dendropyTrees.append((dendropyTree, treeLabel))
            print(f"Completed processing for '{treeLabel}'.")
        outputTreePath = os.path.join(outputDir, "All_Phylogenetic_Trees.nex")
        save_trees_nexus(dendropyTrees, outputTreePath)
        print("\nAll phylogenetic tree generation and saving completed successfully.")

if __name__ == "__main__":
    generalFolderPath = '/users/harry/desktop/Computational_Genetics/Final_project/code/'
    NN(generalFolderPath)
