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

def construct_nj_tree(distanceDf):
    ids = distanceDf.index.tolist()
    distanceArray = distanceDf.values
    dm = DistanceMatrix(distanceArray, ids)
    tree = nj(dm)
    return tree

def skbio_to_dendropy_tree(skbioTree):
    newickIo = io.StringIO()
    skbioTree.write(newickIo, format='newick')
    newickStr = newickIo.getvalue()
    newickIo.close()
    dendropyTree = DendroPyTree.get(data=newickStr, schema="newick", taxon_namespace=TaxonNamespace())
    return dendropyTree

def visualize_tree(dendropyTree, outputImagePath, treeLabel):
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
    njTreeSkbio = construct_nj_tree(distanceDf)
    print(" - Constructed Neighbor-Joining tree using skbio.")
    dendropyTree = skbio_to_dendropy_tree(njTreeSkbio)
    print(" - Converted skbio tree to DendroPy format.")
    outputImagePath = os.path.join(outputDir, f"{treeLabel}_Phylogenetic_Tree.png")
    visualize_tree(dendropyTree, outputImagePath, treeLabel)
    return dendropyTree, treeLabel

def NN(generalFolderPath):
    inputDir = os.path.join(generalFolderPath, 'input')
    outputDir = os.path.join(generalFolderPath, 'result', 'NN')
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
        outputTreePath = os.path.join(outputDir, "All_Phylogenetic_Trees.nex")
        save_trees_nexus(dendropyTrees, outputTreePath)
        print("\nAll phylogenetic tree generation and saving completed successfully.")

if __name__ == "__main__":
    generalFolderPath = '/users/harry/desktop/Computational Genetics/Final project/code/'
    NN(generalFolderPath)
