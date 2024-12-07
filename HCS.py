import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
from io import StringIO
from Bio import Phylo
import os
import glob
import dendropy
from dendropy import Tree as DendroPyTree
from dendropy import TaxonNamespace
from scipy.cluster.hierarchy import to_tree

def readDistanceCsv(filePath):
    distanceDf = pd.read_csv(filePath, index_col=0, na_values=['-'])
    distanceDf = distanceDf.astype(float)
    return distanceDf

def fitchMargoliash(distanceMatrix):
    n = distanceMatrix.shape[0]
    condensedDistMatrix = sch.distance.squareform(distanceMatrix)
    linkageMatrix = sch.linkage(condensedDistMatrix, method='single')
    return linkageMatrix

def linkageToNewick(linkageMatrix, labels):
    tree = to_tree(linkageMatrix, rd=False)

    def buildNewick(node):
        if node.is_leaf():
            return labels[node.id]
        else:
            left = buildNewick(node.left)
            right = buildNewick(node.right)
            leftLength = node.dist - node.left.dist
            rightLength = node.dist - node.right.dist
            return f"({left}:{leftLength:.2f},{right}:{rightLength:.2f})"
    
    return buildNewick(tree) + ";"

def visualizeTree(linkageMatrix, outputImagePath, treeLabel, labels):
    plt.figure(figsize=(12, 8))
    plt.title(f"Phylogenetic Tree: {treeLabel}", fontsize=16)
    plt.xlabel("Distance")
    plt.ylabel("Taxa")
    sch.dendrogram(linkageMatrix, labels=labels, leaf_rotation=90)
    plt.tight_layout()
    plt.savefig(outputImagePath, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Phylogenetic tree visualization saved to '{outputImagePath}'.")

def saveTreesNexus(newickTrees, outputTreePath):
    with open(outputTreePath, 'w') as nexusFile:
        nexusFile.write("#NEXUS\n")
        nexusFile.write("BEGIN TREES;\n")
        for idx, (newickStr, label) in enumerate(newickTrees, 1):
            nexusFile.write(f"    TREE {label} = {newickStr}\n")
        nexusFile.write("END;\n")
    print(f"All phylogenetic trees saved to '{outputTreePath}' in Nexus format.")

def processSingleCsv(filePath, outputDir):
    baseName = os.path.basename(filePath)
    treeLabel = os.path.splitext(baseName)[0]
    print(f"\nProcessing file: '{filePath}' with label '{treeLabel}'.")
    distanceDf = readDistanceCsv(filePath)
    print(" - Loaded distance matrix.")
    linkageMatrix = fitchMargoliash(distanceDf.values)
    print(" - Constructed phylogenetic tree using Fitch-Margoliash method.")
    labels = distanceDf.index.tolist()
    outputImagePath = os.path.join(outputDir, f"{treeLabel}_Phylogenetic_Tree.png")
    visualizeTree(linkageMatrix, outputImagePath, treeLabel, labels)
    newickStr = linkageToNewick(linkageMatrix, labels)
    dendroTree = DendroPyTree.get(data=newickStr, schema="newick", taxon_namespace=TaxonNamespace())
    return dendroTree, treeLabel

def FM(generalFolderPath):
    inputDir = os.path.join(generalFolderPath, 'input')
    outputDir = os.path.join(generalFolderPath, 'result', 'FM')
    os.makedirs(outputDir, exist_ok=True)
    distanceCsvFiles = glob.glob(os.path.join(inputDir, '*.csv'))
    if not distanceCsvFiles:
        print("No CSV files found in the specified directory. Please check the file paths.")
    else:
        newickTrees = []
        for filePath in distanceCsvFiles:
            try:
                dendroTree, treeLabel = processSingleCsv(filePath, outputDir)
                newickStr = dendroTree.as_string(schema='newick').strip()
                newickTrees.append((newickStr, treeLabel))
                print(f"Completed processing for '{treeLabel}'.")
            except Exception as e:
                print(f"Error processing file '{filePath}': {e}")
        
        if newickTrees:
            outputTreePath = os.path.join(outputDir, "All_Phylogenetic_Trees.nex")
            saveTreesNexus(newickTrees, outputTreePath)
            print("\nAll phylogenetic tree generation and saving completed successfully.")
        else:
            print("No trees were processed successfully.")

if __name__ == "__main__":
    generalFolderPath = '/users/harry/desktop/Computational_Genetics/Final_project/code/'
    FM(generalFolderPath)
