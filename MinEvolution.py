import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import glob
import os
import matplotlib.pyplot as plt

def readDistanceMatrix(filePath):
    distanceDf = pd.read_csv(filePath, index_col=0, na_values=['-'])
    distanceDf = distanceDf.astype(float)
    return distanceDf

def temporarilyMaskMatrix(distanceDf):
    originalDf = distanceDf.copy()
    mask = np.triu(np.ones(distanceDf.shape), k=1).astype(bool)
    distanceDf.where(~mask, other=0, inplace=True)
    return distanceDf, originalDf

def restoreMatrix(distanceDf, originalDf):
    distanceDf[:] = originalDf[:]

def constructMinEvolutionTree(distanceDf):
    distanceDf, originalDf = temporarilyMaskMatrix(distanceDf)
    print("Masked Distance Matrix (Lower Triangle Only):")
    print(distanceDf)
    ids = distanceDf.index.tolist()
    print("IDs (Taxon Names):", ids)
    print("Number of taxa (matrix size):", len(ids))
    n = len(ids)
    matrix = []
    for i in range(n):
        row = []
        for j in range(i):
            row.append(distanceDf.iloc[i, j])
        row.append(0)
        matrix.append(row)
    print("Matrix size (rows):", len(matrix))
    print("First few rows of the matrix:")
    for row in matrix[:9]:
        print(row)
    distanceMatrix = DistanceMatrix(names=ids, matrix=matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distanceMatrix)
    restoreMatrix(distanceDf, originalDf)
    return tree

def visualizeTree(tree, outputImagePath, treeLabel):
    fig = plt.figure(figsize=(12, 8))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=axes)
    plt.title(f"Phylogenetic Tree: {treeLabel}", fontsize=16)
    plt.savefig(outputImagePath, format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Phylogenetic tree visualization saved to '{outputImagePath}'.")

def saveTreesNexus(trees, outputTreePath):
    with open(outputTreePath, 'w') as nexusFile:
        nexusFile.write("#NEXUS\n")
        nexusFile.write("BEGIN TREES;\n")
        for idx, (tree, label) in enumerate(trees, 1):
            nexusFile.write(f"    TREE {label} = {tree.format('newick').strip()}\n")
        nexusFile.write("END;\n")
    print(f"All phylogenetic trees saved to '{outputTreePath}' in Nexus format.")

def processSingleCsv(filePath, outputDir):
    baseName = os.path.basename(filePath)
    treeLabel = os.path.splitext(baseName)[0]
    print(f"\nProcessing file: '{filePath}' with label '{treeLabel}'.")
    distanceDf = readDistanceMatrix(filePath)
    print(" - Loaded distance matrix.")
    minEvoTree = constructMinEvolutionTree(distanceDf)
    print(" - Constructed Minimum Evolution tree.")
    outputImagePath = os.path.join(outputDir, f"{treeLabel}_Phylogenetic_Tree.png")
    visualizeTree(minEvoTree, outputImagePath, treeLabel)
    return minEvoTree, treeLabel

def minEvolution(generalFolderPath):
    inputDir = os.path.join(generalFolderPath, 'input')
    outputDir = os.path.join(generalFolderPath, 'result', 'MinE')
    os.makedirs(outputDir, exist_ok=True)
    distanceCsvFiles = glob.glob(os.path.join(inputDir, '*.csv'))
    if not distanceCsvFiles:
        print("No CSV files found in the specified directory. Please check the file paths.")
    else:
        trees = []
        for filePath in distanceCsvFiles:
            tree, treeLabel = processSingleCsv(filePath, outputDir)
            trees.append((tree, treeLabel))
            print(f"Completed processing for '{treeLabel}'.")
        outputTreePath = os.path.join(outputDir, "All_Phylogenetic_Trees.nex")
        saveTreesNexus(trees, outputTreePath)
        print("\nAll phylogenetic tree generation and saving completed successfully.")

if __name__ == "__main__":
    generalFolderPath = '/users/harry/desktop/Computational_Genetics/Final_project/code/'
    minEvolution(generalFolderPath)