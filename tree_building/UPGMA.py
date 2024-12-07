from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from Bio import Phylo
from Bio.Phylo.Newick import Tree, Clade
import glob

def saveToPhylip(distances, labels, filePath):
    n = len(labels)
    with open(filePath, 'w') as f:
        f.write(f"{n}\n")
        for i, label in enumerate(labels):
            formattedLabel = label[:10].ljust(10)
            f.write(f"{formattedLabel} ")
            f.write(" ".join(f"{distances[i, j]:.5f}" for j in range(n)))
            f.write("\n")

def loadDistanceMatrix(filePath):
    df = pd.read_csv(filePath, index_col=0)
    df.replace("-", np.nan, inplace=True)
    distanceMatrix = df.to_numpy(dtype=float)
    iUpper = np.triu_indices_from(distanceMatrix, 1)
    distanceMatrix[(iUpper[1], iUpper[0])] = distanceMatrix[iUpper]
    return distanceMatrix, df.columns.tolist()

def linkageToPhylo(linkageMatrix, labels):
    treeRoot, nodelist = to_tree(linkageMatrix, rd=True)
    def buildClade(node):
        if node.is_leaf():
            return Clade(name=labels[node.id])
        else:
            clade = Clade()
            clade.clades.append(buildClade(node.left))
            clade.clades.append(buildClade(node.right))
            clade.branch_length = node.dist
            return clade
    phyloTree = Tree(root=buildClade(treeRoot))
    return phyloTree

def upgmaTree(distanceMatrix, labels, title, ax):
    condensedDist = distanceMatrix[np.triu_indices_from(distanceMatrix, 1)]
    linkageMatrix = linkage(condensedDist, method='average')
    dendrogram(linkageMatrix, labels=labels, orientation='left', leaf_rotation=0, leaf_font_size=10, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Distance")
    ax.set_ylabel("Clusters")
    tree = linkageToPhylo(linkageMatrix, labels)
    return tree

def UPGMA(generalFolder):
    inputFolder = os.path.join(generalFolder, "input")
    resultFolder = os.path.join(generalFolder, "result", "UMPGA")
    os.makedirs(resultFolder, exist_ok=True)
    inputFilePaths = sorted(glob.glob(os.path.join(inputFolder, "*.csv")))
    trees = []
    allPhyloTrees = []
    numFiles = len(inputFilePaths)
    cols = 2
    rows = (numFiles + 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(16, 6 * rows))
    if rows == 1 and cols == 1:
        axes = np.array([[axes]])
    elif rows == 1 or cols == 1:
        axes = axes.reshape(rows, cols)
    for idx, filePath in enumerate(inputFilePaths):
        matrix, labels = loadDistanceMatrix(filePath)
        baseName = os.path.basename(filePath)
        title = f"{os.path.splitext(baseName)[0].replace('_', ' ').title()} UPGMA Tree"
        row = idx // cols
        col = idx % cols
        tree = upgmaTree(distanceMatrix=matrix, labels=labels, title=title, ax=axes[row, col])
        trees.append(tree)
        allPhyloTrees.append(tree)
        phylipFilename = f"{os.path.splitext(baseName)[0]}_Distance_Matrix.phy"
        phylipPath = os.path.join(resultFolder, phylipFilename)
        saveToPhylip(matrix, labels, phylipPath)
    totalPlots = rows * cols
    if numFiles < totalPlots:
        for emptyIdx in range(numFiles, totalPlots):
            row = emptyIdx // cols
            col = emptyIdx % cols
            axes[row, col].axis('off')
    nexusPath = os.path.join(resultFolder, "UPGMA_Trees.nex")
    Phylo.write(allPhyloTrees, nexusPath, "nexus")
    plotPath = os.path.join(resultFolder, "UPGMA_Trees.png")
    plt.tight_layout()
    plt.savefig(plotPath, dpi=300)

if __name__ == "__main__":
    generalFolderInput = '/users/harry/desktop/Computational_Genetics/Final_project/code/'
    UPGMA(generalFolderInput)
