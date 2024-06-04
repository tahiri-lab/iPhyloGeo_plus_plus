import csv
import os
import shutil

import UserConfig
import pandas as pd
import toyplot
import toytree
from Bio.Phylo.Consensus import *
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix

from aPhyloGeo.Alignement import AlignSequences
from aPhyloGeo.MultiProcessor import Multi

userData = UserConfig.DataConfig()  # object used to store parameters provided by user

bootstrapThreshold = 0
lsThreshold = 60
windowSize = 200
stepSize = 100
# dataNames = ['ALLSKY_SFC_SW_DWN_newick', 'T2M_newick', 'QV2M_newick', 'PRECTOTCORR_newick', 'WS10M_newick']
referenceGeneFile = '../datasets/small_seq.fasta'
# fileName = '../datasets/geo.csv'
specimen = 'id'
# names = ['id', 'ALLSKY_SFC_SW_DWN', 'T2M', 'QV2M', 'PRECEPTOR', 'WS10M']
bootstrapList = []
data = []
bootstrapAmount = 100
referenceGeneDir = '../datasets/'
makeDebugFiles = True


def openCsv(file):
    """
    Open and read the csv file to get the datas

    Args:
        fileName (the file with the content to read from)
    
    Return:
        The content of the file
    """
    df = pd.read_csv(file)
    return df


def getDissimilaritiesMatrix(df, columnWithSpecimenName, columnToSearch):
    """
    Creation of a list containing the names of specimens and minimum
    temperatures

    Args:
        df (content of CSV file)
        columnWithSpecimenName (first column of names)
        columnToSearch (column to compare with the first one)
    
    Return:
        The dissimilarity matrix

    """
    meteoData = df[columnToSearch].tolist()
    nomVar = df[columnWithSpecimenName].tolist()
    nbrSeq = len(nomVar)
    maxValue = 0
    minValue = 0

    # First loop that allows us to calculate a matrix for each sequence
    tempTab = []

    for e in range(nbrSeq):
        # A list that will contain every distance before normalization
        tempList = []
        for i in range(nbrSeq):
            maximum = max(float(meteoData[e]), float(meteoData[i]))
            minimum = min(float(meteoData[e]), float(meteoData[i]))
            distance = maximum - minimum
            tempList.append(float("{:.6f}".format(distance)))

        # Allow finding the maximum and minimum value for the weather value and
        # then to add the temporary list in an array   
        if maxValue < max(tempList):
            maxValue = max(tempList)
        if minValue > min(tempList):
            minValue = min(tempList)
        tempTab.append(tempList)

    # Calculate normalised matrix
    tabDf = pd.DataFrame(tempTab)
    dmDf = (tabDf - minValue) / (maxValue - minValue)
    dmDf = dmDf.round(6)

    matrix = [dmDf.iloc[i, :i + 1].tolist() for i in range(len(dmDf))]
    dm = _DistanceMatrix(nomVar, matrix)
    return dm


def leastSquare(tree1, tree2):
    """
    Method that calculates the least square distance between two trees.
    Trees must have the same number of leaves.
    Leaves must all have a twin in each tree.
    A tree must not have duplicate leaves
     x x
    ╓╫╖ ╓╫╖
    123 312
 
    Args:
        tree1 (distanceTree object from biopython)
        tree2 (distanceTree object from biopython)
    
    Return:
        return result (double) the final distance between the two 
        
    """
    ls = 0.00
    leaves1 = tree1.get_terminals()

    leavesName = list(map(lambda l: l.name, leaves1))

    for i in leavesName:
        leavesName.pop(0)
        for j in leavesName:
            d1 = (tree1.distance(tree1.find_any(i), tree1.find_any(j)))
            d2 = (tree2.distance(tree2.find_any(i), tree2.find_any(j)))
            ls += (abs(d1 - d2))
    return ls


def drawTreesMake(trees, user_provided_names):
    """
    Function that will draw the trees for each climatic variable.
    The DistanceTreeConstructor object is transformed to Newick format and 
    loaded as a toytree MulTitre object.
    Some styling is applied and the
    resulting trees are drawn into a .pdf in the viz/ dir.
    
    Args:
        trees (Dictionary of DistanceTreeConstructor object with climatic
        variable for keys)

    """
    treesNewick = {}
    toyTrees = []
    names = user_provided_names

    for k, v in trees.items():
        treesNewick[k] = v.format('newick')
        ttree = toytree.tree(treesNewick[k], tree_format=1)
        toyTrees.append(ttree)
    mtree = toytree.mtree(toyTrees)

    # Setting up the styling for nodes
    for tree in mtree.treelist:
        tree.style.edge_align_style = {'stroke': 'black', 'stroke-width': 1}
        for node in tree.treenode.traverse():
            if node.is_leaf():
                node.add_feature('color', toytree.colors[7])
            else:
                node.add_feature('color', toytree.colors[1])
    colors = tree.get_node_values('color', show_root=1, show_tips=1)

    # Draw the climatic trees
    canvas, axes, mark = mtree.draw(nrows=round(len(mtree) / 5),
                                    ncols=len(mtree), height=400, width=1000,
                                    node_sizes=8, node_colors=colors,
                                    tip_labels_align=True);

    for i in range(len(mtree)):
        randColor = "#%03x" % random.randint(0, 0xFFF)
        axes[i].text(0, mtree.ntips, names[i + 1], style={'fill': randColor,
                                                          'font-size': '10px', 'font-weight': 'bold'});
    toyplot.png.render(canvas, '/tmp/climatic_trees.png')


def createTree(dm):
    '''
    Create a dna tree from content coming from a fasta file.

    Args:
        dm (content used to create the tree)

    Return:
        tree (the new tree)
    '''
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    return tree


def climaticPipeline(user_provided_file, user_provided_names):
    '''
    Creates a dictionary with the climatic Trees

    Return:
        trees (the climatic tree dictionary)
    '''
    trees = {}
    df = openCsv(user_provided_file)
    names = user_provided_names
    names[len(names) - 1] = names[len(names) - 1].strip()
    for i in range(1, len(names)):
        dm = getDissimilaritiesMatrix(df, names[0], names[i])
        trees[names[i]] = createTree(dm)
    # leastSquare(trees[names[1]],trees[names[2]])
    return trees


def createBoostrap(msaSet):
    '''
    Create a tree structure from sequences given by a dictionary.
    Args:
        msaSet (dictionary with multiple sequences alignment to transform into trees)
    Return:
        *********TO WRITE**********
    '''
    constructor = DistanceTreeConstructor(DistanceCalculator('identity'))

    # creation of an intermediary list
    # each list is a process
    list = []
    for key in msaSet.keys():
        list.append([msaSet, constructor, key])

    # multiprocessing
    print("Creating bootstrap variations with multiplayer of:", bootstrapAmount)
    result = Multi(list, bootSingle).processingSmallData()

    # reshaping the output into a readable dictionary
    consensusTree = {}
    for i in result:
        consensusTree[i[1]] = i[0]

    return consensusTree


def bootSingle(args):
    msaSet = args[0]
    constructor = args[1]
    key = args[2]
    result = bootstrap_consensus(msaSet[key], bootstrapAmount, constructor,
                                 majority_consensus)
    return [result, key]


def calculateAverageBootstrap(tree):
    '''
    Calculate if the average confidence of a tree

    Args:
        tree (The tree to get the average confidence from)
    Return : 
        averageBootstrap (the average Bootstrap (confidence))
    '''
    leaves = tree.get_nonterminals()
    treeConfidences = list(map(lambda l: l.confidence, leaves))
    treeConfidences.pop(0)
    totalConfidence = 0
    for confidences in treeConfidences:
        totalConfidence += confidences
    averageBootstrap = totalConfidence / len(treeConfidences)
    return averageBootstrap


def createGeneticList(geneticTrees):
    '''
    Create a list of Trees if the bootstrap Average is higher than
    the threshold

    Args :
        geneticTrees (a dictionary of genetic trees)
    Return : 
        geneticList (a list with the geneticTrees)
    '''
    geneticList = []
    for key in geneticTrees:
        bootstrap_average = calculateAverageBootstrap(geneticTrees[key])
        if (bootstrap_average >= bootstrapThreshold):
            bootstrapList.append(bootstrap_average)
            geneticList.append(key)
    return geneticList


def createClimaticList(climaticTrees):
    '''
    Create a list of climaticTrees

    Args :
        climaticTrees (a dictionary of climatic trees)
    Return : 
        climaticList(a list with the climaticTrees)
    '''
    climaticList = []
    for key in climaticTrees:
        climaticList.append(key)
    return climaticList


def getData(leavesName, ls, index, climaticList, geneticList):
    '''
    Get data from a csv file a various parameters to store into a list

    Args :
        leavesName (the list of the actual leaves)
        ls (the least square distance between two trees)
        climaticList (the list of climatic trees)
        geneticList: (the list of genetic trees)
    '''
    with open(userData.get_fileName(), 'r') as file:
        csvreader = csv.reader(file)
        for leave in leavesName:
            for row in csvreader:
                if (row[0] == leave):
                    return [referenceGeneFile, climaticList[index],
                            leave, geneticList[0],
                            str(bootstrapList[0]), str(round(ls, 2))]


def writeOutputFile(data):
    '''
    Write the datas from data list into a new csv file

    Args :
        data (the list containing the final data)
    '''
    header = ['Gene', 'Phylogeographic tree', 'Name of species',
              'Position in ASM', 'Bootstrap mean', 'Least-Square distance']
    with open("output.csv", "w", encoding="UTF8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for i in range(len(data)):
            writer.writerow(data[i])
        f.close


def filterResults(climaticTrees, geneticTrees):
    '''
    Create the final datas from the Climatic Tree and the Genetic Tree

    Args :
        climaticTrees (the dictionary containing every climaticTrees)
        geneticTrees (the dictionary containing every geneticTrees)
    '''
    # Create a list of the tree if the bootstrap is superior of the
    # bootstrap threshold
    geneticList = createGeneticList(geneticTrees)

    # Create a list with the climatic trees name
    climaticList = createClimaticList(climaticTrees)

    # Compare every genetic tree with every climatic tree. If the returned
    # value is inferior or equal to the (Least-Square distance) LS threshold, we keep the datas
    while (len(geneticList) > 0):
        leaves1 = geneticTrees[geneticList[0]].get_terminals()
        leavesName = list(map(lambda l: l.name, leaves1))
        i = 0
        for tree in climaticTrees.keys():
            ls = leastSquare(geneticTrees[geneticList[0]],
                             climaticTrees[climaticList[i]])
            if ls == None:
                raise Exception(f'The LS distance is not calculable' +
                                'pour {aligned_file}.')
            if ls <= lsThreshold:
                data.append(getData(leavesName, ls, i, climaticList,
                                    geneticList))
            i += 1
        geneticList.pop(0)
        bootstrapList.pop(0)
    # We write the datas into an output csv file
    writeOutputFile(data)


def geneticPipeline(climaticTrees):
    '''
    Get the genetic Trees from the initial file datas, so we
    can compare every valid tree with the climatic ones. In the 
     end, it calls a method that creates a final csv file with all
    the data that we need for the comparison

    Args:
        climaticTrees (the dictionary of climaticTrees)
    '''
    ####### JUST TO MAKE THE DEBUG FILES ####### 
    if os.path.exists("./debug"):
        shutil.rmtree("./debug")
    if makeDebugFiles:
        os.mkdir("./debug")
    ####### JUST TO MAKE THE DEBUG FILES ####### 

    alignmentObject = AlignSequences()
    # alignedSequences = alignmentObject.aligned
    # heuristicMSA = alignmentObject.heuristicMSA
    # windowedSequences = alignmentObject.windowed
    msaSet = alignmentObject.msaSet
    geneticTrees = createBoostrap(msaSet)
    filterResults(climaticTrees, geneticTrees)


def createSeqAlign(self):
    align = AlignSequences.alignSequences()
    return align


def createGenTree(align_obj):
    '''
    Create a genetic trees dictionary
    Args:
     align_obj (object from AlignSequences() class)
    Return:
     geneticTrees (the genetic trees' dictionary)
    '''
    msaSet = align_obj.msaSet
    geneticTrees = createBoostrap(msaSet)
    return geneticTrees


def createAndSaveTree(user_provided_file, user_provided_names):
    """
    Draw and show climatic trees
    Args: user_provided_file (name of the climatic data file),
        user_provided_names (names of the columns in the climatic file)
    Return:
     draw Trees make (trees, user_provided_names) (call to function that draws and returns the tree)
    """
    trees = climaticPipeline(user_provided_file, user_provided_names)
    return drawTreesMake(trees, user_provided_names)
