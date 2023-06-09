"""
A module to practice clustering.

Instead of using modules for WPGMA and dendogram plot, I implemented them by myself for more practice and, well, Python is fun :)

GitHub: https://github.com/qwerdenkerXD/Practicing-DATR
"""

from matplotlib import pyplot as plot
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import math
import numpy as np
from numpy.typing import NDArray

GENES = 10
PATIENTS = 4
REPLICATES = 3

microarray_cmap = LinearSegmentedColormap('microarray', {
    'green': [(0.0, 1.0, 1.0), (0.5, 0.2, 0.2), (1.0, 0.0, 0.0)],
    'red': [(0.0, 0.0, 0.0), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0)],
    'blue': [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
})


def main():
    """
    The main-method generates fictive measurements of patient tissues concerning breast cancer
    before and after a treatment on different genes, then clusters it with WPGMA-algorithm and
    plots it to a heatmap of the correlations with the dendograms.
    """

    #  simulate data
    np.random.seed(1)

    # generate random values for each gene
    simDataGenes = {}
    cancerGenes = np.random.randint(GENES, size=GENES // 3)
    for g in range(1, GENES+1):
        if g in cancerGenes:
            simDataGenes["gene%d" % g] = getRandoms(30, 70, PATIENTS * REPLICATES * 2)
        else:
            simDataGenes["gene%d" % g] = getRandoms(-70, 30, PATIENTS * REPLICATES * 2)

    # simDataGenes contains per gene values in order [patient_1_before..., patient_1_after..., ...], so simDataPatient gets per patient the values of its index
    patientData = list(zip(*simDataGenes.values()))  # transpose

    simDataPatient = {}
    for p in range(1, 5):
        for r in range(1, 4):
            before, after = np.array(patientData[(p-1)*6 + 2 * r - 2]), np.array(patientData[(p-1)*6 + 2 * r - 1])

            # check if data makes sense (if healthy before, it shouldn't be diseased after)
            for i in range(10):
                if after[i] > 0 and before[i] <= 0:
                    after[i] = -after[i]
                    simDataGenes["gene%d" % (i+1)][(p-1) * (r-1)] = after[i]

            simDataPatient["patient_%d_after_replicate_%d" % (p, r)] = after
            simDataPatient["patient_%d_before_replicate_%d" % (p, r)] = before

    # calc correlations
    corrMatPatient = calcCorrMat(simDataPatient)
    corrMatGenes = calcCorrMat(simDataGenes)

    # gridspec is necessary to adjust the distances between plots, so the dendograms can fit to the borders of the heatmap
    specs = GridSpec(2, 2, width_ratios=[1, .1], height_ratios=[.5, .5], wspace=0, hspace=0.0, top=0.95, bottom=0.34, left=-.07, right=.97)

    fig = plot.figure(figsize=(20, 10))
    plots = [plot.subplot(specs[0, 0]), plot.subplot(specs[1, 0]), plot.subplot(specs[1, 1])]
    topDendo, heat, rightDendo = plots

    # clustering and plotting the patients' dendogram to the top -> topDendo
    newickPat, depth = WPGMA(corrMatPatient)
    labels = treeNewick(newickPat + ";", topDendo, root_loc="top", leafLabels=True, showScale=False)

    # now get the keys and values from simDataPatient sorted by the labels to fit to the dendogram
    keys, values = simDataPatient.keys(), simDataPatient.values()
    keys, values = zip(*sorted(zip(keys, values), key=lambda k: labels.index(k[0])))

    # clustering and plotting the genes' dendogram to the right -> rightDendo
    newickGen, depth = WPGMA(corrMatGenes)
    labels = treeNewick(newickGen + ";", rightDendo, root_loc="right", leafLabels=False, showScale=False)

    # stretch the top dendogram to fit to the heatmap's columns
    lim1, lim2 = topDendo.get_xlim()
    topDendo.set_xlim([lim1-7.4, lim2-.68])

    # now sort the values by the genes' labels to fit to the dendogram
    genes, values = zip(*sorted(zip(["gene%d" % i for i in range(1, GENES + 1)], zip(*values)), key=lambda k: labels.index(k[0])))
    values = zip(*values)

    # update simDataPatient with new sorted entries
    simDataPatient = {key: value for key, value in zip(keys, values)}

    plotHeatmap(heat, simDataPatient, genes)
    plot.savefig("heatmap.png")

    # now plot the dendograms itself to a single extra file
    specs = GridSpec(1, 2, width_ratios=[1, 1], height_ratios=[1], wspace=.1, hspace=0, top=0.95, bottom=0.34, left=0.1, right=.97)
    fig = plot.figure(figsize=(20, 8))
    plots = [plot.subplot(specs[0, 0]), plot.subplot(specs[0, 1])]
    for newick, subP in zip([newickPat, newickGen], plots):
        treeNewick(newick + ";", subP, root_loc="top")
    plot.savefig("trees.png")


def getRandoms(mean, stddev, num) -> 'NDArray[np.float64]':
    """
    Generates num random values between +-100
    """
    res = np.array([])
    while len(res) < num:
        random = np.random.normal(mean, stddev)
        if abs(random) <= 100:
            res = np.append(res, random)
    return res


def calcCorrMat(data: dict) -> "dict[str, dict[str, int]]":
    """
    Calculates a correlation matrix from data generated in main()
    Values are transformed to 1-correlation -> interval [0,2] with 0 as high correlation(1) and 2 as anti-correlation(-1)
    """
    corrMat = dict.fromkeys(data.keys(), {})
    for i in corrMat:
        corrMat[i] = dict.fromkeys(data.keys(), 0)
    for i in corrMat:
        for j in corrMat[i]:
            corrMat[i][j] = 1-sum(data[i] * data[j] / (sum(data[i]**2)**.5 * sum(data[j]**2)**.5))
    return corrMat


def WPGMA(corrMat: dict) -> "tuple[str, float]":
    """
    WPGMA to cluster, input will be changed, so make a copy

    corrMat like {"a":{"a":0, "b":17, "c":21, "d":31, "e":23},
                  "b":{"a":17, "b":0, "c":30, "d":34, "e":21},
                  "c":{"a":21, "b":30, "c":0, "d":28, "e":39},
                  "d":{"a":31, "b":34, "c":28, "d":0, "e":43},
                  "e":{"a":23, "b":21, "c":39, "d":43, "e":0}
                  }

    #Return:
    tuple of Newick-formatted string of clustered tree and its depth
    """
    while len(corrMat) > 1:
        maxCorr = (float("+inf"), None, None)  # correl and its indices in corrMat

        #  find best correlation
        for i in corrMat:
            for j in corrMat[i]:
                if i != j and maxCorr[0] > corrMat[i][j]:
                    maxCorr = (corrMat[i][j], i, j)

        correl, child_1, child_2 = maxCorr

        # in corrMat[maxCorr[i]][maxCorr[i]] is the current distance to the root
        # this can be done because in original WPGMA is this index ignored (because 0)
        child_1_distance = np.abs(correl/2-corrMat[child_1][child_1])
        child_2_distance = np.abs(correl/2-corrMat[child_2][child_2])

        #  newick is the new node name
        newick = "(%s:%f,%s:%f)" % (child_1, child_1_distance, child_2, child_2_distance)

        #  merge best pair of tissues by arithmetic mean
        for i in corrMat[child_1]:  # calc new distances to nodes from new node
            corrMat[child_1][i] = (corrMat[child_1][i] + corrMat[child_2][i]) / 2

        # replace old node name with new node name
        corrMat[newick] = corrMat[child_1]
        corrMat[newick][newick] = maxCorr[0]/2
        corrMat.pop(child_1)
        corrMat.pop(child_2)

        # replace old node pair with new joined node in all other nodes
        for i in corrMat:
            corrMat[i][newick] = corrMat[newick][i]
            corrMat[i].pop(child_1)
            corrMat[i].pop(child_2)

    newick, depth = list(corrMat.keys())[0], list(list(corrMat.values())[0].values())[0]
    return newick, depth


def plotHeatmap(plot, data, dataLabels=[]) -> None:
    """
    dataLabels: because the measurements are stored in arrays, there are labels needed for each index, so here the genes
    """
    plot.set_xticks(ticks=np.arange(len(data)), labels=data.keys(), rotation=90, fontsize=15)
    plot.set_yticks(ticks=np.arange(GENES), labels=dataLabels, fontsize=15)
    vals = np.array(list(data.values()))
    rotatedVals = [[vals[j][i]for j in range(len(vals))]for i in range(GENES)]
    show = plot.imshow(rotatedVals, cmap=microarray_cmap, interpolation="nearest", aspect="auto")
    bar = globals()["plot"].colorbar(show, location="left", ax=plot)
    bar.ax.tick_params(labelsize=15)


def treeNewick(newick: str, plot, root_loc="left", leafLabels=True, showScale=True):
    """
    Plots a dendogram from a newick-formatted tree.
    This method is currently developed for this script only and works well with it,
    but for general use, there is no check, if the newick string is correct, and no
    method to make it compatible.

    Works properly with the WPGMA-returned newick with following EBNF:
    tree => "(" (tree|label) ":" distance "," (tree|label) ":" distance ")" ";"
    label => str
    distance => float

    other possible strings aren't tested
    """

    def recursive(newick: str, plot, x=0., **kwargs) -> "tuple[float, float, float]":
        """
        Plots a dendogram recursively from a given newick string, root is left.
        This is done in 3 steps:
        1. reach child-level: There are being the child labels plotted and their horizontal branches drawn
        loop till end:
            2. connection: horizontal branches of the next level are drawn recursively, so here will they be connected by the returned y-positions of these branches
            3. root-drawing: if two branches were connected horizontally, there is the root branch following, fitting to the middle of the vertical connection branch

        x: the x-position where to plot the root
        kwargs: treeDepth=0, lowestLeaf=x, leafLabels=True, nodeLabels=True, treeColor="black"

        #Return:
        tuple of current tree depth, y-value of lowest leaf and self y-value
        """
        if newick[0] == "(" and newick[-1] == ")":
            newick = newick[1:-1]

        lowestLeaf = kwargs.get("lowestLeaf", 15)  # important to have the same distance between the children/leafs
        treeColor = kwargs.get("treeColor", "black")
        nodeLabels = kwargs.get("nodeLabels", True)  # not used for this script, but can be implemented
        leafLabels = kwargs.get("leafLabels", True)
        kwargs["treeDepth"] = kwargs.get("treeDepth", 0)
        myPosY = 0  # necessary to return it to center the previous root properly

        if newick[0] == "(":

            # detect end of the subtree by finding its closing bracket
            splitIndex, openBrackets = 1, 1  # splitIndex is one behind the character
            while openBrackets > 0:
                openBrackets += newick[splitIndex] == "("
                openBrackets -= newick[splitIndex] == ")"
                splitIndex += 1
            splitIndex += newick[splitIndex:].find(",")  # if there is a subtree following like in (subtree)A:1,(subtree)B:1 then the following comma is the end

            # splitIndex is now the index of the comma after the tree or (if no comma there) the last bracket's index

            if newick[splitIndex] == ")":  # if it's sth. like (subTree)C:1             -> root-drawing (horizontal)
                label, length = newick[splitIndex+1:].split(":")
                kwargs["treeDepth"], kwargs["lowestLeaf"], myPosY = recursive(newick[1: splitIndex], plot, x + float(length), **kwargs)
                # if nodeLabels:
                #    plot.text(x + float(length), myPosY, " " + label)
                draw([x, x + float(length)], [myPosY, myPosY], treeColor=treeColor)

            else:  # if it's sth. like (subTree)R:1,(subTree)Q:1                        -> connection (vertical)
                upper, lower = newick[:splitIndex], newick[splitIndex+1:]
                kwargs["treeDepth"], kwargs["lowestLeaf"], upperPos = recursive(upper, plot, x, **kwargs)
                kwargs["treeDepth"], kwargs["lowestLeaf"], lowerPos = recursive(lower, plot, x, **kwargs)
                draw([x, x], [upperPos, lowerPos], treeColor=treeColor)
                myPosY = (upperPos - lowerPos) / 2 + lowerPos  # centering for the previous root-line

        elif "," in newick:  # if it's sth. like A:1,(subtree)B:1                       -> connection (vertical)
            upper, lower = newick.split(",", 1)
            kwargs["treeDepth"], kwargs["lowestLeaf"], upperPos = recursive(upper, plot, x, **kwargs)
            kwargs["treeDepth"], kwargs["lowestLeaf"], lowerPos = recursive(lower, plot, x, **kwargs)
            draw([x, x], [upperPos, lowerPos], treeColor=treeColor)
            myPosY = (upperPos - lowerPos) / 2 + lowerPos

        elif newick:  # if it's sth. like A:1  # will only be reached if A:1 is a leaf  -> child-level (horizontal)
            label, length = newick.split(":")
            myPosY = lowestLeaf-1
            kwargs["treeDepth"] = max(kwargs.get("treeDepth", 0), x + float(length))
            draw(x + float(length), myPosY, label)  # plot label
            draw([x, x + float(length)], [myPosY, myPosY], treeColor=treeColor)
            kwargs["lowestLeaf"] = myPosY

        return kwargs["treeDepth"], kwargs["lowestLeaf"], myPosY

    def draw(x, y, label=None, treeColor="black"):
        """
        Simplifies plotting in right location and alignment
        """
        if label is None:  # -> draw a line
            if root_loc in ("left", "right"):
                plot.plot(x, y, color=treeColor)
            else:
                plot.plot(y, x, color=treeColor)

        else:  # -> write text into the plot (for children)
            label = " %s " % label
            alpha = float(leafLabels)  # labels are always plotted to return the order of children, alpha makes them invisible if wanted
            if root_loc == "left":
                plot.text(x, y, label, alpha=alpha, verticalalignment="center")
            elif root_loc == "right":
                plot.text(x, y, label, alpha=alpha, verticalalignment="center", horizontalalignment="right")
            elif root_loc == "bottom":
                plot.text(y, x, label, alpha=alpha, horizontalalignment="center", rotation=90)
            else:
                plot.text(y, x, label, alpha=alpha, verticalalignment="top", horizontalalignment="center", rotation=90)

    newick = newick.replace(" ", "")
    if newick[-1] != ";":
        raise ValueError("Newick string doesn't end with a semicolon")

    depth, _, _ = recursive(newick[:-1], plot, leafLabels=leafLabels)

    plot.spines['bottom'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.spines['left'].set_visible(False)
    plot.spines['right'].set_visible(False)

    plot.set_xlim([0, depth])
    plot.set_xticks([0, round(depth/2, 2),  round(depth, 2)])

    plot.yaxis.set_visible(False)
    plot.xaxis.set_visible(False)

    if showScale:
        plot.xaxis.set_visible(True)
        plot.spines['bottom'].set_visible(True)

    allLabels = [text.get_text()[1:-1] for text in plot.texts]

    # the following mirrors the plot to locate the root correctly
    if root_loc == "right":
        plot.set_xticklabels([round(depth, 2),  round(depth/2, 2), 0])
        plot.invert_xaxis()

    elif root_loc != "left":
        allLabels = allLabels[::-1]
        if showScale:
            plot.spines["left"].set_visible(True)
            plot.spines["bottom"].set_visible(False)
            plot.yaxis.set_visible(True)
            plot.xaxis.set_visible(False)
        plot.autoscale("x")
        plot.set_ylim([0, depth])
        plot.set_yticks([0, round(depth/2, 2),  round(depth, 2)])

        if root_loc == "top":
            plot.set_yticklabels([round(depth, 2),  round(depth/2, 2), 0])
            plot.invert_yaxis()

    return allLabels


if __name__ == '__main__':
    main()
