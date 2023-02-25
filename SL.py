from matplotlib import pyplot as plot
from matplotlib.colors import LinearSegmentedColormap
import math
import numpy as np

GENES = 10

microarray_cmap = LinearSegmentedColormap('microarray', {
    'green': [(0.0, 1.0, 1.0), (0.5, 0.2, 0.2), (1.0, 0.0, 0.0)],
    'red': [(0.0, 0.0, 0.0), (0.5, 0.2, 0.2), (1.0, 1.0, 1.0)],
    'blue': [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0)],
})


def main():
    #  simulate data
    np.random.seed(1)
    simData = {"patient_%d_%s_replicate_%d" % (i//6+1, ["before", "after"][i//3%2], i%3+1): getRandoms(1, 1) for i in range(24)}
    corrMat = calcCorrMat(simData)
    newick, length = WPGMA(corrMat)
    print("Tree in Newick: \n%s;\n\nLength: %f" % (newick, length))
    fig, plots = plot.subplots(1, 1)
    fig.set_size_inches(20,20)
    treeNewick(newick + ";", plots)


def getRandoms(stdDev: float, mean: float, num=GENES) -> 'list of floats':
    return np.random.normal(mean, stdDev, num)


def calcCorrMat(simData: list) -> "dict[dict[list]]":
    #  calc correlation matrix, values are transformed to 1-correlation -> interval [0,2] with 0 as good correl(1) and 2 as bad(-1) and 
    corrMat = dict.fromkeys(simData.keys(), 0)
    for i in corrMat:
        corrMat[i] = dict.fromkeys(simData.keys(), 0)
    for i in corrMat:
        for j in corrMat[i]:
            # if j != i:
            corrMat[i][j] = 1-sum(simData[i] * simData[j] / (sum(simData[i]**2)**.5 * sum(simData[j]**2)**.5))
    return corrMat


def WPGMA(corrMat: dict) -> "(Newick string, length tree)":
    # WPGMA to cluster 
    # corrMat={"a":{"a":0, "b":17, "c":21, "d":31, "e":23},
    #          "b":{"a":17, "b":0, "c":30, "d":34, "e":21},
    #          "c":{"a":21, "b":30, "c":0, "d":28, "e":39},
    #          "d":{"a":31, "b":34, "c":28, "d":0, "e":43},
    #          "e":{"a":23, "b":21, "c":39, "d":43, "e":0}
    #          }
    while len(corrMat) > 1:
        maxCorr = (float("+inf"), None, None)

        #  find best correlation
        for i in corrMat:
            for j in corrMat[i]:
                if i != j and maxCorr[0] > corrMat[i][j]:
                    maxCorr = (corrMat[i][j], i, j)

        #  newick is the new node name
        newick = "(%s:%f,%s:%f)" % (maxCorr[1], np.abs(maxCorr[0]/2-corrMat[maxCorr[1]][maxCorr[1]]), maxCorr[2], np.abs(maxCorr[0]/2-corrMat[maxCorr[2]][maxCorr[2]]))

        #  merge best pair of tissues by arithmetic mean
        for i in corrMat[maxCorr[1]]:  # calc new distances to nodes from new node
            corrMat[maxCorr[1]][i] = (corrMat[maxCorr[1]][i] + corrMat[maxCorr[2]][i]) / 2

        # replace old node name with new node name
        corrMat[newick] = corrMat[maxCorr[1]]
        corrMat[newick][newick] = maxCorr[0]/2
        corrMat.pop(maxCorr[1])
        corrMat.pop(maxCorr[2])

        # replace old node pair with new joined node in all other nodes
        for i in corrMat:
            corrMat[i][newick] = corrMat[newick][i]
            corrMat[i].pop(maxCorr[1])
            corrMat[i].pop(maxCorr[2])

    newick, length = list(corrMat.keys())[0], list(list(corrMat.values())[0].values())[0]
    return newick, length


def plotHeatmap(plot) -> None:
    #  plot simData as heatmap
    plot.clf()
    plot.figure(figsize=(20, 10))
    plot.xticks(ticks=np.arange(len(simData)), labels=simData.keys(), rotation=90)
    plot.yticks(ticks=np.arange(GENES), labels=["gene%d" % (i+1) for i in range(GENES)])
    vals = np.array(list(simData.values()))
    rotatedVals = [[vals[j][i]for j in range(len(vals))]for i in range(GENES)]
    cmap = Cmap("rg")
    cmap.set_extremes(under="green", over="red")
    show = plot.imshow(rotatedVals, cmap=microarray_cmap, interpolation="nearest")
    plot.colorbar(show, location="left")
    plot.savefig("plot.svg")
    plot.clf()


def treeNewick(newick: str, plot):

    def plotIt(newick: str, plot, x=0, **kwargs) -> "float[tree depth]":
        """
        kwargs: treeDepth=0, lowestLeaf=x, leafLabels=True, nodeLabels=True, treeColor="black"
        """
        if newick[0] == "(" and newick[-1] == ")":
            newick = newick[1:-1]

        lowestLeaf = kwargs.get("lowestLeaf", 15)
        treeColor = "black"
        nodeLabels = kwargs.get("nodeLabels", True)
        leafLabels = kwargs.get("leafLabels", True)
        kwargs["treeDepth"] = kwargs.get("treeDepth", 0)
        myPosY = None

        if newick[0] == "(":
            splitIndex, openBrackets = 1, 1
            while openBrackets > 0:
                openBrackets += newick[splitIndex] == "("
                openBrackets -= newick[splitIndex] == ")"
                splitIndex += 1
            splitIndex += newick[splitIndex:].find(",")

            if newick[splitIndex] == ")":  # if it's sth. like (subTree)C:1
                label, length = newick[newick.rfind(")")+1:].split(":")
                kwargs["treeDepth"], kwargs["lowestLeaf"], myPosY = plotIt(newick[1: newick.rfind(")")], plot, x + float(length), **kwargs)
                if nodeLabels:
                    plot.text(x + float(length), myPosY, " " + label)
                plot.plot([x, x + float(length)], [myPosY, myPosY], color=treeColor)

            else:  # if it's sth. like (subTree)R:1,(subTree)Q:1
                upper, lower = newick[:splitIndex], newick[splitIndex+1:]
                kwargs["treeDepth"], kwargs["lowestLeaf"], upperPos = plotIt(upper, plot, x, **kwargs)
                kwargs["treeDepth"], kwargs["lowestLeaf"], lowerPos = plotIt(lower, plot, x, **kwargs)
                plot.plot([x, x], [upperPos, lowerPos], color=treeColor)
                myPosY = (upperPos - lowerPos) / 2 + lowerPos

        elif "," in newick:  # if it's sth. like A:1,(subtree)B:1
            upper, lower = newick.split(",", 1)
            kwargs["treeDepth"], kwargs["lowestLeaf"], upperPos = plotIt(upper, plot, x, **kwargs)
            kwargs["treeDepth"], kwargs["lowestLeaf"], lowerPos = plotIt(lower, plot, x, **kwargs)
            plot.plot([x, x], [upperPos, lowerPos], color=treeColor)
            myPosY = (upperPos - lowerPos) / 2 + lowerPos

        elif newick:  # if it's sth. like A:1
            label, length = newick.split(":")
            myPosY = lowestLeaf-1
            kwargs["treeDepth"] = max(kwargs.get("treeDepth", 0), x + float(length))
            if leafLabels:
                plot.text(x + float(length), myPosY, " " + label)
            plot.plot([x, x + float(length)], [myPosY, myPosY], color=treeColor)
            kwargs["lowestLeaf"] = myPosY

        return kwargs["treeDepth"], kwargs["lowestLeaf"], myPosY

    newick = newick.replace(" ", "")
    childs = newick.count(":")
    if newick[-1] != ";":
        raise ValueError("Newick string doesn't end with a semicolon")

    depth, _, _ = plotIt(newick[:-1], plot)
    plot.set_xlim([0, depth])
    plot.set_xticks([0, depth/2,  depth])
    plot.spines['bottom'].set_visible(True)
    plot.spines['top'].set_visible(False)
    plot.spines['left'].set_visible(False)
    plot.spines['right'].set_visible(False)
    plot.yaxis.set_visible(False)
    globals()["plot"].savefig("tree.png")


if __name__ == '__main__':
    main()
    # a = np.random.normal(0, 1, 10)
    # print("len=%d, avg=%f, <100: %s" % (len(a), sum(a)/len(a), np.all(abs(a) < 100)))
    # print(abs(a))