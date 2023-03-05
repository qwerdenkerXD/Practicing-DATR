from matplotlib import pyplot as plot
from matplotlib.gridspec import GridSpec
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
    simDataPatient = {"patient_%d_%s_replicate_%d" % (i//6+1, ["before", "after"][i//3%2], i%3+1): getRandoms() for i in range(24)}

    simDataGenes = {key: value[0] for key, value in zip(["gene%d" % i for i in range(1, GENES+1)], zip(simDataPatient.values()))}
    corrMatPatient = calcCorrMat(simDataPatient)
    corrMatGenes = calcCorrMat(simDataGenes)
    specs = GridSpec(2, 2, width_ratios=[1, .1], height_ratios=[.2, .5], wspace=0.0, hspace=0, top=0.95, bottom=0.34, left=-.07, right=.97)

    fig = plot.figure(figsize=(20, 10))
    plots = [plot.subplot(specs[0, 0]), plot.subplot(specs[1, 0]), plot.subplot(specs[1, 1])]
    topDendo, heat, rightDendo = plots

    newickPat, depth = WPGMA(corrMatPatient)
    labels = treeNewick(newickPat + ";", topDendo, root_loc="top", leafLabels=True, showScale=False)
    keys, values = simDataPatient.keys(), simDataPatient.values()
    keys, values = zip(*sorted(zip(keys, values), key=lambda k: labels.index(k[0])))

    newickGen, depth = WPGMA(corrMatGenes)
    labels = treeNewick(newickGen + ";", rightDendo, root_loc="right", leafLabels=False, showScale=False)
    lim1, lim2 = topDendo.get_xlim()
    topDendo.set_xlim([lim1-7.4, lim2-.68])
    genes, values = zip(*sorted(zip(["gene%d" % i for i in range(1, GENES + 1)], zip(*values)), key=lambda k: labels.index(k[0])))
    values = zip(*values)
    simDataPatient = {key: value for key, value in zip(keys, values)}

    plotHeatmap(heat, simDataPatient, genes) # gene labels false
    plot.savefig("heatmap.png")
    specs = GridSpec(1, 2, width_ratios=[1, 1], height_ratios=[1], wspace=.1, hspace=0, top=0.95, bottom=0.34, left=0.1, right=.97)
    fig = plot.figure(figsize=(20, 8))
    plots = [plot.subplot(specs[0, 0]), plot.subplot(specs[0, 1])]
    for newick, subP in zip([newickPat, newickGen], plots):
        treeNewick(newick + ";", subP, root_loc="top")
    plot.savefig("trees.png")


def getRandoms(num=GENES) -> 'list of floats':
    res = np.array([])
    while len(res) < num:
        random = np.random.normal(0, 100)
        if abs(random) <= 100:
            res = np.append(res, random)
    return res


def calcCorrMat(data: dict) -> "dict[dict[list]]":
    #  calc correlation matrix, values are transformed to 1-correlation -> interval [0,2] with 0 as good correl(1) and 2 as bad(-1) and 
    corrMat = dict.fromkeys(data.keys(), 0)
    for i in corrMat:
        corrMat[i] = dict.fromkeys(data.keys(), 0)
    for i in corrMat:
        for j in corrMat[i]:
            corrMat[i][j] = 1-sum(data[i] * data[j] / (sum(data[i]**2)**.5 * sum(data[j]**2)**.5))
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


def plotHeatmap(plot, data, dataLabels=[]) -> None:
    plot.set_xticks(ticks=np.arange(len(data)), labels=data.keys(), rotation=90, fontsize=15)
    plot.set_yticks(ticks=np.arange(GENES), labels=dataLabels, fontsize=15)
    vals = np.array(list(data.values()))
    rotatedVals = [[vals[j][i]for j in range(len(vals))]for i in range(GENES)]
    show = plot.imshow(rotatedVals, cmap=microarray_cmap, interpolation="nearest", aspect="auto")
    bar = globals()["plot"].colorbar(show, location="left", ax=plot)
    bar.ax.tick_params(labelsize=15)


def treeNewick(newick: str, plot, root_loc="left", leafLabels=True, showScale=True):

    def plotIt(newick: str, plot, x=0, **kwargs) -> "float[tree depth]":
        """
        kwargs: treeDepth=0, lowestLeaf=x, leafLabels=True, nodeLabels=True, treeColor="black"
        """
        if newick[0] == "(" and newick[-1] == ")":
            newick = newick[1:-1]

        lowestLeaf = kwargs.get("lowestLeaf", 15)
        treeColor = kwargs.get("treeColor", "black") 
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
                #if nodeLabels:
                #    plot.text(x + float(length), myPosY, " " + label)
                draw([x, x + float(length)], [myPosY, myPosY], treeColor=treeColor)

            else:  # if it's sth. like (subTree)R:1,(subTree)Q:1
                upper, lower = newick[:splitIndex], newick[splitIndex+1:]
                kwargs["treeDepth"], kwargs["lowestLeaf"], upperPos = plotIt(upper, plot, x, **kwargs)
                kwargs["treeDepth"], kwargs["lowestLeaf"], lowerPos = plotIt(lower, plot, x, **kwargs)
                draw([x, x], [upperPos, lowerPos], treeColor=treeColor)
                myPosY = (upperPos - lowerPos) / 2 + lowerPos

        elif "," in newick:  # if it's sth. like A:1,(subtree)B:1
            upper, lower = newick.split(",", 1)
            kwargs["treeDepth"], kwargs["lowestLeaf"], upperPos = plotIt(upper, plot, x, **kwargs)
            kwargs["treeDepth"], kwargs["lowestLeaf"], lowerPos = plotIt(lower, plot, x, **kwargs)
            draw([x, x], [upperPos, lowerPos], treeColor=treeColor)
            myPosY = (upperPos - lowerPos) / 2 + lowerPos

        elif newick:  # if it's sth. like A:1
            label, length = newick.split(":")
            myPosY = lowestLeaf-1
            kwargs["treeDepth"] = max(kwargs.get("treeDepth", 0), x + float(length))
            draw(x + float(length), myPosY, label)  # plot label
            draw([x, x + float(length)], [myPosY, myPosY], treeColor=treeColor)
            kwargs["lowestLeaf"] = myPosY

        return kwargs["treeDepth"], kwargs["lowestLeaf"], myPosY

    def draw(x, y, label=None, treeColor="black"):
        if label is None:
            if root_loc in ("left", "right"):
                plot.plot(x, y, color=treeColor)
            else:
                plot.plot(y, x, color=treeColor)
        else:
            label = " %s " % label
            alpha = float(leafLabels)
            if root_loc == "left":
                plot.text(x, y, label, alpha=alpha, verticalalignment="center")
            elif root_loc == "right":
                plot.text(x, y, label, alpha=alpha, verticalalignment="center", horizontalalignment="right")
            elif root_loc == "bottom":
                plot.text(y, x, label, alpha=alpha, horizontalalignment="center", rotation=90)
            else:
                plot.text(y, x, label, alpha=alpha, verticalalignment="top", horizontalalignment="center", rotation=90)

    newick = newick.replace(" ", "")
    childs = newick.count(":")
    if newick[-1] != ";":
        raise ValueError("Newick string doesn't end with a semicolon")

    depth, _, _ = plotIt(newick[:-1], plot, leafLabels=leafLabels)
    plot.spines['bottom'].set_visible(False)
    plot.spines['top'].set_visible(False)
    plot.spines['left'].set_visible(False)
    plot.spines['right'].set_visible(False)
    plot.set_xlim([0, depth])
    plot.set_xticks([0, depth/2,  depth])
    plot.yaxis.set_visible(False)
    plot.xaxis.set_visible(False)
    if showScale:
        plot.xaxis.set_visible(True)
        plot.spines['bottom'].set_visible(True)
    allLabels = [text.get_text()[1:-1] for text in plot.texts]
    if root_loc == "right":
        plot.set_xticklabels([depth, depth/2, 0])
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
        plot.set_yticks([0, depth/2,  depth])
        if root_loc == "top":
            plot.set_yticklabels([depth, depth/2, 0])
            plot.invert_yaxis()

    return allLabels


if __name__ == '__main__':
    main()
    # a = np.random.normal(0, 1, 10)
    # print("len=%d, avg=%f, <100: %s" % (len(a), sum(a)/len(a), np.all(abs(a) < 100)))
    # print(abs(a))