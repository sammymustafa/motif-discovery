from setup import *

def assignPoints(tbl, ctrs):
    """Assign each of the points in tbl to the cluster with center in ctrs"""
    ptsAsgn = []
    for i in tbl:
        dists = [np.linalg.norm(np.array(i) - np.array(center)) for center in ctrs]
        min_distance_index = np.argmin(dists)
        ptsAsgn.append(min_distance_index)
    return ptsAsgn


def recalculateCtrs(tbl, ctrs, ptsAsgn):
    """Update the centroids based on the points assigned to them"""
    newCtrs = []
    for i in range(len(ctrs)):
        points = [tbl[j] for j in range(len(tbl)) if ptsAsgn[j] == i]
        if points:
            avg = np.mean(points, axis=0)
            newCtrs.append(avg)
        else:
            newCtrs.append(ctrs[i])
    return newCtrs

def run_kmeans(dataTable):
    """initializes centroids, stop criterion and step counting for clustering"""
    newCtrs = [[5, 0], [5, 40], [5, 80]]
    ptMemb = assignPoints(dataTable, newCtrs)
    stopCrit = False
    stepCount = 0

    """performs k-means clustering, plotting the clusters at each step"""
    while not stopCrit:
        stepCount += 1
        plotClusters(dataTable, ptMemb, newCtrs, stepCount)
        ptMemb = assignPoints(dataTable, newCtrs)
        oldCtrs = newCtrs
        newCtrs = recalculateCtrs(dataTable, newCtrs, ptMemb)
        """stop criterion - when centroids' total movement after a step is below
            the threshold, stop the algorithm"""
        stopDist = 0;
        for i in range(len(newCtrs)):
            stopDist += np.linalg.norm(np.array(oldCtrs[i]) - np.array(newCtrs[i]))
        if stopDist < 5:
            stopCrit = True

def plotClusters(tbl, ptMemb, cntrs, stepCnt):
    """Generate a scatterplot of the current k-means cluster assignments"""
    pt_colors = ["salmon", "lightgreen", "lightblue"]
    ctr_colors = ["red", "green", "blue"]
    for i in range(len(cntrs)):
        pts = [tbl[j] for j in range(len(tbl)) if ptMemb[j] == i]
        plt.scatter([pt[0] for pt in pts], [pt[1] for pt in pts], color = pt_colors[i], s = 10)
        plt.scatter([cntrs[i][0]], [cntrs[i][1]], color = ctr_colors[i], s = 100, facecolors = "none")
    plt.title("Step " + str(stepCnt))
    plt.xlabel("BRCA1 Gene Expression")
    plt.ylabel("SOX2 Gene Expression")
    plt.show()