from setup import *
from kmeans_clust import run_kmeans

!wget -c https://www.dropbox.com/sh/u73ktzeykydjsos/AABtGwiaaE-iMOd-FWNXJE6Da?dl=0 -O tissue_data.zip
!unzip -o tissue_data.zip -d tissue_data
os.unlink("tissue_data.zip")
tissue_data = {}
for fname in os.listdir("tissue_data"):
    tissue_data[fname] = pd.read_csv(os.path.join("tissue_data", fname), sep="\t", header=None).values.tolist()

# Evaluation on Tissue 1
run_kmeans(tissue_data["tissue1_data.txt"])

# Evaluation on Tissue 2
run_kmeans(tissue_data["tissue2_data.txt"])

# Other Method
def run_kmeans(data_table):
    """initializes centroids, stop criterion and step counting for clustering"""
    newCtrs = [[5, 0], [5, 40], [5, 80]]
    ptMemb = assignPoints(data_table, newCtrs)
    stopCrit = False
    stepCount = 0

    """performs k-means clustering, plotting the clusters at each step"""
    while not stopCrit:
        stepCount += 1

        plotClusters(data_table, ptMemb, newCtrs, stepCount)

        ### YOUR CODE HERE ###
        ptMemb = assignPoints(data_table, newCtrs)
        if len(newCtrs) >= 3:
            unique_clusters = list(set(ptMemb))
            cluster_sizes = [sum(1 for p in ptMemb if p == cluster) for cluster in unique_clusters]
            max_size, min_size = max(cluster_sizes, default = 0), min(cluster_sizes, default = 0)

            if max_size > 1.1 * min_size:
                newCtrs.append([np.random.uniform(), np.random.uniform()])

        oldCtrs = newCtrs
        updated_centroids = recalculateCtrs(data_table, newCtrs, ptMemb)

        movement = sum(np.linalg.norm(np.array(prev) - np.array(curr)) for prev, curr in zip(oldCtrs, updated_centroids))

        if movement < 5:
            stopCrit = True
        else:
            newCtrs = updated_centroids

    return newCtrs

# Evaluation of Other Method
run_kmeans(tissue_data["tissue2_data.txt"])
