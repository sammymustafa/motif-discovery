from setup import *

def compute_overlap(estimations, actual):
    common_area = 0
    for estimate in estimations:
        for truth in actual:
            start_common = max(estimate[0], truth[0])
            end_common = min(estimate[1], truth[1])
            common_area += max(0, end_common - start_common)
    return common_area


estimated_regions = [[300, 600], [2100, 3000]]  # Example output from the Viterbi algorithm
true_regions = [[400, 900]]  # Example actual regions

common_overlap = compute_overlap(estimated_regions, true_regions)

TruePositive = sum(1 for estimation in estimated_regions if compute_overlap([estimation], true_regions) >= (estimation[1] - estimation[0]) * 0.5)
FalsePositive = len(estimated_regions) - TruePositive
FalseNegative = sum(1 for truth in true_regions if compute_overlap([truth], estimated_regions) < (truth[1] - truth[0]) * 0.5)

FalsePositiveRate = FalsePositive / (FalsePositive + TruePositive)
FalseNegativeRate = FalseNegative / (FalseNegative + TruePositive)

print(FalsePositiveRate, FalseNegativeRate)
