from setup import *
import pyranges

def scale_matrix(mat):
    sums_rows = mat[1:, 1:].sum(axis=1)
    scaled_mat = mat.copy()
    scaled_mat[1:, 1:] = mat[1:, 1:] / sums_rows[:, np.newaxis]
    return scaled_mat

letters = ["", "A", "C", "G", "T"]
categories = ["", "notCpG", "CpG"]

emit_matrix = np.full((3, 5), 1e-3, dtype=object)
trans_matrix = np.full((3, 3), 1e-3, dtype=object)

emit_matrix[0, :] = letters
emit_matrix[:, 0] = categories
trans_matrix[0, :] = categories
trans_matrix[:, 0] = categories

sequence_file = Fasta("chr21.fa.gz")
cpg_data = pd.read_csv("cpgIslandExt.txt.gz", sep="\t", header=None)
cpg_data = cpg_data.rename(columns={1: "Chromosome", 2: "Start", 3: "End"}).iloc[:, [1,2,3]]
cpg_ranges = pyranges.PyRanges(cpg_data)

for sequence in sequence_file:
    annotations = np.zeros(len(sequence), dtype=int)
    loc_rng = pyranges.from_dict({"Chromosome": [sequence.name], "Start": [0], "End": [len(sequence)]})
    sequence = ''.join(random.choice('ACGT') if base == 'N' else base for base in str(sequence))
    overlap = loc_rng.intersect(cpg_ranges).df

    for idx in range(overlap.shape[0]):
        start_idx, end_idx = overlap["Start"].values[idx], overlap["End"].values[idx]
        annotations[start_idx:end_idx] = 1

    current_label = annotations[0]
    for idx in range(len(sequence) - 1):
        for letter in letters[1:]:
            if sequence[idx] == letter:
                emit_matrix[current_label + 1, letters.index(letter)] += 1

        next_label = annotations[idx + 1]
        trans_matrix[current_label + 1, next_label + 1] += 1
        current_label = next_label

emit_matrix = scale_matrix(emit_matrix)
trans_matrix = scale_matrix(trans_matrix)

print("3. Final Transition Probability Matrix:")
print(trans_matrix)

print("4. Final Emission Probability Matrix:")
print(emit_matrix)