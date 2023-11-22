from setup import *

!wget -q -c https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr21.fa.gz -O chr21.fa.gz
!gunzip -f chr21.fa.gz
!samtools faidx chr21.fa
!wget -q -c https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr22.fa.gz -O chr22.fa.gz
!gunzip -f chr22.fa.gz
!samtools faidx chr22.fa
!wget -q -c https://www.encodeproject.org/files/ENCFF029THO/@@download/ENCFF029THO.bed.gz -O ENCFF029THO.bed.gz
!wget -q -c https://www.encodeproject.org/files/ENCFF237RIC/@@download/ENCFF237RIC.bed.gz -O ENCFF237RIC.bed.gz


### ENCODE narrowPeak format
def load_narrow_peak(bed, chr:str=None):
    df = pd.read_csv(bed, header=None, sep="\t")
    df.columns = ["Chromosome", "Start", "End", "name", "score", "Strand", "signalValue", "pValue", "qValue", "peak"]
    if chr is not None:
      df = df.loc[df["Chromosome"]==chr, :]
    return pyranges.PyRanges(df)

tbl = {}
for sample in ["ENCFF029THO", "ENCFF237RIC"]:
    for chr in ["chr21", "chr22"]:
      seq = pyranges.get_sequence(load_narrow_peak("%s.bed.gz" % sample, chr), "%s.fa" % chr)
      tbl["%s_%s" % (sample, chr)] = seq
      print(sample, chr, len(seq))


### To Implement as well: get background

# Get sequences from a FASTA file for given ranges
def extract_sequence(fasta_path, ranges):
    fasta_seq = Fasta(fasta_path)
    extracted_sequences = []
    for _, entry in ranges.iterrows():
        chrom = entry['Chromosome']
        begin = entry['Start']
        stop = entry['End']
        extracted_sequences.append(str(fasta_seq[chrom][begin:stop]))
    return extracted_sequences

# Retrieve narrowPeak data and the corresponding sequences
def retrieve_narrow_peak(peak_file, chr:str=None, fasta:str=None):
    df = pd.read_csv(peak_file, sep="\t", header=None)
    df.columns = ["Chromosome", "Start", "End", "name", "score", "Strand", "signalValue", "pValue", "qValue", "peak"]
    if chr:
        df = df[df["Chromosome"] == chr]
    if fasta:
        return extract_sequence(fasta, df)
    return []

# Collecting data
data_table = {}
for sample_id in ["ENCFF029THO", "ENCFF237RIC"]:
    for chr_id in ["chr21", "chr22"]:
        fasta_file_name = f"{chr_id}.fa"
        bed_file_name = f"{sample_id}.bed.gz"
        sequences_list = retrieve_narrow_peak(bed_file_name, chr_id, fasta_file_name)
        data_table[f"{sample_id}_{chr_id}"] = sequences_list
        print(sample_id, chr_id, len(sequences_list))

# Sequences are already loaded for further computations
input_sequences = data_table["ENCFF029THO_chr21"]
motif_width = 14
MatrixM = (np.random.rand(4, motif_width) + 1e-8) / np.sum(np.random.rand(4, motif_width) + 1e-10, axis=0)
MatrixB = np.random.rand(4, 1)
max_iter = 1000

# Map nucleotide to corresponding index
def nucleotide_index(nucleotide):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return mapping.get(nucleotide, -1)

# Estimating likelihood
def estimate_likelihood(subsequence, MatrixM):
    product = 1
    for idx, nucleotide in enumerate(subsequence):
        nucleotide_idx = nucleotide_index(nucleotide)
        product *= MatrixM[nucleotide_idx, idx]
    return product

# E
def expectation_step(input_sequences, MatrixM, motif_width):
    optimal_positions = []
    for sequence in input_sequences:
        optimal_pos = np.argmax([estimate_likelihood(sequence[i:i + motif_width], MatrixM) for i in range(len(sequence) - motif_width + 1)])
        optimal_positions.append(optimal_pos)
    return optimal_positions

# M
def maximization_step(input_sequences, optimal_positions, motif_width):
    UpdatedMatrixM = np.zeros((4, motif_width))
    UpdatedMatrixB = np.zeros(4)
    for idx, sequence in enumerate(input_sequences):
        for pos in range(len(sequence) - motif_width + 1):
            is_optimal = pos == optimal_positions[idx]
            for offset, nucleotide in enumerate(sequence[pos:pos + motif_width]):
                index = nucleotide_index(nucleotide)
                if is_optimal:
                    UpdatedMatrixM[index, offset] += 1
                else:
                    UpdatedMatrixB[index] += 1
    UpdatedMatrixM = (UpdatedMatrixM + 1e-10) / np.sum(UpdatedMatrixM + 1e-10, axis=0)
    UpdatedMatrixB /= np.sum(UpdatedMatrixB)
    return UpdatedMatrixM, UpdatedMatrixB

# Convergence
def verify_convergence(previous_m, current_m, epsilon=1e-6):
    return np.linalg.norm(previous_m - current_m) < epsilon

# Position Weight Matrix
def generate_pwm(MatrixM, MatrixB):
    log_ratio = np.log((MatrixM + 1e-10) / (MatrixB[:, np.newaxis] + 1e-10))
    return np.where(log_ratio == -np.inf, -1e10, log_ratio)

# Iterations
for i in range(max_iter):
    PreviousMatrixM = MatrixM.copy()
    optimal_positions = expectation_step(input_sequences, MatrixM, motif_width)
    MatrixM, MatrixB = maximization_step(input_sequences, optimal_positions, motif_width)
    if verify_convergence(PreviousMatrixM, MatrixM):
        break

PWM = generate_pwm(MatrixM, MatrixB)
print(f"Updated MatrixM: {MatrixM}")
print(f"Updated MatrixB: {MatrixB}")
print(f"Position Weight Matrix: {PWM}")