from setup import *

alphabet = ["A", "C", "G", "T"]

### GibbsSampler:
### INPUTS:	S - list of sequences
###		    L - length of motif
###	OUTPUT:	PWM - 4 x L list with frequencies of each base at each position
###               Order of bases should be consistent with alphabet variable
def GibbsSampler(S, L):
    PWM = [[0 for _ in range(L)] for _ in range(len(alphabet))]

    motif_positions = [random.randint(0, len(S[i]) - L) for i in range(len(S))]

    num_iterations = 1000
    for _ in range(num_iterations):
        excluded_sequence_index = random.randint(0, len(S) - 1)
        excluded_sequence = S[excluded_sequence_index]
        updated_PWM = [[0 for _ in range(L)] for _ in range(len(alphabet))]
        for i in range(len(S)):
            if i != excluded_sequence_index:
                for j in range(L):
                    base_index = alphabet.index(S[i][motif_positions[i] + j])
                    updated_PWM[base_index][j] += 1
        scores = []
        for pos in range(len(excluded_sequence) - L + 1):
            score = 1.0
            for j in range(L):
                base_index = alphabet.index(excluded_sequence[pos + j])
                score *= updated_PWM[base_index][j]
            scores.append(score)
        if all(score == 0 for score in scores):
            sampled_position = random.randint(0, len(excluded_sequence) - L)
        else:
            total_score = sum(scores)
            sampled_position = random.choices(
                range(len(scores)),
                weights=[score / total_score if total_score > 0 else 1 / len(scores) for score in scores]
            )[0]
        motif_positions[excluded_sequence_index] = sampled_position
        for i in range(len(S)):
            for j in range(L):
                base_index = alphabet.index(S[i][motif_positions[i] + j])
                PWM[base_index][j] += 1
        for j in range(L):
            col_sum = sum(PWM[i][j] for i in range(len(alphabet)))
            for i in range(len(alphabet)):
                PWM[i][j] /= col_sum
    return PWM

def get_motif_seq(PWM):
    motif_seq = ""
    for i in range(len(PWM[0])):
        comp = [row[i] for row in PWM]
        ind = max(range(len(comp)), key = comp.__getitem__)
        motif_seq += alphabet[ind]
    return motif_seq

def print_PWM(PWM, L):
    PWM_comp = []
    PWM_comp.append("||" + "|".join([str(i) for i in range(1, L + 1)]) + "|")
    PWM_comp.append("|-" * (L + 1) + "|")
    for i in range(4):
        PWM_comp.append("|" + alphabet[i] + "|" + "|".join([str(round(val, 2)) for val in PWM[i]]) + "|")
    print("\n".join(PWM_comp) + "\n")

def print_logo(PWM):
    import logomaker
    PWM = np.array(PWM).T
    PWM = pd.DataFrame(PWM, columns = alphabet)
    logomaker.Logo(PWM, fade_below = 0.8)

def run_GibbsSampler(S, L, n):
    motif_seqs, motif_PWM = {}, {}
    for i in range(n):
        PWM = GibbsSampler(S, L)
        motif_seq = get_motif_seq(PWM)
        motif_seqs.setdefault(motif_seq, 0)
        motif_seqs[motif_seq] += 1
        motif_PWM[motif_seq] = PWM
    best_motif = max(motif_seqs.keys(), key = lambda x: motif_seqs[x])
    print("Most consistent motif: \n" + best_motif + "\n")
    print("PWM (paste into text block):")
    print_PWM(motif_PWM[best_motif], L)
    print("Sequence logo:")
    print_logo(motif_PWM[best_motif])