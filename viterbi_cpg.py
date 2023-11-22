from setup import *

initial_probs = {'notCpG': 0.5, 'CpG': 0.5}

transition_probs = {
    'notCpG': {'notCpG': 0.9999903799268958, 'CpG': 9.620073104211598e-06},
    'CpG': {'notCpG': 0.0012798320652894705, 'CpG': 0.9987201679347105}
}

emission_probs = {
    'notCpG': {'A': 0.28681602060341194, 'C': 0.21308951118173802, 'G': 0.21285794813745923, 'T': 0.2872365200773908},
    'CpG': {'A': 0.14731887312907363, 'C': 0.3483368693986412, 'G': 0.3449527280473017, 'T': 0.1593915294249835}
}

# Hypothetical sample
sequence_sample = 'ACGT' * 20000

def run_viterbi(seq, state_list, init_probs, trans_probs, emit_probs):
    v_matrix = np.zeros((len(state_list), len(seq)))
    paths = {state: [] for state in state_list}

    for st in state_list:
        v_matrix[state_list.index(st), 0] = init_probs[st] * emit_probs[st][seq[0]]
        paths[st] = [st]

    for idx in range(1, len(seq)):
        temp_path = {}
        for st in state_list:
            (max_prob, max_state) = max(
                [(v_matrix[state_list.index(prev_st), idx-1] * trans_probs[prev_st][st] * emit_probs[st][seq[idx]], prev_st) for prev_st in state_list])
            v_matrix[state_list.index(st), idx] = max_prob
            temp_path[st] = paths[max_state] + [st]

        paths = temp_path

    (max_prob, max_state) = max([(v_matrix[state_list.index(st), len(seq) - 1], st) for st in state_list])
    return max_prob, paths[max_state]

# Viterbi algorithm
state_list = ['CpG', 'notCpG']
max_prob, optimal_path = run_viterbi(sequence_sample, state_list, initial_probs, transition_probs, emission_probs)

# Get CpG Island Regions
cpg_islands = []
cpg_start = None
for idx, st in enumerate(optimal_path):
    if st == 'CpG':
        if cpg_start is None: cpg_start = idx
    elif cpg_start is not None:
        cpg_islands.append([cpg_start, idx - 1])
        cpg_start = None
if cpg_start is not None: cpg_islands.append([cpg_start, len(optimal_path) - 1])

print(cpg_islands)
