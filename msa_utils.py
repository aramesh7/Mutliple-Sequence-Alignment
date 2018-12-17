import numpy as np
import json


def compute_profile(msa, alphabet):
    """

    :param msa: The multiple sequence alignment to be converted into profile
    :param alphabet: The alphabet from which the sequence is constructed
    :return: a dictionary of size
    :runtime:
    """
    if not msa:
        return {}
    align_size = len(msa[0])
    profile = {x: [sum([(1.0 if msa[k][j] == x else 0.0)
                        for k in range(len(msa))
                        ]) / len(msa)
                   for j in range(align_size)
                   ]
               for x in alphabet
               }
    return profile


def compute_sigma_profile_profile(msa1, msa2, scoring, alphabet):
    scoring['-']['-'] = 0
    p = compute_profile(msa1, alphabet)
    q = compute_profile(msa2, alphabet)
    align1_size = len(msa1[0])
    tao1 = {x: [sum([p[y][j] * scoring[x][y]
                     for y in alphabet
                     ])
                for j in range(align1_size)
                ]
            for x in alphabet
            }
    align2_size = len(msa2[0])
    tao2 = {x: [sum([q[y][j] * scoring[x][y]
                     for y in alphabet
                     ])
                for j in range(align2_size)
                ]
            for x in alphabet
            }
    sigma = [[sum([sum([p[x][i] * q[y][j] * scoring[x][y]
                        for y in alphabet])
                   for x in alphabet])
              for j in range(align2_size)]
             for i in range(align1_size)]
    return tao1, tao2, sigma


def traceback(msa1, msa2, backtrace, start_level = 2):
    n = start_level
    i = len(backtrace[0]) - 1
    j = len(backtrace[0][0]) - 1
    new_al1 = [s[:] for s in msa1]
    new_al2 = [w[:] for w in msa2]
    while True:
        n, di, dj = backtrace[n][i][j]
        if not di and dj:
            for seq1 in new_al1:
                seq1.insert(i, '-')
        if not dj and di:
            for seq2 in new_al2:
                seq2.insert(j, '-')
        i, j = i + di, j + dj
        if (i <= 0 and j <= 0):
            break
    new_alignment = new_al1 + new_al2
    return new_alignment


def align_profile_profile(msa1, msa2, scoring, alphabet):
    """
    Aligns two multiple sequence alignments using profile-profile alignment in quadratic
    time and space. This function can also align two sequences or a sequence to a profile.

    :param msa1: A list of lists of characters in the alphabet. Represents the first MSA
         (multiple sequence alignment)
    :param msa2: Same as msa1 but represents the second MSA.
    :param scoring: a dictionary representing the scoring function
    :return:
    """
    # Compute scoring functions
    tao1, tao2, sigma = compute_sigma_profile_profile(msa1, msa2, scoring, alphabet)

    # Affine Gap coefficients
    gap_open = scoring['gap_open']

    # Allocate tables
    S = np.zeros((len(msa1[0]) + 1, len(msa2[0]) + 1))
    X = np.zeros((len(msa1[0]) + 1, len(msa2[0]) + 1))
    Y = np.zeros((len(msa1[0]) + 1, len(msa2[0]) + 1))

    # 3-D grid for back-tracking
    backtrace = [[[(-1, 0, 0) for _ in range(S.shape[1])] for _ in range(S.shape[0])] for _ in range(3)]
    backtrace[0][0][0], backtrace[1][0][0], backtrace[2][0][0] = (0, -1, -1), (1, -1, -1), (2, -1, -1)

    for j in range(S.shape[1]):
        for i in range(S.shape[0]):
            if i == 0 and j == 0:
                continue
            temp_S, temp_X, temp_Y = [], [], []
            # First compute
            if i > 0:
                temp_Y.append(((2, -1, 0), S[i - 1, j] + tao1['-'][i - 1] + gap_open))
            if j > 0:
                temp_X.append(((2, 0, -1), S[i, j - 1] + tao2['-'][j - 1] + gap_open))
            if i > 1:
                temp_Y.append(((1, -1, 0), Y[i - 1, j] + tao1['-'][i - 1]))
            if j > 1:
                temp_X.append(((0, 0, -1), X[i, j - 1] + tao2['-'][j - 1]))

            # Fill in X and Y
            for t, (matrix, temp) in enumerate([(X, temp_X), (Y, temp_Y)]):
                if temp:
                    max_i = np.argmax([temp[x][1] for x in range(len(temp))])
                    matrix[i, j] = temp[max_i][1]
                    backtrace[t][i][j] = temp[max_i][0]

            # Now calculate S
            if i > 0:
                temp_S.append(((1, 0, 0), Y[i, j]))
            if j > 0:
                temp_S.append(((0, 0, 0), X[i, j]))
            if i > 0 and j > 0:
                temp_S.append(((2, -1, -1), S[i - 1, j - 1] + sigma[i - 1][j - 1]))

            max_i = np.argmax([temp_S[x][1] for x in range(len(temp_S))])
            S[i, j] = temp_S[max_i][1]
            backtrace[2][i][j] = temp_S[max_i][0]

    scores = [X[-1, -1], Y[-1, -1], S[-1, -1]]
    start_level = np.argmax(scores)
    score = scores[start_level]
    return score, traceback(msa1, msa2, backtrace, start_level)


def align_sequence_sequence_affine(v, w, delta):
    delta['-']['-'] = -np.inf
    gap_open = delta['gap_open']
    gap_extend = delta['gap_extend']
    # We now have 3 different recurrences
    M = np.zeros((len(v) + 1, len(w) + 1))
    X = np.zeros((len(v) + 1, len(w) + 1))
    Y = np.zeros((len(v) + 1, len(w) + 1))

    # Initialize gap matrices
    X[1:, 0] = -np.inf
    Y[0, 1:] = -np.inf


    # Initialize the backtracking matrix
    backtrace = [[[(-1,0, 0) for _ in range(len(w) + 1)] for _ in range(len(v) + 1)] for _ in range(3)]
    backtrace[0][0][0], backtrace[1][0][0], backtrace[2][0][0] = (0, -1, -1), (1, -1, -1), (2, -1, -1)

    for j in range(len(w) + 1):
        for i in range(len(v) + 1):
            if i == 0 and j == 0:
                continue

            # X[i,j] and Y[i,j] must be calculated first
            temp_M, temp_X, temp_Y = [], [], []
            if i > 0:
                temp_Y.append(((2, -1, 0), M[i - 1, j] + gap_open + gap_extend))
            if j > 0:
                temp_X.append(((2, 0, -1), M[i, j - 1] + gap_open + gap_extend))
            if i > 1:
                temp_Y.append(((1, -1, 0), Y[i - 1, j] + gap_extend))
            if j > 1:
                temp_X.append(((0, 0, -1), X[i, j - 1] + gap_extend))

            for t, (matrix, temp) in enumerate([(X, temp_X), (Y, temp_Y)]):
                if temp:
                    max_i = np.argmax([temp[x][1] for x in range(len(temp))])
                    matrix[i, j] = temp[max_i][1]
                    backtrace[t][i][j] = temp[max_i][0]

            # Now evaluate M[i,j]
            if i > 0:
                temp_M.append(((1, 0, 0), Y[i, j]))
            if j > 0:
                temp_M.append(((0, 0, 0), X[i, j]))
            if i > 0 and j > 0:
                temp_M.append(((2, -1, -1), M[i - 1, j - 1] + delta[v[i - 1]][w[j - 1]]))
            max_i = np.argmax([temp_M[x][1] for x in range(len(temp_M))])
            M[i, j] = temp_M[max_i][1]
            backtrace[2][i][j] = temp_M[max_i][0]

    scores = [X[len(v),len(w)], Y[len(v), len(w)], M[len(v), len(w)]]
    start_level = np.argmax(scores)
    score = scores[start_level]
    return score, traceback([v], [w], backtrace, start_level)


def print_alignment(name, aln, scoring, labels=None, file=None):
    if file:
        json.dump(scoring, file)
        print('', file=file)
    print(name, file=file)
    print('SP-score: %f'%sp_score(aln, scoring), file=file)
    if labels:
        for i,seq in enumerate(aln):
            print(''.join(seq)+ '   <--   ' + labels[i], file=file)
    else:
        for seq in aln:
            print(''.join(seq), file=file)



def get_scoring_function(alphabet, match, mismatch, gap_open, gap_extend):
    size = len(alphabet)
    scoring = {}
    for i in range(size):
        scoring[alphabet[i]] = {alphabet[j]: (match if alphabet[i] == alphabet[j] else mismatch) for j in range(size)}
    scoring['-'] = {x: gap_extend for x in alphabet}
    for x in alphabet:
        scoring[x]['-'] = gap_extend
    scoring['gap_open'] = gap_open
    scoring['gap_extend'] = gap_extend
    return scoring


def sp_score(msa, scoring=None):
    if not scoring:
        alphabet = ['A', 'C', 'T', 'G', '-']
        size = len(alphabet)

        # Create scoring function for sequence alignment
        scoring = {}
        for i in range(size):
            scoring[alphabet[i]] = {alphabet[j]: (1 if alphabet[i] == alphabet[j] else -1) for j in range(size)}
        scoring['gap_open'] = 0
        scoring['gap_extend'] = -1
    scoring['-']['-'] = 2

    sp = 0
    for i in range(len(msa)):
        for j in range(i + 1, len(msa)):
            gap_cost = 0
            gap_opened_i = False
            gap_opened_j = False
            for k in range(len(msa[i])):
                if msa[i][k] == '-' and not gap_opened_i:
                    gap_opened_i = True
                    gap_cost += (scoring['gap_open'] + scoring['gap_extend'])
                elif msa[i][k] == '-':
                    gap_cost += scoring['gap_extend']
                elif gap_opened_i:
                    gap_opened_i = False
                if msa[j][k] == '-' and not gap_opened_j:
                    gap_opened_j = True
                    gap_cost += (scoring['gap_open'] + scoring['gap_extend'])
                elif msa[j][k] == '-':
                    gap_cost += scoring['gap_extend']
                elif gap_opened_j:
                    gap_opened_j = False
            match_score = sum([scoring[msa[i][k]][msa[j][k]] if not ((msa[i][k] == '-') ^ (msa[j][k] == '-')) else 0.0  for k in range(len(msa[i]))])
            sp += (gap_cost + match_score)
    return sp

def sp_score_vanilla(msa, scoring=None):
    if not scoring:
        alphabet = ['A', 'C', 'T', 'G', '-']
        size = len(alphabet)

        # Create scoring function for sequence alignment
        scoring = {}
        for i in range(size):
            scoring[alphabet[i]] = {alphabet[j]: (1 if alphabet[i] == alphabet[j] else -1) for j in range(size)}
        scoring['gap_open'] = 0
        scoring['gap_extend'] = -1
    scoring['-']['-'] = 0
    sp = 0
    for i in range(len(msa)):
        for j in range(i + 1, len(msa)):
            sp += sum([scoring[msa[i][k]][msa[j][k]] for k in range(len(msa[i]))])
    return sp

#X, Y, M, score, new_aln = align_profile_profile([list(words[0])], [list(words[1])], scoring, alphabet)