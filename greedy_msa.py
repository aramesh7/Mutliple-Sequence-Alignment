import numpy as np

def compute_profile(msa, alphabet):
    if not msa:
        return {}
    align_size = len(msa[0])
    profile = {x:[sum([(1.0 if msa[k][j] == x else 0.0) 
                        for k in range(len(msa))
                    ])/len(msa)
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
    tao1 = {x:[sum([p[y][j]*scoring[x][y] 
                    for y in alphabet
                ]) 
                for j in range(align1_size)
            ] 
            for x in alphabet
        }
    align2_size = len(msa2[0])
    tao2 = {x:[sum([q[y][j]*scoring[x][y] 
                    for y in alphabet
                ]) 
                for j in range(align2_size)
            ] 
            for x in alphabet
        }
    sigma = [[sum([sum([p[x][i]*q[y][j]*scoring[x][y] 
                        for y in alphabet]) 
                    for x in alphabet]) 
                for j in range(align2_size)] 
            for i in range(align1_size)]
    return tao1, tao2, sigma 


def traceback(msa1, msa2, backtrace):
    i = len(backtrace)-1
    j = len(backtrace[0])-1
    new_al1 = [s[:] for s in msa1]
    new_al2 = [w[:] for w in msa2]
    while True:
        di, dj = backtrace[i][j]
        if not di:
            for seq1 in new_al1:
                seq1.insert(i,'-')
        if not dj:
            for seq2 in new_al2:
                seq2.insert(j,'-')
        i,j = i+di, j+dj
        if (i<=0 and j<=0):
            break
    new_alignment = new_al1 + new_al2
    return new_alignment

def align_profile_profile(msa1, msa2, scoring, alphabet):
    """
    Aligns two multiple sequence alignments using profile-profile alignment in quadratic
    time and space. 

    @msa1 A list of lists of characters in the alphabet. Represents the first MSA 
         (multiple sequence alignment)
    @msa2 Same as msa1 but represents the second MSA.

    @scoring a dictionary representing the scoring function

    """
    tao1, tao2, sigma = compute_sigma_profile_profile(msa1, msa2, scoring, alphabet)
    S = np.zeros((len(msa1[0])+1,len(msa2[0])+1))
    backtrace = [[(0,0) for j in range(S.shape[1])] for i in range(S.shape[0])]
    backtrace[0][0] = (-1,-1)

    for j in range(S.shape[1]):
        for i in range(S.shape[0]):
            if i == 0 and j == 0:
                continue
            temp = []
            if i > 0:
                temp.append(((-1,0), S[i-1][j] + tao1['-'][i-1]))
            if j > 0:
                temp.append(((0,-1), S[i][j-1] + tao2['-'][j-1]))
            if i > 0 and j > 0:
                temp.append(((-1,-1), S[i-1][j-1] + sigma[i-1][j-1]))
            max_i = np.argmax([temp[x][1] for x in range(len(temp))])
            S[i,j] = temp[max_i][1]
            backtrace[i][j] = temp[max_i][0]
    score = S[S.shape[0]-1, S.shape[1]-1]
    return score, traceback(msa1, msa2, backtrace)


def print_alignment(name, aln):
    print name, aln[0], sp_score(aln[1], None)
    for seq in aln[1]:
        print seq

def get_scoring_function(alphabet, match, mismatch, gap_open, gap_extend):
    size = len(alphabet)
    scoring  = {}
    for i in range(size):
        scoring[alphabet[i]] = {alphabet[j]:(match if alphabet[i] == alphabet[j] else mismatch) for j in range(size)}
    scoring['-'] = {x: gap_open for x in alphabet}
    for x in alphabet:
        scoring[x]['-'] = gap_open
    return scoring, gap_extend

def sp_score(msa,scoring):
    if not scoring:
        alphabet = ['A', 'C', 'T', 'G', '-']
        size = len(alphabet)

        # Create scoring function for sequence alignment
        scoring = {}
        for i in range(size):
            scoring[alphabet[i]] = {alphabet[j]:(1 if alphabet[i] == alphabet[j] else -1) for j in range(size)}
    scoring['-']['-'] = 0
    sp = 0
    for i in range(len(msa)):
        for j in range(i+1, len(msa)):
            sp += sum([scoring[msa[i][k]][msa[j][k]] for k in range(len(msa[i]))])
    return sp

alphabet = ['A', 'C', 'T', 'G', '-']

# Create scoring function for sequence alignment
scoring, gap_extend = get_scoring_function(alphabet, 1, -1, -1, 0)

# Do pairwise alignments with Needleman-Wunsch
words = list(map(list, ['GATTCA', 'GTCTGA', 'GATATT','GTCAGC']))
aligned_inputs = []
for i in range(len(words)):
    for j in range(i+1, len(words)):
        aligned_inputs.append(align_profile_profile([list(words[i])], [list(words[j])], scoring, alphabet))
for aln in aligned_inputs:
    print_alignment('Intermediate Pairs: ', aln)

# find the closest pair of sequences to form a profile
max_i = np.argmax([aligned_inputs[i][0] for i in range(len(aligned_inputs))])
most_similar_pair = aligned_inputs[max_i]

# Remove these sequences from consideration
for word in most_similar_pair[1]:
    words.remove(word)

print_alignment('Most Similar Pair: ',most_similar_pair)
alignments = [most_similar_pair[1]]
seq_profile_aligns = []

# Now get the score of the remaining two sequences
aligned_inputs = []
for i in range(len(words)):
    for j in range(i+1, len(words)):
        aligned_inputs.append(align_profile_profile([list(words[i])], [list(words[j])], scoring, alphabet))

next_best_pair = aligned_inputs[0]
print_alignment('Next Pair: ', next_best_pair)

# Now align every sequence-profile pair of 
for alignment in alignments:
    for seq in words:
        seq_profile_aligns.append(align_profile_profile(alignment, [seq], scoring, alphabet))

for aln in seq_profile_aligns:
    print_alignment('Intermediate sequence-profile alignments', aln)

final_alns = []
for i in range(len(seq_profile_aligns)):
    final_alns.append(align_profile_profile(seq_profile_aligns[i][1], [words[1-i]], scoring, alphabet))

for aln in final_alns:
    print_alignment('Final sequence-profile alignment', aln)

# Align the two profiles
final_aln = align_profile_profile(alignments[0], next_best_pair[1], scoring, alphabet)
print_alignment('Final Profile-Profile alignment', final_aln)
print sp_score(final_aln[1], scoring)
# Compute the scoring function for profile profile alignme