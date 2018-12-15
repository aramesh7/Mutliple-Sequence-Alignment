#!/usr/bin/python
import msa_utils as utils
import numpy as np
import getopt
import sys

# Define the sequence alphabet
alphabet = ['A', 'C', 'T', 'G', '-']

# Create scoring function for sequence alignment
scoring = utils.get_scoring_function(alphabet, 1, -1, 0, -1)

# Do pairwise alignments with Needleman-Wunsch
words = list(map(list, ['GATTCA', 'GTCTGA', 'GATATT', 'GTCAGC']))
aligned_inputs = []
for i in range(len(words)):
    for j in range(i + 1, len(words)):
        aligned_inputs.append(utils.align_profile_profile([list(words[i])], [list(words[j])], scoring, alphabet))
for aln in aligned_inputs:
    utils.print_alignment('Intermediate Pairs: ', aln, scoring)

# find the closest pair of sequences to form a profile
max_i = np.argmax([aligned_inputs[i][0] for i in range(len(aligned_inputs))])
most_similar_pair = aligned_inputs[max_i]

# Remove these sequences from consideration
for word in most_similar_pair[1]:
    words.remove(word)

utils.print_alignment('Most Similar Pair: ', most_similar_pair, scoring)
alignments = [most_similar_pair[1]]
seq_profile_aligns = []

# Now get the score of the remaining two sequences
aligned_inputs = []
for i in range(len(words)):
    for j in range(i + 1, len(words)):
        aligned_inputs.append(utils.align_profile_profile([list(words[i])], [list(words[j])], scoring, alphabet))

next_best_pair = aligned_inputs[0]
utils.print_alignment('Next Pair: ', next_best_pair, scoring)

# Now align every sequence-profile pair of
for alignment in alignments:
    for seq in words:
        seq_profile_aligns.append(utils.align_profile_profile(alignment, [seq], scoring, alphabet))

for aln in seq_profile_aligns:
    utils.print_alignment('Intermediate sequence-profile alignments', aln, scoring)

final_alns = []
for i in range(len(seq_profile_aligns)):
    final_alns.append(utils.align_profile_profile(seq_profile_aligns[i][1], [words[1 - i]], scoring, alphabet))

for aln in final_alns:
    utils.print_alignment('Final sequence-profile alignment', aln, scoring)

# Align the two profiles
final_aln = utils.align_profile_profile(alignments[0], next_best_pair[1], scoring, alphabet)
utils.print_alignment('Final Profile-Profile alignment', final_aln, scoring)
# Compute the scoring function for profile profile alignment
#words = ['AATCGAGTCG', 'GTACGATG']


def main(argv):







if __name__ == "__main__":
    main(sys.argv[1:])
