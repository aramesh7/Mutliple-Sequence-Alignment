import msa_utils as utils
import numpy as np
import argparse
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt

def progressive_alignment(sequences, scoring, alphabet, show_tree=False):

    # Initialize a distance matrix
    n = len(sequences)
    D = np.zeros((n,n))

    # Fill the pairwise distance matrix
    keys = list(sequences.keys())
    for i in range(n):
        for j in range(i):
            D[i,j], alignment = utils.align_profile_profile([list(sequences[keys[i]])], [list(sequences[keys[j]])], scoring, alphabet)
            D[j,i] = D[i,j]

    # Create Linkage
    Z = linkage(D,method='complete')

    # Create actual alignment
    clusters = [[list(sequences[k])] for k in keys]
    labels = [[k] for k in keys]
    for i in range(Z.shape[0]):
        score, alignment = utils.align_profile_profile(clusters[int(Z[i][0])], clusters[int(Z[i][1])], scoring, alphabet)
        clusters.append(alignment)
        labels.append(labels[int(Z[i][0])] + labels[int(Z[i][1])])

    fig = None
    # Display guide tree
    if show_tree:
        fig = plt.figure(figsize=(10,10))
        dendrogram(Z, labels=keys, leaf_rotation=90)
        plt.xlabel("Samples")
        plt.ylabel("Distance/Scores")
        plt.title("Guide Tree for Given Sequence Data")
        plt.tight_layout()
        plt.show()
    return clusters[-1], labels[-1], fig

def main(args, scoring, alphabet):
    input_file = args.input
    # Read in files
    sequences = {}
    with open(input_file) as f:
        for i,filename in enumerate(f.readlines()):
            with open(filename.rstrip()) as seq_file:
                sequences[filename.rsplit('.', 1)[0].split('/')[-1]] = seq_file.read().replace('\n', '')[args.seq_start:args.seq_end]

    if args.gen_seq_list:
        with open('greedy_msa_intermediate.txt', 'w') as f:
            for name, seq in sequences.items():
                f.write('>' + name + '\n' + seq + '\n')


    # Run the main functions
    aln, labels, figure = progressive_alignment(sequences, scoring, alphabet, args.show_guide_tree)
    utils.print_alignment("Final Alignment:", aln, scoring, labels=labels)
    if args.output:
        with open(args.output, 'w') as f:
            utils.print_alignment("Final Alignment:", aln, scoring, labels=labels, file=f)

    if args.save_guide_tree and figure:
        figure.savefig(args.save_guide_tree)

if __name__ == "__main__":

    # Parse args and call main
    parser = argparse.ArgumentParser(description='A script that takes an input of 2 or more sequences and outputs a pogressive multiple sequence alignment')
    parser.add_argument('input', help='file path for .txt of list of input sequences')
    parser.add_argument('-o', '--output', '--output-file', help='file path for writing the output')
    parser.add_argument('-affine', '--affine', help='Specifies that affine gap penalties must be used when aligning sequences', action='store_true')
    parser.add_argument('scoring', help='file specifying the scoring matrix values as a .csv', default=None)

    parser.add_argument('-match', '--match', help='value overwrites the score for a character match during alignment', type=float)
    parser.add_argument('-mismatch', '--mismatch', help='value overwrites the score for a character mismatch during alignment', type=float)
    parser.add_argument('-gap-open', '--gap-open',
                        help='overwrites the affine penalty coefficient for opening a gap during alignment', type=float)
    parser.add_argument('-gap-extend', '--gap-extend',
                        help='overwrites the affine penalty coefficient for opening a gap during alignment', type=float)
    parser.add_argument('-alphabet', '--alphabet', help='The sequence alphabet. Must choose one of the supported options', choices=['dna','rna'], default='dna')
    parser.add_argument('-seq-start', help='The starting position in the sequences provided for segment alignment', type=int, default=0)
    parser.add_argument('-seq-end', help='The ending position in the sequences provided for segment alignment', type=int, default=100)
    parser.add_argument('-gen-seq-list', '--gen-seq-list', help='Generates an intermediate FASTA format file with all sequence segments being aligned', action='store_true')
    parser.add_argument('-show-tree', '--show-guide-tree', help='When present, the code will generate a dendrogram of the clustered sequences with filenames as labels', action='store_true')
    parser.add_argument('-save-guide-tree', '--save-guide-tree', help='Save the dendrogram of the clustered sequences to the filename provided with this flag')
    args = parser.parse_args()
    print(args)
    alphabet_map = {
        'dna': ['A', 'C', 'G', 'T', '-'],
        'rna': ['A', 'C', 'G', 'U', '-']
    }
    alphabet = alphabet_map[args.alphabet]
    with open(args.scoring) as f:
        [match, mismatch, gap_open, gap_extend] = list(map(float, f.readlines()[1].split(',')))
        scoring = utils.get_scoring_function(alphabet, match, mismatch, gap_open, gap_extend)

    main(args, scoring, alphabet)

