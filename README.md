# Multiple Sequence Alignment
Requirements
------------

All the requirements can be installed via `pip` or `pip3` by running:
```
> pip3 install -r requirements.txt
```
from the root ```Multiple-Sequence-Alignment``` directory.

Instructions
------------

In order to run Greedy Progressive Multiple Sequence Alignment, we need to provide
the script with 2 positional arguments.

1. An `inputs.txt` - a text file where each line represents the filepath
of a `.txt` or `fasta` format file. There is a provided `example_input.txt`.

2. A scoring function as either a text file input to the `--scoring` argument or as 
values to each of the following arguments.
    - `--match` : takes a float value to score matches in the sequences.
    - `--mismatch` : takes a float value to score mismatches.
    - `--affine` : If this argument is present, alignments will be scored using affine 
    gaps, i.e. gap open penalties and gap extension penalties are separated. This argumen
    will require the next two arguments as well.
    - `--gap-open` : The penalty for opening a gap if `--affine` is specified.
    - `--gap-extend` : Gap extension penalty if `--affine` is specified.
    
Optionally, we can provide the following arguments: 

1. An `alphabet` with `-alphabet`, `--alphabet` flag. Values passed to this must be amongst ones currenly
supported (V1: [`DNA`,`RNA`]). This value *dna* by default and represents **[A,C,T,G,-]**.

2. An `output.txt` with `-o`, `--output`, or `--output-file`: The file to write out the final multiple sequence alignment as well as the sp-score
under the given scoring function. The code will dump the provided scoring matrix as a json onto the first line
of the output file, followed by the title, the SP-score of the final alignment then the individual sequences in 
the alignment.  
    
We can now run the script with 

```
> python3 greedy_msa.py examples/example_inputs.txt example_scoring.csv --affine
```

You can say all modes of operation by running:

```
> python3 greedy_msa.py -h
```

- Other flags:
    - `-seq-start`, `--seq-start`: Allows you to specify the start index of the segment of the sequences to be aligned. `0` by default.
    - `-seq-end`, `--seq-end`: Allows you to specify the end index of the segment of the sequences to be aligned. `100` by default.
    - `-gen-seq-list`, `--gen-seq-list` : allows you to create a FASTA formatted .txt file of the segments 
        of sequences you are aligning. This is useful when cross-testing the program with other programs
        like CLUSTAL or MUSCLE.