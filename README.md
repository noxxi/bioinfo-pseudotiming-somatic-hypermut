This is the program which was used to generate the figures 5B and S3A in the
paper "Regulation of the germinal center reaction and somatic hypermutation
dynamics by homologous recombination" by Gianna Hirth, Carl-Magnus Svensson,
Katrin Böttcher, Steffen Ullrich, Marc Thilo Figge and Berit Jungnickel.

Input files are in FASTA format. The first sequence is used as a consensus
while the following sequences are variants of this consensus sequence
containing mutations. All sequences have the same length and there are no
indels.
There are two kinds of sequences: control sequence  (Ctrl*.txt) and experimental
data (Brca2*.txt).

We count a) the total number of mutations compared to consensus and b) the kind
of mutations.
The result will be graphics which show the mutation frequency for a specific
kind of mutation after a specific "pseudo-time", which is equivalent to the
total number of mutations in a sequence.

Usage:

    $ python3 relbase.py --ctrl in/Ctrl_* --exp in/Brca2_*

    Resulting output in out/:
      RelBase_AC.png, RelBase_AT.png ... - specific kind of mutation
      RelBase.png - all kinds of mutations combined in one image


