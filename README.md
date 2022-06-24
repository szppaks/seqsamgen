# seqsamgen
DNA sequence generator for AI training


Sequence generator for AI training
required PYTHON modules :biopython, seqlogo, pandas,numpy

It has two usage modes: command line mode with paramter file or gui mode
command line mode: python seqsamgen.py  [controlfile.cfg]
GUI mode:          python seqsamgen.py

You can set: motif file, order, gap, seq length, number of sequenties, output directory name (if not set, it generated randomly)
Known motif files: JASPAR PFM, HOMER PPM
Output format: FASTA
The control (negative) set generated automically.

Rules of generating
Motif number  |   ordered    | Gaps       |            Control set
1             |     -        |    -       |  Random
N             |     N        |    N       |  Random, N-1 motif
N             |     Y        |    N       |  Random, N-1 motif, swapped
N             |     N        |    Y       |  Random, N-1 motif, different gap
N             |     Y        |    Y       |  N-1 motif, different gap, swapped
Different gap: if you give gaps: 1,3,5 the different gap :2,4,6.. maxpossiblegaplen


## seqsamgen config file sample
[sequence generating]
#the number of the samples
seqnum=1000
#MOTIF file name (JASPAR PFM or HOMER PPM format)
motiffilename=example.motif
#If you want to use only one motif from file, vou can set the name.
#motifname=
#Gap sizes separated with comma.
gaps=1,3,5
#The length of the generated sequences
sequence length=120
#Motifs are ordered or not Y/N
ordered=Y
[other]
#Verbose mode: set Y if qou want to print the generated motifs
verbose mode=N
#output directory name. If not set, it will be generated automatically
output name=seqoutdir
