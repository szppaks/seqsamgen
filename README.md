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
<i>
[sequence generating]
#the number of the samples<br>
seqnum=1000<br>
#MOTIF file name (JASPAR PFM or HOMER PPM format)<br>
motiffilename=example.motif<br>
#If you want to use only one motif from file, vou can set the name.<br>
#motifname=<br>
#Gap sizes separated with comma.<br>
gaps=1,3,5<br>
#The length of the generated sequences<br>
sequence length=120<br>
#Motifs are ordered or not Y/N<br>
ordered=Y<br>
[other]<br>
#Verbose mode: set Y if qou want to print the generated motifs<br>
verbose mode=N<br>
#output directory name. If not set, it will be generated automatically<br>
output name=seqoutdir<br>
