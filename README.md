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


