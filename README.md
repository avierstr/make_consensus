# make_consensus
Makes the consensus sequence of several DNA sequences in a file.

(I only tested this on Linux, probably works on Windows and Mac)
Requirements:
- Python 3
- python3-dev, python3-setuptools, python3-pip, python3-wheel
  (`sudo apt-get install python3-dev python3-setuptools python3-pip python3-wheel`)
- c-implementation of Levenshtein: https://pypi.org/project/python-Levenshtein/
  (`python3 -m pip install python-Levenshtein`)

### Options:

`-i --input`: Input file in fastq of fasta format

 `-s --similar`: Similarity to add a read to the consensus (value between 50 and 100). Default=85.0 for nanopore reads

The script compares every sequence (in forward and reverse direction) with the first sequence in the file.  If a sequence is in reverse order, it takes the complement reverse for the consensus.  By default, the similarity between the sequences has to be >= 85%, otherwise it is excluded.  For this reason, it is important that the first sequence in the file is one that you want for the consensus and not a very different one.

### Command examples:

`python3 make_consensus.py -i inputset1.fasta`

`python3 make_consensus.py --input inputset2.fastq --similar 90`
