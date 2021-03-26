Oriented TNF
============

This project uses the *oriented* TNF to detect and analyze inversions in bacterial genomes. The orientation matrix is the pairwise orientations between small windows evenly spaced along the genome. A simple spectral clustering method is used to find inversions in the matrix. 

Available features (all ran by default):

1) Compute and plot Orientation Matrix heatmap.

2) Locate inversions in the input genome. If only one large inversion exists, it corresponds to the Origin/Terminus replication sites.

3) Detect and correct an inverted misassembly.

Requirements
-------------

- Python 3.5+
- MUMmer 3+ (with nucmer)

Check **requirements.txt** for python package requirements or run 

``python3 -m pip install -r requirements.txt``

Run code
------------

Must be run in the Oriented-TNF directory! Description of inputs and options can be found using
    
``python3 detect.py -h``

Example on preinstalled  *L. borgpetersenii* genome (assuming nucmer executable is in PATH):

``python3 detect.py --window 50000 \
                    --stride 25000 \
                    --pad 10000 \
                    --genome data/L_borgpetersenii.fasta.gz \
                    --out data/test``
