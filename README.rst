Oriented TNF
============

This project uses the *oriented* TNF to detect and analyze inversions in bacterial genomes. The orientation matrix is the pairwise orientations between small windows evenly spaced along the genome. A simple spectral clustering method is used to find inversions in the matrix. 

Available features (all ran by default):

1) Compute and plot Orientation Matrix heatmap.

2) Locate inversions in the input genome. If only one large inversion exists, it corresponds to the Origin/Terminus replication sites.

3) Detect and correct an inverted misassembly.


Run code
------------

Python 3.5 or newer is required.

Check **requirements.txt** for package requirements or run 

``python3 -m pip install -r requirements.txt``

Description of inputs and options can be found using
    
``python3 detect.py -h``
