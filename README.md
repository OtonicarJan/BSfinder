# BSfinder
Find a potential binding sites for a transcription factor in a small phage genome (<100000bp).

Input files: fasta file of a genome and list of known binding sites.
Output file: file with scores for binding sites. Cutoff is already set to the top 80%. 

You must have installed numpy and biopython. 
To install them: 

```pip3 install numpy biopython```

To run BSfinder: 

```python3 BSfinder.py genome.fasta list_of_binding_sites.txt outfile.txt```
