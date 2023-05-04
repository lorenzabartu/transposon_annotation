# transposon_annotation
Final project for Computational Genomics

# Requirements
Biopython 
- pip install biopython
DNA features viewer
- pip install dna-features-viewer

# Usage 
python find_transposon_final.py -f plasmid.fasta -o plasmid 

# Optional flags:
-i: sets cutoff for percent identity (default is 0.85)
-c: blasts the genetic context around mobile genetic elements
-m: constructs plasmid map

# Inputs
FASTA file of bacterial plasmid sequence

# Find_transposon_final.py

Using a blastn approach, this script blasts the inputted sequence against a transposon and insertion sequence (IS) database. Following the blast, the results are interpreted using features like percent identity, location, e-value, and subject ID to construct expected mobile genetic element maps. 

# Outputs
- GenBank of plasmid annotated with found mobile genetic elements
- CSV of mobile genetic elements found (Subject ID, Start, Stop)
- Plasmid map of the mobile genetic elements located on plasmid
