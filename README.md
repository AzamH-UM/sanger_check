# Sanger Check: Package for Automating Analysis of Sanger Sequencing Data

## Examples:
See `examples/sanger_check_demo.ipynb` for introduction to using sanger_check.

See `https://colab.research.google.com/drive/1HsShDkgk7yC9oluZwzf2G2QnpACqnUzF?usp=sharing` for interactive colab notebook.

## Availible Functions
`sc.parse_fasta(fasta_file)`
 - Parses a Fasta file and returns a pandas DataFrame

`sc.auto_assign(sanger_df, seq_df)`
 - Parsed multifasta for sanger sequences and dna sequences can be automatically matched together using pairwise alignment
 - Returns a dataframe that maps sanger sequence names to reference dna and DNA orientation (forward/reverse)

`sc.check_sanger_sequence(dna_seq_name, dna_seq sanger_seq_name, sanger_seq, plot=False)`
 - Prints information and plots alignment between reference sequence and sanger sequence