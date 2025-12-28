# Simple-dimple-biotools
**Simple-dimple-biotools** is a tool that allows you to perform simple procedures on DNA, RNA, and fastq-sequences.
## Functional description
File main.py contains head functions: `run_dna_rna_tools` - for working with DNA and RNA sequences and `filter_fastq` - for working with fastq-sequences

**run_dna_rna_tools**

The function get any number of sequences and the name of one of the foloowinf procedures as input:
+ `is_nucleic_acid` - returns the Boolean result of the sequence check.
+ `transcribe` — return the transcribed sequence
+ `reverse` — return the expanded sequence
+ `complement` — return the complementary sequence
+ `reverse_complement` — return the reverse complementary sequence

<ins>Features</ins>
+ Regardless of the chosen procedure, the tool verifies the validity of the sequences. If the sequence is not DNA or RNA, the tool will show the positions of these sequences.
+ Any number of sequences can be input.
+ Only one procedure can be performed at a time.


**run_dna_rna_tools**

The function allows you to select all sequences that satisfy the specified conditions. It takes 4 arguments as input:
+ `seqs` - is a dictionary consisting of fastq sequences, where
  the `key` is the name of the sequence, the `value` is a tuple consisting of two strings: sequence and quality.
Example:
```
EXAMPLE_FASTQ = {
    # 'name' : ('sequence', 'quality')
    '@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079803': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD')
    }  
```
+ `gc_bounds` - is the GC-content interval (in percent) for filtering (by default, it is (0, 100)
+ `length_bounds` - the length interval for filtering
+ `quality_threshold` - the threshold value of the average read quality for filtering is 0 by default (phred33 scale). Reads with average quality for all nucleotides below the threshold are discarded.

## Example of usage
```
run_dna_rna_tools("ATGCAG", "complement")
run_dna_rna_tools("ata", "UAG", "ATGCAG", "reverse")
```
```
filter_fastq(EXAMPLE_FASTQ)
filter_fastq(EXAMPLE_FASTQ, (30, 60), (70, 80), 33)
```
