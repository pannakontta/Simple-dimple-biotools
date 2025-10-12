# Simple-dimple-biotools
**Simple-dimple-biotools** is a tool that allows you to perform simple procedures on with different biological data.
## Functional description
File main.py contains head functions: `run_dna_rna_tools` - for working with DNA and RNA sequences and `filter_fastq` - for working with fastq-sequences
File bio_files_processor.py contain functions: `convert_multiline_fasta_to_oneline` and `parse_blast_output`.

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


**filter_fastq**

The function allows you to select all sequences that satisfy the specified conditions. It takes 4 arguments as input:
+ `input_fastq` - is name of the standard `.fastq` file to be processed
+ `output_fastq` - is name of the standard `.fastq` file with relevant secuences (default - `output_file.fastq`)
+ `gc_bounds` - is the GC-content interval (in percent) for filtering (by default, it is (0, 100)
+ `length_bounds` - the length interval for filtering
+ `quality_threshold` - the threshold value of the average read quality for filtering is 0 by default (phred33 scale). Reads with average quality for all nucleotides below the threshold are discarded.

**convert_multiline_fasta_to_oneline**

The function get fasta-file, where the sequence can be split into several lines. The function joins these lines into single line. It takes 2 arguments as input:

+ `input_fasta` - is name of the standard `.fasta` file to be processed
+ `output_fasta` - is name of the standard `.fasta` file with sequences after joining (default - `output_fasta.fasta`)

**parse_blast_output**

The function get information parsed from the blast in txt-format and find the bestmatches with the database and save all the results to a file. It takes 2 arguments as input:

+ `input_file` - is name of the `.txt` file to be processed
+ `output_file` - is name of the `.txt` file with Description of the best matches (default - `best-proteins.txt`)

## Example of usage
```
run_dna_rna_tools("ATGCAG", "complement")
run_dna_rna_tools("ata", "UAG", "ATGCAG", "reverse")
```
```
filter_fastq("example_fastq.fastq", "example_output.fastq")
filter_fastq("example_fastq.fastq", "example_output.fastq", (30, 60), (70, 80), 33)
```
