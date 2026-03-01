# Simple-dimple-biotools

Simple-dimple-biotools is a small collection of utilities for working with
biological sequences (DNA, RNA and amino-acid sequences) and for
filtering FASTQ reads.

## Functional description

`main.py` script provides core sequence data structures and a FASTQ filtering utility. It defines a compact object model for
biological sequences and a convenience function for filtering FASTQ files.

**Key classes**:
  - `BiologicalSequence` — abstract base class defining the sequence interface.
  - `ManageableBiologicalSequence` — concrete base with common implementations
    for length, indexing, iteration and membership.
  - `NucleicAcidSequence` — base for DNA/RNA sequences with `reverse`,
    `complement` and `reverse_complement` methods.
  - `DNASequence` — DNA-specific subclass that validates absence of `U` and
    provides `transcribe()` to produce RNA.
  - `RNASequence` — RNA-specific subclass that validates absence of `T`.
  - `AminoAcidSequence` — protein sequence class with validation and helper
    methods such as `get_start_codon_idxs()` to locate methionine residues in zero-based format.


**filter_fastq function**:
  reads a FASTQ file using Biopython, filters records by GC content, length and average PHRED quality, and writes the passing records to `output_fastq`.

  The function takes 5 arguments as input:
    
  + `input_fastq` - is name of the standard `.fastq` file to be processed
  + `output_fastq` - is name of the standard `.fastq` file with relevant secuences (default - `output_file.fastq`)
  + `gc_bounds` - is the GC-content interval (in percent) for filtering (by default, it is (0, 100)); can be a number as the max bound
  + `length_bounds` - the length interval for filtering; can be a number as the max bound
  + `quality_threshold` - the threshold value of the average read quality for filtering is 0 by default (phred33 scale). Reads with average quality for all nucleotides below the threshold are discarded.

## Usage example

```py
from main import filter_fastq, DNASequence

# Filter FASTQ file
filter_fastq('reads.fastq', 'filtered.fastq', gc_bounds=(30,70), length_bounds=(50,300), quality_threshold=20)

# Work with sequences
dna = DNASequence('ATGCGT')
print(dna.transcribe())  # outputs RNA sequence
> AUGCGU

aaseq = AminoAcidSequence('QMEGHILMK')
print(aaseq.get_start_codon_idxs())
> [1, 7]
```

