from typing import Union
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from numbers import Number


class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.
    
    This class defines the interface that all biological sequence types must implement.
    """

    @abstractmethod
    def __len__(self) -> int:
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __iter__(self):
        pass

    @abstractmethod
    def __contains__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class ManageableBiologicalSequence(BiologicalSequence):
    """
    Concrete implementation of BiologicalSequence with basic sequence operations.
    
    This class provides common functionality for managing biological sequences,
    including length, indexing, string conversion, and iteration.
    """
    
    def __init__(self, seq: str):
        self.seq = seq

    def __len__(self):
        return len(self.seq)
    
    def __getitem__(self, index):
        if isinstance(index, int) or isinstance(index, slice):
            return self.seq[index]

    def __str__(self):
        return self.seq
    
    def __iter__(self):
        return iter(self.seq)

    def __contains__(self, item):
        return item in self.seq
    
    def _check_alphabet(self, alphabet):
        return set(self.seq.upper()).issubset(alphabet)


class NucleicAcidSequence(ManageableBiologicalSequence):
    def __init__(self, seq):
        super().__init__(seq)
        alphabet = 'ATGCU'
        if not self._check_alphabet(alphabet):
            raise ValueError("Invalid characters in the sequence. Only A, T, G, C, U are allowed.")

    def reverse(self):
        return self.seq[::-1]

    def complement(self):
        if 'U' in self.seq.upper():
            complement_pairs = str.maketrans("AUGCaugc", "UACGuacg")
        else:
            complement_pairs = str.maketrans("ATGCatgc", "TACGtacg")
        return self.seq.translate(complement_pairs)
    
    def reverse_complement(self):
        complementary_seq = self.complement()
        return complementary_seq[::-1]


class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence composed of A, T, G, C nucleotides.

    Inherits from NucleicAcidSequence and provides a transcription method to
    convert DNA into its RNA counterpart.

    Such methods as reverse, complement, and reverse_complement are inherited from NucleicAcidSequence.
    """

    def __init__(self, seq):
        super().__init__(seq)
        if 'U' in self.seq.upper():
            raise ValueError("RNA instead of DNA at the input")

    def transcribe(self):
        return self.seq.replace('T', 'U').replace('t', 'u')


class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence composed of A, U, G, C nucleotides.
    
    Such methods as reverse, complement, and reverse_complement are inherited from NucleicAcidSequence.
    """
    def __init__(self, seq):
        super().__init__(seq)
        if 'T' in self.seq.upper():
            raise ValueError('DNA instead of RNA at the input')
        

class AminoAcidSequence(ManageableBiologicalSequence):
    """
    Represents an amino acid sequence composed of 20 standard amino acids.
    
    Provides a method to identify the indices of start codons (Methionine, 'M') in the sequence."""
    def __init__(self, seq):
        super().__init__(seq)
        alphabet = 'ARNDCQEGHILKMFPSTWYV'
        if not self.check_alphabet(alphabet):
            raise ValueError("Invalid characters in the Amino Acid sequence.")

    def get_start_codon_idxs(self):
        idxs = []
        for i in range(len(self.seq)):
            if self.seq[i] == 'M':
                idxs.append(i)
        return idxs




def filter_fastq (input_fastq : str, output_fastq : str = 'output.fastq', gc_bounds: Union[int, float, tuple] =(0, 100), 
                  length_bounds: Union[int, tuple]=(0, 2**32), quality_threshold: int = 0):
    """
    Filter sequences from a FASTQ file by GC content, length, and quality.
    
    Reads records from `input_fastq` and writes those meeting the provided
    criteria to `output_fastq` in FASTQ format.
    
    Args:
        input_fastq (str): Path to input FASTQ file.
        output_fastq (str): Path where filtered records will be written.
        gc_bounds (Union[int, float, tuple]): GC percentage boundaries. If a
            number, lower bound is 0; if a tuple, it should be (min, max).
        length_bounds (Union[int, tuple]): Sequence length boundaries. If an
            int, lower bound is 0; if a tuple, it should be (min, max).
        quality_threshold (int): Minimum average PHRED quality score.
    """

    def parse_bounds(bounds):
        if isinstance(bounds, Number):
            lower_bound, upper_bound = 0, bounds 
        else:
            lower_bound, upper_bound = bounds
        return lower_bound, upper_bound
    
    lower_gc_bound, upper_gc_bound = parse_bounds(gc_bounds)
    lower_length_bound, upper_length_bound = parse_bounds(length_bounds)
          
    filtered_records = [] 
    with open(input_fastq) as handle:
        for seq_record in SeqIO.parse(handle, "fastq"):
            gc_content = gc_fraction(seq_record.seq)*100
            seq_length = len(seq_record)
            qualities = seq_record.letter_annotations["phred_quality"]
            avg_quality = sum(qualities) / len(qualities)
            if (lower_gc_bound <= gc_content <= upper_gc_bound and
                lower_length_bound <= seq_length <= upper_length_bound and
                quality_threshold <= avg_quality):
                filtered_records.append(seq_record)

    SeqIO.write(filtered_records, output_fastq, "fastq")

