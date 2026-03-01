from typing import Union
import os
from os import path
from abc import ABC, abstractmethod


script_dir = os.path.dirname(os.path.abspath(__file__))


class BiologicalSequence(ABC):

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
    def __init__(self, seq):
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
        if not self.check_alphabet(alphabet):
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
    def __init__(self, seq):
        super().__init__(seq)
        if 'U' in self.seq.upper():
            raise ValueError("RNA instead of DNA at the input")

    def transcribe(self):
        return self.seq.replace('T', 'U').replace('t', 'u')


class RNASequence(NucleicAcidSequence):
    def __init__(self, seq):
        super().__init__(seq)
        if 'T' in self.seq.upper():
            raise ValueError('DNA instead of RNA at the input')
        

class AminoAcidSequence(ManageableBiologicalSequence):
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




def filter_fastq (input_fastq, output_fastq = 'output_file.fastq', gc_bounds: Union[int, float, tuple] =(0, 100), 
                  length_bounds: Union[int, tuple]=(0, 2**32), quality_threshold: int = 0) -> dict:

    input_fastq = os.path.join(script_dir, input_fastq)
    output_fastq = os.path.join(script_dir, output_fastq)

    with open(input_fastq) as fastq_file:
        while True:
            name = fastq_file.readline()
            if not name:  # the end of file
                break
            sequence = fastq_file.readline()
            comment = fastq_file.readline()
            phred = fastq_file.readline()

    	    # check each fastq-sequence for compliance with the specified conditions
            gc_result = fq.is_relevant_gc(sequence, gc_bounds)
            length_result = fq.is_relevant_length(sequence, length_bounds)
            quality_result = fq.is_relevant_quality(phred, quality_threshold)
	    
	    # compile a file with relevant sequences
            if gc_result and length_result and quality_result:
                processed_sequence = []
                processed_sequence.extend((name, sequence, comment, phred))
                fq.write_relevant_fastq(output_fastq, processed_sequence)

