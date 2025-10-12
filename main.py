from typing import Union
import os
from os import path
import modules.run_dna_rna as rdr
import modules.fastq as fq

script_dir = os.path.dirname(os.path.abspath(__file__))

def run_dna_rna_tools(*args):
    process = args[-1]          
    sequences = args[:-1]       # get any number of sequences
    checking = rdr.check_sequnces(sequences)[0]

    if process == "is_nucleic_acid":
        result = checking
    elif False in checking:
        result = f"Invalid sequences: {rdr.check_sequnces(sequences)[1]}"
    else:
        result = []
        # applying the tool to each sequence
        for seq in sequences:
            if "U" in seq.upper():
                is_rna = True
            else:
                is_rna = False
            # Procedure selection
            if process == "transcribe":
                inter_res = rdr.transcribe(seq, is_rna)
            elif process == "reverse":
                inter_res = rdr.reverse(seq)
            elif process == "complement":
                inter_res = rdr.complement(seq, is_rna)
            elif process == "reverse_complement":
                inter_res = rdr.reverse(seq)
                inter_res = rdr.complement(inter_res, is_rna)
            result.append(inter_res)

    if len(result) == 1:
        return result[0]
    else:
        return result


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

