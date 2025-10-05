from typing import Union
import modules.run_dna_rna as rdr
import modules.fastq as fq


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


def filter_fastq (seqs: dict[str, tuple], gc_bounds: Union[int, float, tuple] =(0, 100), 
                  length_bounds: Union[int, tuple]=(0, 2**32), quality_threshold: int = 0) -> dict:
    """
    Seqs is dictionary of fastq sequences
        key : str - name of the sequence
        value : tuple of two strings - sequence and quality
    """

    filtered_fastq = dict()
    # check each fastq-sequence for compliance with the specified conditions
    for name, read in seqs.items():
        sequence = read[0].upper() 
        quality = read[1].upper()
        gc_result = fq.is_relevant_gc(sequence, gc_bounds)
        length_result = fq.is_relevant_length(sequence, length_bounds)
        quality_result = fq.is_relevant_quality(phred, quality_threshold)
        if gc_result & length_result & quality_result: 
        # fastq-sequences are included in the dictionary only if they satisfy all the conditions
            filtered_fastq[name] = (sequence, quality)
        else:
            continue
    return(filtered_fastq)
