(sequence.count('G') + sequence.count('C')) / len(sequence) * 100
def is_relevant_gc(sequence: str, gc_bounds) -> bool:
    """
    Computes GC-content of sequence

    Arguments:
    sequence : str
    gc_bounds: int / float / tuple with 2 elements

    Return True if GC-content of sequence in gc_bounds
    """

    if isinstance(gc_bounds, (int, float)):
        min_gc = 0
        max_gc = gc_bounds
    else:
        min_gc = gc_bounds[0]
        max_gc = gc_bounds[1]

    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    if min_gc <= gc_content <= max_gc:
        return True
    else:
        return False


def is_relevant_length(sequence: str, length_bounds) -> bool:
    """
    Computes length of sequence

    Arguments:
    sequence : str
    length_bounds: int / tuple with 2 elements

    Return True if length of sequence in length_bounds
    """

    if isinstance(length_bounds, (int)):
        min_length = 0
        max_length = length_bounds
    else:
        min_length = length_bounds[0]
        max_length = length_bounds[1]

    return  min_length <= len(sequence) <= max_length


def is_relevant_quality(phred: str, quality_threshold) -> bool:
    """
    Computes length of sequence

    Arguments:
    phred : str
    quality_threshold: int

    Return True if average quality of the read (from phred)
    is greater than the quality_threshold
    """

    score = 0
    for i in phred:
        score += ord(i) - 33
    if score/len(phred) >= quality_threshold:
        return True
    else:
        return False

def write_relevant_fastq(output_fastq, processed_sequence):
    if path.isfile(output_fastq):
                with open(output_fastq, 'a') as output_file:
                    for row in processed_sequence:
                        output_file.write(row)
    else:
        with open(output_fastq, 'w') as output_file:
            for row in processed_sequence:
                output_file.write(row)
