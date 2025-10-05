# Checking the validity of sequences
def check_sequnces(sequences):
    alphabet = {"A", "T", "G", "C", "U", "a", "t", "g", "c", "u"}
    checking = []
    invalid_sequence = []
    for seq in sequences:
        unique_chars = set(seq.upper())
        if {"U", "T"}.issubset(unique_chars):
            check = False
        else:
            check = unique_chars.issubset(alphabet) & bool(unique_chars)
        if check is False:
            invalid_sequence.append(seq)
        checking.append(check)
    return checking, invalid_sequence


def transcribe(seq, is_rna):
    if is_rna is False:
        seq = seq.replace("T", "U")
        seq = seq.replace("t", "u")
    return seq


def reverse(seq):
    return seq[::-1]


def complement(seq, is_rna):
    if is_rna is True:
        complement_pairs = str.maketrans("AUGCaugc", "UACGuacg")
    else:
        complement_pairs = str.maketrans("ATGCatgc", "TACGtacg")
    result = seq.translate(complement_pairs)
    return result
