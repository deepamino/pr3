
def check_nucleotides_DNA(sequence):
    sequence = sequence.upper()
    for nucleotide in sequence:
        if nucleotide not in "ACGT":
            return False
    return True

def check_no_spaces(sequence):
    if " " in sequence:
        sequence = sequence.replace(" ", "")
    return sequence

def sequence_upper(sequence):
    return sequence.upper()

def check_sequence(sequence):
    sequence = check_no_spaces(sequence)
    sequence = sequence_upper(sequence)
    if not check_nucleotides_DNA(sequence):
        raise ValueError("Invalid DNA sequence")
    return sequence