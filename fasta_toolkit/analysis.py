from .sequence import Sequence

def sequence_length(sequence: Sequence) -> int:
    """
    Returns the length of a given sequence.

    Args:
        sequence (Sequence): The sequence object.

    Returns:
        int: The length of the sequence.
    """
    return len(sequence.sequence)


def nucleotide_count(sequence: Sequence) -> dict[str, int]:
    """
    Returns a dictionary with the count of each nucleotide in the sequence.

    Args:
        sequence (Sequence): The sequence object.

    Returns:
        dict: A dictionary with nucleotides as keys and their counts as values.
    """
    counts = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
    for nucleotide in sequence.sequence:
        if nucleotide in counts:
            counts[nucleotide] += 1
    return counts

def gc_content(sequence: Sequence) -> float:
    """
    Calculates the GC content of a given sequence.

    Args:
        sequence (Sequence): The sequence object.

    Returns:
        float: The GC content as a percentage.
    """
    counts = nucleotide_count(sequence)
    total_count = sequence_length(sequence)
    if total_count == 0:
        return 0.0
    gc_count = counts["G"] + counts["C"]
    return round((gc_count / total_count) * 100, 2)


def reverse_complement(sequence: Sequence) -> str:
    """
    Return the reverse complement of a given sequence.
    
    Args:
        sequence (Sequence): The sequence object
    
    Returns:
        str: The reverse complement of the sequence.
    
    """
    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }

    return "".join(complement[base] for base in reversed(sequence.sequence))