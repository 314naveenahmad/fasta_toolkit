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


def transcribe(sequence: Sequence) -> str:
    """
    Transcribe a DNA sequence to RNA by replacing thymine (T) with uracil (U).
    
    Args:
        sequence (Sequence): The DNA sequence object.
        
    Returns:
        str: The transcribed RNA sequence.
        
    """

    return sequence.sequence.replace("T", "U")


def find_motif(sequence: Sequence, motif: str) -> list[int]:
    """
    Find all occurrences of a motif in a DNA sequnce.

    Args:
        sequence (Sequence): The DNA sequence object.
        motif (str): The motif to search for.

    Returns:
        list: A list of starting positions (1-based) where motif occurs.

    """

    motif = motif.upper()
    if not motif:
        raise ValueError("Motif cannot be empty.")

    positions = []
    start = 0
    while True:
        start = sequence.sequence.find(motif, start)
        if start == -1:
            break
        positions.append(start + 1)
        start += 1  # Move past the last found motif

    return positions

