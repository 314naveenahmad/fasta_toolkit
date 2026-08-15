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
