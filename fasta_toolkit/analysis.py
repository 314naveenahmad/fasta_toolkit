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


CODON_TABLE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y",
    "TAA": "*", "TAG": "*", "TGA": "*",
    "TGT": "C", "TGC": "C",
    "TGG": "W",

    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(sequence: Sequence, frame: int = 1) -> str:
    """
    Translate a DNA sequence into a protein sequence.

    Args:
        sequence: The DNA Sequence object.
        frame: Reading frame (1, 2, or 3).

    Returns:
        The translated protein sequence.

    """

    if frame not in [1, 2, 3]:
        raise ValueError("Frame must be 1, 2, or 3")

    dna = sequence.sequence[frame - 1:] # Adjust for frame
    protein = []

    for i in range(0, len(dna) -2, 3):
        codon = dna[i:i + 3]
        protein.append(CODON_TABLE.get(codon, "X"))

    return "".join(protein)


def find_orfs(sequence: Sequence) -> list[str]:
    """
    Find open reading frames in both strands of a DNA sequence.

    An ORF starts with ATG and ends at the first in-frame
    stop codon (TAA, TAG, or TGA).

    All six reading frames are searched:
    three on the forward strand and three on the reverse
    complement.

    Args:
        sequence: The DNA Sequence object.

    Returns:
        A list of ORF DNA sequences.
    """

    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    strands = [
        sequence.sequence,
        reverse_complement(sequence)
    ]

    for dna in strands:
        for frame in range(3):
            start = frame

            while start <= len(dna) - 3:
                codon = dna[start:start + 3]

                if codon == "ATG":
                    for i in range(start + 3, len(dna) - 2, 3):
                        stop = dna[i:i + 3]

                        if stop in stop_codons:
                            orfs.append(dna[start:i + 3])
                            break

                start += 3

    return orfs


