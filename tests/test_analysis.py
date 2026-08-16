from fasta_toolkit.analysis import sequence_length, nucleotide_count, gc_content, reverse_complement
from fasta_toolkit.sequence import Sequence


def test_sequence_length():
    seq = Sequence(
        id="seq1",
        description="Example sequence",
        sequence="ATGCGTAC"
    )

    assert sequence_length(seq) == 8


def test_sequence_length_single_base():
    seq = Sequence(
        id="seq1",
        description="Single base",
        sequence="A"
    )

    assert sequence_length(seq) == 1


def test_sequence_length_with_ambiguous_bases():
    seq = Sequence(
        id="seq1",
        description="Sequence containing N",
        sequence="ATGCNN"
    )

    assert sequence_length(seq) == 6



def test_sequence_length():
    seq = Sequence(
        id="seq1",
        description="Example sequence",
        sequence="ATGCGTAC"
    )

    assert sequence_length(seq) == 8


def test_sequence_length_single_base():
    seq = Sequence(
        id="seq1",
        description="Single base",
        sequence="A"
    )

    assert sequence_length(seq) == 1


def test_sequence_length_with_ambiguous_bases():
    seq = Sequence(
        id="seq1",
        description="Sequence containing N",
        sequence="ATGCNN"
    )

    assert sequence_length(seq) == 6


def test_nucleotide_count():
    seq = Sequence(
        id="seq1",
        description="Normal sequence",
        sequence="ATGCGT"
    )

    assert nucleotide_count(seq) == {
        "A": 1,
        "T": 2,
        "G": 2,
        "C": 1,
        "N": 0
    }


def test_nucleotide_count_with_ambiguous_bases():
    seq = Sequence(
        id="seq1",
        description="Sequence containing N",
        sequence="ATGCNN"
    )

    assert nucleotide_count(seq) == {
        "A": 1,
        "T": 1,
        "G": 1,
        "C": 1,
        "N": 2
    }


def test_nucleotide_count_single_nucleotide():
    seq = Sequence(
        id="seq1",
        description="Only adenine",
        sequence="AAAA"
    )

    assert nucleotide_count(seq) == {
        "A": 4,
        "T": 0,
        "G": 0,
        "C": 0,
        "N": 0
    }


def test_nucleotide_count_lowercase_sequence():
    seq = Sequence(
        id="seq1",
        description="Lowercase sequence",
        sequence="atgcNN"
    )

    assert nucleotide_count(seq) == {
        "A": 1,
        "T": 1,
        "G": 1,
        "C": 1,
        "N": 2
    }

def test_gc_content():
    seq = Sequence(
        id="seq1",
        description="Normal sequence",
        sequence="ATGCGT"
    )

    assert gc_content(seq) == 50.0


def test_gc_content_100_percent():
    seq = Sequence(
        id="seq1",
        description="100 percent GC",
        sequence="GGCC"
    )

    assert gc_content(seq) == 100.0


def test_gc_content_0_percent():
    seq = Sequence(
        id="seq1",
        description="0 percent GC",
        sequence="AATT"
    )

    assert gc_content(seq) == 0.0


def test_gc_content_with_ambiguous_bases():
    seq = Sequence(
        id="seq1",
        description="Sequence containing N",
        sequence="ATGCNN"
    )

    assert gc_content(seq) == 33.33


def test_reverse_complement():
    seq = Sequence(
        id="seq1",
        description="Example sequence",
        sequence="ATGC"
    )

    assert reverse_complement(seq) == "GCAT"


def test_reverse_complement_palindrome():
    seq = Sequence(
        id="seq1",
        description="Palindromic sequence",
        sequence="ATGCAT"
    )

    assert reverse_complement(seq) == "ATGCAT"


def test_reverse_complement_single_base():
    seq = Sequence(
        id="seq1",
        description="Single nucleotide",
        sequence="A"
    )

    assert reverse_complement(seq) == "T"


def test_reverse_complement_with_ambiguous_base():
    seq = Sequence(
        id="seq1",
        description="Sequence containing N",
        sequence="ATGCN"
    )

    assert reverse_complement(seq) == "NGCAT"


def test_reverse_complement_lowercase_sequence():
    seq = Sequence(
        id="seq1",
        description="Lowercase sequence",
        sequence="atgc"
    )

    assert reverse_complement(seq) == "GCAT"


