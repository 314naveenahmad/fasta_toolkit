from fasta_toolkit.analysis import sequence_length
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