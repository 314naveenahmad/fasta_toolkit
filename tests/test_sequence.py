import pytest

from fasta_toolkit.sequence import Sequence

def test_valid_sequence():
    seq = Sequence(
        id="seq1",
        description="Example sequence",
        sequence="ACGTACGT"
    )

    assert seq.id == "seq1"
    assert seq.description == "Example sequence"
    assert seq.sequence == "ACGTACGT"


def test_lowercase_sequence_is_converted_to_uppercase():
    seq = Sequence(
        id="seq1",
        description="Example sequence",
        sequence="acgtacgt"
    )

    assert seq.sequence == "ACGTACGT"


def test_whitespace_is_removed():
    seq = Sequence(
        id="seq1",
        description="Example sequence",
        sequence="ACGT\nAC GT"
    )

    assert seq.sequence == "ACGTACGT"


def test_empty_sequence_raises_error():
    with pytest.raises(ValueError, match="Sequence cannot be empty"):
        Sequence(
            id="seq1",
            description="Example sequence",
            sequence=""
        )


def test_empty_id_raises_error():
    with pytest.raises(ValueError, match="ID cannot be empty"):
        Sequence(
            id="",
            description="Example sequence",
            sequence="ACGTACGT"
        )


def test_invalid_bases_raise_error():
    with pytest.raises(ValueError, match="Invalid bases found"):
        Sequence(
            id="seq1",
            description="Example sequence",
            sequence="ACGTXYZ"
        )


def test_empty_description_is_allowed():
    seq = Sequence(
        id="seq1",
        description="",
        sequence="ACGTACGT"
    )

    assert seq.description == ""