import pytest
from fasta_toolkit.analysis import sequence_length, nucleotide_count, gc_content, reverse_complement, transcribe, find_motif, translate
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


def test_transcribe():
    seq = Sequence(
        id="seq1",
        description="Example DNA",
        sequence="ATGC"
    )

    assert transcribe(seq) == "AUGC"


def test_transcribe_all_thymine():
    seq = Sequence(
        id="seq1",
        description="Thymine sequence",
        sequence="TTTT"
    )

    assert transcribe(seq) == "UUUU"


def test_transcribe_without_thymine():
    seq = Sequence(
        id="seq1",
        description="No thymine",
        sequence="AGCG"
    )

    assert transcribe(seq) == "AGCG"


def test_transcribe_with_ambiguous_base():
    seq = Sequence(
        id="seq1",
        description="Sequence containing N",
        sequence="ATGCN"
    )

    assert transcribe(seq) == "AUGCN"


def test_transcribe_lowercase_sequence():
    seq = Sequence(
        id="seq1",
        description="Lowercase DNA",
        sequence="atgc"
    )

    assert transcribe(seq) == "AUGC"


def test_find_motif():
    seq = Sequence(
        id="seq1",
        description="Example sequence",
        sequence="ATGCGTAC"
    )

    assert find_motif(seq, "CGT") == [4]


def test_find_motif_multiple_occurrences():
    seq = Sequence(
        id="seq1",
        description="Repeated motif",
        sequence="ATGCGATGACCTGATG"
    )

    assert find_motif(seq, "ATG") == [1, 6, 14]


def test_find_motif_not_found():
    seq = Sequence(
        id="seq1",
        description="No motif",
        sequence="AAAAAA"
    )

    assert find_motif(seq, "CG") == []


def test_find_motif_overlapping():
    seq = Sequence(
        id="seq1",
        description="Overlapping motif",
        sequence="ATATAT"
    )

    assert find_motif(seq, "ATA") == [1, 3]


def test_find_motif_lowercase_motif():
    seq = Sequence(
        id="seq1",
        description="DNA sequence",
        sequence="ATGCGTATG"
    )

    assert find_motif(seq, "atg") == [1, 7]


def test_find_motif_empty_motif():
    seq = Sequence(
        id="seq1",
        description="Example sequence",
        sequence="ATGC"
    )

    with pytest.raises(ValueError, match="Motif cannot be empty"):
        find_motif(seq, "")


def test_translate():
    seq = Sequence(
        id="seq1",
        description="Basic coding sequence",
        sequence="ATGGCCATT"
    )

    assert translate(seq) == "MAI"


def test_translate_with_stop_codon():
    seq = Sequence(
        id="seq1",
        description="Sequence with stop",
        sequence="ATGTAA"
    )

    assert translate(seq) == "M*"


def test_translate_frame_2():
    seq = Sequence(
        id="seq1",
        description="Frame 2",
        sequence="ATGGCCATT"
    )

    assert translate(seq, frame=2) == "WP"


def test_translate_frame_3():
    seq = Sequence(
        id="seq1",
        description="Frame 3",
        sequence="ATGGCCATT"
    )

    assert translate(seq, frame=3) == "GH"


def test_translate_incomplete_codon():
    seq = Sequence(
        id="seq1",
        description="Incomplete codon",
        sequence="ATGGCCA"
    )

    assert translate(seq) == "MA"


def test_translate_ambiguous_codon():
    seq = Sequence(
        id="seq1",
        description="Ambiguous codon",
        sequence="ATGGCN"
    )

    assert translate(seq) == "MX"


def test_translate_invalid_frame():
    seq = Sequence(
        id="seq1",
        description="Invalid frame",
        sequence="ATGGCC"
    )

    with pytest.raises(ValueError, match="Frame must be 1, 2, or 3"):
        translate(seq, frame=4)



