import pytest

from fasta_toolkit.parser import parse_header, parse_fasta


def test_parse_header_with_id_and_description():
    header = ">seq1 Human TP53 gene fragment"

    result = parse_header(header)

    assert result == ("seq1", "Human TP53 gene fragment")


def test_parse_header_with_id_only():
    header = ">seq1"

    result = parse_header(header)

    assert result == ("seq1", "")


def test_parse_header_with_extra_whitespace():
    header = "   >seq1 Human TP53 gene fragment   "

    result = parse_header(header)

    assert result == ("seq1", "Human TP53 gene fragment")


def test_parse_header_without_greater_than_symbol():
    header = "seq1 Human TP53 gene fragment"

    with pytest.raises(ValueError, match="Invalid FASTA header line"):
        parse_header(header)


def test_parse_empty_header():
    header = ">"

    with pytest.raises(ValueError, match="FASTA header cannot be empty"):
        parse_header(header)


def test_parse_header_with_multiple_spaces():
    header = ">seq1    Human    TP53    gene fragment"

    result = parse_header(header)

    assert result == ("seq1", "Human    TP53    gene fragment")


def test_parse_single_fasta_record(tmp_path):
    fasta_file = tmp_path / "single.fasta"

    fasta_file.write_text(
        ">seq1 Example sequence\n"
        "ATGCATGC\n"
    )

    sequences = parse_fasta(str(fasta_file))

    assert len(sequences) == 1
    assert sequences[0].id == "seq1"
    assert sequences[0].description == "Example sequence"
    assert sequences[0].sequence == "ATGCATGC"


def test_parse_multiple_fasta_records(tmp_path):
    fasta_file = tmp_path / "multiple.fasta"

    fasta_file.write_text(
        ">seq1 First sequence\n"
        "ATGCATGC\n"
        ">seq2 Second sequence\n"
        "GCTAGCTA\n"
    )

    sequences = parse_fasta(str(fasta_file))

    assert len(sequences) == 2

    assert sequences[0].id == "seq1"
    assert sequences[0].description == "First sequence"
    assert sequences[0].sequence == "ATGCATGC"

    assert sequences[1].id == "seq2"
    assert sequences[1].description == "Second sequence"
    assert sequences[1].sequence == "GCTAGCTA"


def test_parse_multiline_sequence(tmp_path):
    fasta_file = tmp_path / "multiline.fasta"

    fasta_file.write_text(
        ">seq1 Multiline sequence\n"
        "ATGC\n"
        "GCTA\n"
        "CGAT\n"
    )

    sequences = parse_fasta(str(fasta_file))

    assert len(sequences) == 1
    assert sequences[0].sequence == "ATGCGCTACGAT"


def test_parse_fasta_with_blank_lines(tmp_path):
    fasta_file = tmp_path / "blank_lines.fasta"

    fasta_file.write_text(
        ">seq1 First sequence\n"
        "ATGC\n"
        "\n"
        "GCTA\n"
        "\n"
        ">seq2 Second sequence\n"
        "CGAT\n"
    )

    sequences = parse_fasta(str(fasta_file))

    assert len(sequences) == 2
    assert sequences[0].sequence == "ATGCGCTA"
    assert sequences[1].sequence == "CGAT"


def test_parse_fasta_normalizes_sequence(tmp_path):
    fasta_file = tmp_path / "normalize.fasta"

    fasta_file.write_text(
        ">seq1 Example sequence\n"
        "atgc\n"
        "AT GC\n"
    )

    sequences = parse_fasta(str(fasta_file))

    assert sequences[0].sequence == "ATGCATGC"


def test_parse_fasta_without_description(tmp_path):
    fasta_file = tmp_path / "no_description.fasta"

    fasta_file.write_text(
        ">seq1\n"
        "ATGCATGC\n"
    )

    sequences = parse_fasta(str(fasta_file))

    assert len(sequences) == 1
    assert sequences[0].id == "seq1"
    assert sequences[0].description == ""
    assert sequences[0].sequence == "ATGCATGC"


def test_sequence_before_header_raises_error(tmp_path):
    fasta_file = tmp_path / "invalid.fasta"

    fasta_file.write_text(
        "ATGCATGC\n"
        ">seq1 Example sequence\n"
        "GCTA\n"
    )

    with pytest.raises(
        ValueError,
        match="sequence data found before any header line"
    ):
        parse_fasta(str(fasta_file))


def test_empty_fasta_file_raises_error(tmp_path):
    fasta_file = tmp_path / "empty.fasta"
    fasta_file.write_text("")

    with pytest.raises(
        ValueError,
        match="FASTA file contains no sequences"
    ):
        parse_fasta(str(fasta_file))