import pytest

from fasta_toolkit.parser import parse_header


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