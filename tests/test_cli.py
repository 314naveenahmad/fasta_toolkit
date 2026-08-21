import pytest

from fasta_toolkit.cli import select_sequence, display_menu, analysis_menu
from fasta_toolkit.sequence import Sequence


def test_select_sequence_valid_choice(monkeypatch):
    sequences = [
        Sequence(
            id="seq1",
            description="First sequence",
            sequence="ATGC"
        ),
        Sequence(
            id="seq2",
            description="Second sequence",
            sequence="GGCC"
        ),
    ]

    monkeypatch.setattr("builtins.input", lambda _: "2")

    selected = select_sequence(sequences)

    assert selected == sequences[1]


def test_select_sequence_exit(monkeypatch):
    sequences = [
        Sequence(
            id="seq1",
            description="First sequence",
            sequence="ATGC"
        )
    ]

    monkeypatch.setattr("builtins.input", lambda _: "0")

    selected = select_sequence(sequences)

    assert selected is None


def test_select_sequence_invalid_input(monkeypatch, capsys):
    sequences = [
        Sequence(
            id="seq1",
            description="First sequence",
            sequence="ATGC"
        )
    ]

    inputs = iter(["abc", "1"])
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))

    selected = select_sequence(sequences)

    captured = capsys.readouterr()

    assert "Please enter a valid number." in captured.out
    assert selected == sequences[0]


def test_select_sequence_out_of_range(monkeypatch, capsys):
    sequences = [
        Sequence(
            id="seq1",
            description="First sequence",
            sequence="ATGC"
        )
    ]

    inputs = iter(["5", "1"])
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))

    selected = select_sequence(sequences)

    captured = capsys.readouterr()

    assert "Please enter a number between 1 and 1, or 0 to exit." in captured.out
    assert selected == sequences[0]


def test_display_menu(capsys):
    display_menu()

    captured = capsys.readouterr()

    assert "1. Show sequence length" in captured.out
    assert "2. Show nucleotide count" in captured.out
    assert "3. Show GC content" in captured.out
    assert "4. Show reverse complement" in captured.out
    assert "5. Transcribe DNA to RNA" in captured.out
    assert "6. Find motif" in captured.out
    assert "7. Translate sequence" in captured.out
    assert "8. Find ORFs" in captured.out
    assert "9. Exit" in captured.out


def test_analysis_menu_sequence_length(monkeypatch, capsys):
    sequence = Sequence(
        id="seq1",
        description="Test sequence",
        sequence="ATGC"
    )

    inputs = iter(["1", "9"])
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))

    analysis_menu(sequence)

    captured = capsys.readouterr()

    assert "Sequence length: 4 bp" in captured.out
    assert "Exiting FASTA Toolkit." in captured.out


def test_analysis_menu_nucleotide_count(monkeypatch, capsys):
    sequence = Sequence(
        id="seq1",
        description="Test sequence",
        sequence="ATGC"
    )

    inputs = iter(["2", "9"])
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))

    analysis_menu(sequence)

    captured = capsys.readouterr()

    assert "A: 1" in captured.out
    assert "T: 1" in captured.out
    assert "G: 1" in captured.out
    assert "C: 1" in captured.out


def test_analysis_menu_gc_content(monkeypatch, capsys):
    sequence = Sequence(
        id="seq1",
        description="Test sequence",
        sequence="ATGC"
    )

    inputs = iter(["3", "9"])
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))

    analysis_menu(sequence)

    captured = capsys.readouterr()

    assert "GC content: 50.00%" in captured.out


def test_analysis_menu_invalid_choice(monkeypatch, capsys):
    sequence = Sequence(
        id="seq1",
        description="Test sequence",
        sequence="ATGC"
    )

    inputs = iter(["99", "9"])
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))

    analysis_menu(sequence)

    captured = capsys.readouterr()

    assert "Invalid choice. Please select a number from 1 to 9." in captured.out


def test_analysis_menu_motif(monkeypatch, capsys):
    sequence = Sequence(
        id="seq1",
        description="Test sequence",
        sequence="ATGCGATG"
    )

    inputs = iter(["6", "ATG", "9"])
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))

    analysis_menu(sequence)

    captured = capsys.readouterr()

    assert "Motif found at positions: [1, 6]" in captured.out


def test_analysis_menu_translation(monkeypatch, capsys):
    sequence = Sequence(
        id="seq1",
        description="Test sequence",
        sequence="ATGGCC"
    )

    inputs = iter(["7", "1", "9"])
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))

    analysis_menu(sequence)

    captured = capsys.readouterr()

    assert "Protein sequence (frame 1):" in captured.out
    assert "MA" in captured.out


def test_analysis_menu_orfs(monkeypatch, capsys):
    sequence = Sequence(
        id="seq1",
        description="Test sequence",
        sequence="ATGAAATAG"
    )

    inputs = iter(["8", "9"])
    monkeypatch.setattr("builtins.input", lambda _: next(inputs))

    analysis_menu(sequence)

    captured = capsys.readouterr()

    assert "ORFs found:" in captured.out
    assert "ATGAAATAG" in captured.out