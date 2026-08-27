# FASTA Toolkit

A Python command-line toolkit for parsing and analyzing FASTA sequences.

FASTA Toolkit provides a collection of common DNA sequence analysis functions through
a simple interactive command-line interface. It is designed as a learning-focused
bioinformatics project while following software engineering practices such as modular
code, type hints, documentation, and automated testing.

---

## Features

### FASTA Parsing
- Parse FASTA files containing one or multiple sequences
- Extract sequence IDs and descriptions
- Handle multiline FASTA sequences
- Normalize sequences to uppercase
- Ignore blank lines
- Validate FASTA input

### Sequence Analysis
- Calculate sequence length
- Count nucleotides (A, T, G, C)
- Calculate GC content
- Generate reverse complements
- Transcribe DNA into RNA
- Search for DNA motifs, including overlapping occurrences
- Translate DNA into protein sequences
- Support reading frames 1, 2, and 3
- Detect open reading frames (ORFs)
- Detect ORFs on both forward and reverse-complement strands

### Command-Line Interface

The toolkit provides an interactive CLI where users can:

1. Load a FASTA file
2. Select a sequence
3. Perform sequence analysis
4. Search for motifs
5. Translate DNA
6. Detect ORFs
7. Exit the program

---

## Project Structure

```text
fasta_toolkit/
│
├── fasta_toolkit/
│   ├── __init__.py
│   ├── sequence.py
│   ├── parser.py
│   ├── analysis.py
│   └── cli.py
│
├── tests/
│   ├── test_sequence.py
│   ├── test_parser.py
│   ├── test_analysis.py
│   └── test_cli.py
│
├── examples/
│   └── sample.fasta
│
├── README.md
├── requirements.txt
└── .gitignore

## Installation

1- Clone the repository 
git clone https://github.com/314naveenahmad/fasta_toolkit.git

2- Navigate into the project
cd fasta_toolkit

3- Create a virtual environment
python -m venv venv

4- Activate the virtual environment
venv\Scripts\Activate.ps1

5- Install the dependencies
pip install -r requirements.txt

## Usage
Run the toolkit using
python -m fasta_toolkit.cli

The program will ask for the path to a FASTA file

FASTA Toolkit
=============

Enter FASTA file path:

After loading the file, available sequences are displayed

Available sequences:
1. MH011443.1 - Homo sapiens TP53 gene
2. seq1 - Human TP53 fragment
3. seq2 - Example sequence

Select a sequence (0 to exit):

After selecting a sequence, the analysis menu becomes available

What would you like to do?
1. Show sequence length
2. Show nucleotide count
3. Show GC content
4. Show reverse complement
5. Transcribe DNA to RNA
6. Find motif
7. Translate sequence
8. Find ORFs
9. Exit

## Testing

Run the complete test suite with:

pytest -v

Current test suite:

80 passed

The tests cover:

Sequence validation
FASTA parsing
Sequence length
Nucleotide counting
GC content
Reverse complement
DNA transcription
Motif searching
DNA translation
ORF detection
CLI functionality
Error handling
Technologies Used
Python 3
pytest
Git
GitHub

## Project Goals

This project was developed to combine bioinformatics concepts with Python programming and software engineering practices.

The main goals are:

Build practical Python programming skills
Implement fundamental DNA sequence analysis algorithms
Learn how bioinformatics tools process FASTA data
Practice modular software design
Write automated tests for scientific code
Build a usable command-line bioinformatics tool
Future Improvements

## Possible future versions may include:

Support for additional FASTA formats
Command-line arguments using argparse
Protein sequence analysis
More advanced ORF analysis
Export analysis results to files
Batch analysis of multiple FASTA files
Performance improvements for very large sequences
Web-based interface


## Repository

GitHub:

https://github.com/314naveenahmad/fasta_toolkit

## Author

Naveen Ahmad

MSc Biotechnology | Bioinformatics & Computational Biology