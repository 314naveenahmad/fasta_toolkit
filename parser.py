from .sequence import Sequence

def parse_header(header_line: str) -> tuple[str, str]:
    """
    Parses the header line of a FASTA file and returns the sequence ID and description.

    Args:
        header_line (str): The header line from a FASTA file, starting with '>'.

    Returns:
        tuple: A tuple containing the sequence ID (str) and description (str).
    """
    header_line = header_line.strip()

    if not header_line.startswith('>'):
        raise ValueError("Invalid FASTA header line")

    # Remove the '>' character and split by whitespace
    header_content = header_line[1:].strip()

    if not header_content:
        raise ValueError("FASTA header cannot be empty")
    
    parts = header_content.split(None, 1)  # Split into at most two parts

    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], parts[1]


def parse_fasta(filepath: str) -> list[Sequence]:
    """
    Parses a FASTA file and returns a list of Sequence objects.

    Args:
        filepath (str): The path to the FASTA file.

        Returns:
        A list of Sequence objects parsed from the file.

    Raises:
        ValueError: If sequence data occurs before a header or
            if the FASTA file contains no sequences.
        FileNotFoundError: If the FASTA file does not exist.
        
        """

    sequences = []
    with open(filepath, "r", encoding="utf-8") as file:
        current_id = None
        current_description = None
        current_sequence_lines = []

        for line in file:
            line = line.strip()

            if not line:
                continue  # Skip empty lines

            if line.startswith('>'):
                # If we have a current sequence, save it before starting a new one
                if current_id is not None:
                    sequence_str = ''.join(current_sequence_lines)
                    sequences.append(Sequence(id=current_id, description=current_description, sequence=sequence_str))
                    current_sequence_lines = []

                # Parse the new header
                current_id, current_description = parse_header(line)
            else:
                if current_id is None:
                    raise ValueError("FASTA file format error: sequence data found before any header line")
                # Accumulate sequence lines
                current_sequence_lines.append(line)

        # Don't forget to add the last sequence after the loop ends
        if current_id is not None:
            sequence_str = ''.join(current_sequence_lines)
            sequences.append(Sequence(id=current_id, description=current_description, sequence=sequence_str))

    if not sequences:
        raise ValueError("FASTA file contains no sequences")
    
    return sequences