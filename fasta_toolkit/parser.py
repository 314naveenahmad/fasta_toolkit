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

    