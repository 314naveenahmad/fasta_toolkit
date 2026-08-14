from dataclasses import dataclass


@dataclass
class Sequence:
    id: str
    description: str
    sequence: str

    def __post_init__(self):
        self.sequence = "".join(self.sequence.split()).upper()

        if not self.sequence:
            raise ValueError("Sequence cannot be empty")

        if not self.id:
            raise ValueError("ID cannot be empty")


        valid_bases = set("ACTGN")
        invalid_bases = set(self.sequence) - valid_bases
        if invalid_bases:
            raise ValueError(f"Invalid bases found in sequence: {','.join(sorted(invalid_bases))}")
        