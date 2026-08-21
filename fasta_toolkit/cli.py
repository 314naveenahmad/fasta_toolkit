from fasta_toolkit.parser import parse_fasta
from fasta_toolkit.analysis import (
    sequence_length,
    nucleotide_count,
    gc_content,
    reverse_complement,
    transcribe,
    find_motif,
    translate,
    find_orfs,
)


def select_sequence(sequences):
    """Display available sequences and return the selected Sequence object."""

    print("\nAvailable sequences:")

    for index, seq in enumerate(sequences, start=1):
        print(f"{index}. {seq.id} - {seq.description}")

    while True:
        choice = input("\nSelect a sequence (0 to exit): ").strip()

        try:
            choice = int(choice)
        except ValueError:
            print("Please enter a valid number.")
            continue

        if choice == 0:
            return None

        if 1 <= choice <= len(sequences):
            return sequences[choice - 1]

        print(
            f"Please enter a number between 1 and {len(sequences)}, "
            "or 0 to exit."
        )


def display_menu():
    """Display the sequence analysis menu."""

    print("\nWhat would you like to do?")
    print("1. Show sequence length")
    print("2. Show nucleotide count")
    print("3. Show GC content")
    print("4. Show reverse complement")
    print("5. Transcribe DNA to RNA")
    print("6. Find motif")
    print("7. Translate sequence")
    print("8. Find ORFs")
    print("9. Exit")


def analysis_menu(sequence):
    """Display the analysis menu for a selected sequence."""

    while True:
        display_menu()

        choice = input("\nEnter your choice: ").strip()

        if choice == "1":
            length = sequence_length(sequence)
            print(f"\nSequence length: {length} bp")

        elif choice == "2":
            counts = nucleotide_count(sequence)
            print("\nNucleotide count:")
            print(f"A: {counts['A']}")
            print(f"T: {counts['T']}")
            print(f"G: {counts['G']}")
            print(f"C: {counts['C']}")

        elif choice == "3":
            gc = gc_content(sequence)
            print(f"\nGC content: {gc:.2f}%")

        elif choice == "4":
            result = reverse_complement(sequence)
            print(f"\nReverse complement:\n{result}")

        elif choice == "5":
            result = transcribe(sequence)
            print(f"\nRNA sequence:\n{result}")

        elif choice == "6":
            motif = input("\nEnter motif to search: ").strip()

            try:
                positions = find_motif(sequence, motif)

                if positions:
                    print(f"\nMotif found at positions: {positions}")
                else:
                    print("\nMotif not found.")

            except ValueError as error:
                print(f"\nError: {error}")

        elif choice == "7":
            frame_input = input(
                "\nEnter reading frame (1, 2, or 3): "
            ).strip()

            try:
                frame = int(frame_input)
                protein = translate(sequence, frame)

                print(f"\nProtein sequence (frame {frame}):")
                print(protein)

            except ValueError as error:
                print(f"\nError: {error}")

        elif choice == "8":
            orfs = find_orfs(sequence)

            if orfs:
                print("\nORFs found:")
                for index, orf in enumerate(orfs, start=1):
                    print(f"{index}. {orf}")
            else:
                print("\nNo ORFs found.")

        elif choice == "9":
            print("\nExiting FASTA Toolkit.")
            return

        else:
            print("\nInvalid choice. Please select a number from 1 to 9.")


def main():
    """Run the FASTA Toolkit CLI."""

    print("FASTA Toolkit")
    print("=============")

    filepath = input("Enter FASTA file path: ").strip()

    try:
        sequences = parse_fasta(filepath)

    except FileNotFoundError:
        print(f"Error: File not found: {filepath}")
        return

    except ValueError as error:
        print(f"Error: {error}")
        return

    print(f"\nSequences found: {len(sequences)}")

    for index, seq in enumerate(sequences, start=1):
        print(f"\n{index}. {seq.id}")
        print(f"   Description: {seq.description}")
        print(f"   Length: {sequence_length(seq)} bp")

    selected_sequence = select_sequence(sequences)

    if selected_sequence is None:
        print("\nExiting FASTA Toolkit.")
        return

    print(f"\nSelected: {selected_sequence.id}")
    print(f"Description: {selected_sequence.description}")
    print(f"Length: {sequence_length(selected_sequence)} bp")

    analysis_menu(selected_sequence)


if __name__ == "__main__":
    main()