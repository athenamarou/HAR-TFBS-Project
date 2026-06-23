def parse_meme_file(jaspar_input):
    """Meme file parser because Bio motif could not parse this version of jaspar meme file"""
    motif_list = []
    with open(jaspar_input, "r") as f:
        lines = f.read().strip().split("\n")
        i = 0
        background_letter_freqs = {
            "A": 0.25,
            "C": 0.25,
            "G": 0.25,
            "T": 0.25,
        }  # default background frequencies given by the file

        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("Background letter frequencies"):
                i += 1
                parts = lines[i].strip().split()

                # catch background freqs directly from file if present
                background_letter_freqs = {
                    parts[j]: float(parts[j + 1]) for j in range(0, len(parts), 2)
                }

            # get matrix_id and name
            elif line.startswith("MOTIF"):
                parts = line.split()
                matrix_id = parts[1]
                name = parts[2] if len(parts) > 2 else matrix_id

                # letter probability matrix header
                i += 1
                while i < len(lines) and "letter-probability matrix" not in lines[i]:
                    i += 1
                if i >= len(lines):
                    break

                # parse w= and nsites= from the header line
                header = lines[i]
                w = int(header.split("w=")[1].split()[0])
                nsites = (
                    int(header.split("nsites=")[1].split()[0])
                    if "nsites=" in header
                    else 0
                )

                # read rows of probabilities
                matrix = []
                i += 1  # moving to the next line
                while i < len(lines) and len(matrix) < w:
                    row = lines[i].strip()
                    if row and not row.startswith(("MOTIF", "URL", "//")):
                        vals = [float(x) for x in row.split()]
                        if (
                            len(vals) == 4
                        ):  # each row in the matrix represents one position in the motif with probabilities for 4 nucleotides
                            matrix.append(vals)
                    i += 1

                # Store ACGT probability lists
                motif_list.append(
                    {
                        "matrix_id": matrix_id,
                        "name": name,
                        "width": w,
                        "nsites": nsites,
                        "probs": {
                            "A": [row[0] for row in matrix],
                            "C": [row[1] for row in matrix],
                            "G": [row[2] for row in matrix],
                            "T": [row[3] for row in matrix],
                        },
                    }
                )
                continue
            i += 1
        return motif_list, background_letter_freqs


def calculate_background_frequencies(motif_list):
    #  Initialize counters for each nucleotide
    totals = {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}
    total_positions = 0

    for motif in motif_list:
        probs = motif["probs"]
        # The number of positions is the length of any probability list (e.g., 'A')
        num_positions = len(probs["A"])
        total_positions += num_positions

        #  Sum up probabilities for each letter
        for letter in totals:
            totals[letter] += sum(probs[letter])

    #  Divide by total_positions to get the average (frequency)
    # This ensures A + C + G + T = 1.0
    bg_frequencies = {
        letter: count / total_positions for letter, count in totals.items()
    }

    return bg_frequencies


if __name__ == "__main__":
    from pathlib import Path
    motif_list, _ = parse_meme_file(
        str(Path(__file__).resolve().parent.parent / "data" / "meme_file.txt")
    )
    bg = calculate_background_frequencies(motif_list)
    print("Background frequencies:")
    for letter, freq in bg.items():
        print(f"  {letter}: {freq:.4f}")
