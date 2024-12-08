"""
Utility functions for the project
"""


def convert_to_fasta(input_file, output_file):
	"""
	Convert gdoc text to FASTA

	Converts the google docs I was given into more usable fasta
	"""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        label = None
        sequence = []
        in_block = False  # Flag to indicate we are in a block

        for line in infile:
            line = line.strip()
            
            if line.startswith(">") and not in_block:  # Only capture the first label in each block
                if label and sequence:  # If a sequence was already collected, write it to file
                    outfile.write(f">{label}\n")
                    outfile.write(''.join(sequence) + "\n")
                label = line[1:]  # Store the first label (without ">")
                sequence = []  # Reset the sequence
                in_block = True  # Indicate we are processing a block
            
            elif line.startswith(">") and in_block:
                continue  # Skip additional labels in the same block
            
            elif line.startswith("#"):  # Skip comment lines
                continue
            
            elif line:  # If it's a sequence line, add to the sequence
                sequence.append(line)
            
            else:  # Blank line or end of block, reset for the next block
                in_block = False

        # Write the last sequence to the file
        if label and sequence:
            outfile.write(f">{label}\n")
            outfile.write(''.join(sequence) + "\n")


