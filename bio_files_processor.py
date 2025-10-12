def convert_multiline_fasta_to_oneline(input_fasta, output_fasta = 'output_fasta.fasta'):
    with open(input_fasta) as fasta_file, open(output_fasta, 'a') as processed_fasta:
        sequence = ''
        for line in fasta_file:
            if line.startswith(">"):
                if sequence:
                        processed_fasta.write(f'{sequence}\n')
                        sequence = ''
                processed_fasta.write(line)
            else:
                sequence += line.strip()
        if sequence:
            processed_fasta.write(sequence.strip())
