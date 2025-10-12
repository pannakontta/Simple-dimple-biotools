def convert_multiline_fasta_to_oneline(input_fasta, output_fasta = 'output_fasta.fasta'):
    with open(input_fasta) as fasta_file, open(output_fasta, 'a') as processed_fasta:
        sequence = ''
        for line in fasta_file:
            if line.startswith(">"):
                if sequence:
                        processed_fasta.write(f'{sequence}\n')  # write the sequence in one line to output file
                        sequence = ''
                processed_fasta.write(line)  # write the secuence id to output file
            else:
                sequence += line.strip()
        if sequence:
            processed_fasta.write(sequence.strip())  # write the last sequence


def parse_blast_output(input_file, output_file = 'best_proteins'):
    with open(input_file) as blast_results, open(output_file, 'a') as best_proteins:
        query = False
        for line in blast_results:
            if query:  # process the next line after "Description"
                if "]" not in line:
                    line = line.replace("...", "]")
                    gene = line.split(']')[0]
                else:
                    gene = line.split(']')[0] + "]"
                best_proteins.write(f'{gene}\n')
            if line.startswith("Description"):
                query = True
            else:
                query = False
