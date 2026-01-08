import sys
from Bio import SeqIO

def Split_fasta(input_file, read_len=100, step=50):
    for record in SeqIO.parse(input_file, "fasta"):
        seq_id = record.id
        metadata = record.description[len(seq_id):].strip()
        sequence = str(record.seq)

        n = 1
        for i in range(0, len(sequence) - read_len + 1, step):
            subseq = sequence[i:i + read_len]
            new_header = f"{seq_id}.{n} {metadata}"
            print(f">{new_header}\n{subseq}")
            n += 1


def main():
    # Asegurarse del uso correcto del programa
    if len(sys.argv) != 2:
        print("Uso:\t python3 assembly2reads.py ruta/archivo.fasta", 
              file=sys.stderr)
        sys.exit(1)

    # Split fasta into pseudo-reads
    fasta_file = sys.argv[1]
    Split_fasta(fasta_file)
    
    return 0


if __name__ == "__main__":
    main()
