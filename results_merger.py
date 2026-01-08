"""
╭────────────────────────────────────────────────────────────────────────────╮
│                           ~ results_merger.py ~                            │
│                                                                            │
│ Author:     Oleanna Guerra Font                                            │
│                                                                            │
│ Date:       27/10/2025                                                     │
│                                                                            │
│ Version:    1.3                                                            │
│                                                                            │
│ Description:                                                               │
│    Program that merges results from other programs, using contig names as  │
│        keys, unifying everything into a single tab-delimeted file.         │
╰────────────────────────────────────────────────────────────────────────────╯
"""
import argparse
import csv
import os
import re
import sys

# Increase the CSV field size limit
csv.field_size_limit(sys.maxsize)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Merge all contig results")
    
    # Required args
    parser.add_argument('-i', '--input', required=True, 
                        help='Path to contigs.txt file')
    parser.add_argument('-r', '--results', required=True, 
                        help='Path to results file')
    parser.add_argument('-p', '--program', required=True, 
                        help='Program that generated the results file. Options:\
                            \n\t"kraken2", "barrnap", "samtools_cov", "gtdbtk"\
                            \n\t"whokaryote", "metabat2"')

    # Optional args
    parser.add_argument('-o', '--output', default='contigs_results.csv', 
                        help='Output CSV file')
    parser.add_argument('-d', '--directory', default='.', 
                        help='Output directory')

    args = parser.parse_args()

    return args


def read_contigs_file(contigs_path):
    contigs = []
    with open(contigs_path, 'r') as f:
        # skip header (first line)
        header = f.readline()
        
        # Save each contig id into a list
        for line in f:
            contig_id = line.strip().split()[0]
            contigs.append(contig_id)
    return contigs


def parse_kraken2_results(results_path):
    dicc_kraken2 = {}
    # open results_path and read it
    with open(results_path, 'r') as f:
        print(f'Proccessing results from "{results_path}"')
        
        for line in f:
            # Skips empty or incomplete lines
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            
            # Extract contig and taxid and add to dictionary
            contig, taxid = parts[1].strip(), parts[2].strip()
            # skip if unclassified (taxid = 0)
            if taxid == "0":
                continue
            
            seqid = re.search(r"(.*\.\d)_|(.*\.\d)", contig).group(1)
            if not seqid:
                seqid = re.search(r"(.*\.\d)_|(.*\.\d)", contig).group(2)
            
            if seqid not in dicc_kraken2:
                entry = []
                entry.append(f"{taxid}:{contig}")
            else:
                # else: append
                entry.append(f"{taxid}:{contig}")
            
            dicc_kraken2[seqid] = entry
            
    return dicc_kraken2


def parse_barrnap_results(results_path):
    dicc_barrnap = {}
    with open(results_path, 'r') as f:
        print(f'Proccessing results from "{results_path}"')
        
        for line in f:
            # Skips empty lines and comments
            if line.startswith('#') or not line.strip():
                continue
            
            # Divides each line into a list of its components
            gff_cols = line.strip().split('\t')
            
            # Each line should contain 9 columns
            if len(gff_cols) < 9:
                continue
            
            # Extract important attributes from the gff file
            seqid = gff_cols[0]
            location = f"{gff_cols[3]}:{gff_cols[4]}({gff_cols[6]})"
            
            # Extract product and add to dictionary
            product = re.search(r"product=([^;\n]+)", line).group(1)
            entry = f"{location};{product}"
            if seqid not in dicc_barrnap:
                arn = []
                arn.append(entry)
            else:
                arn.append(entry)
            dicc_barrnap[seqid] = arn   
    return dicc_barrnap


def parse_samtoolsCov_results(results_path):
    dicc_samcov = {}
    
    with open(results_path, 'r') as f:
        print(f'Proccessing results from "{results_path}"')
        
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            # Divides each line into a list of its components
            parts = line.strip().split('\t')
            
            # Each line should contain 9 columns
            if len(parts) < 9:
                continue
            
            # Extract info and fill dictionary
            seqid, coverage, mdepth = parts[0].strip(), parts[5].strip(), parts[6].strip()
            entry = f"coverage={coverage}%;mean_depth={mdepth}x"
            dicc_samcov[seqid] = entry

    return dicc_samcov


def parse_metabat_results(results_path):
    dicc_metabat = {}
    # open results_path and read it
    with open(results_path, 'r') as f:
        print(f'Proccessing results from "{results_path}"')
        
        for line in f:
            # Skips empty or uncomplete lines
            if not line.strip():
                continue

            if line.startswith('>'):
                parts = line[1:].strip().split('\t')
                contig = parts[0]
                bin = parts[1]
            else:
                continue
            
            dicc_metabat[contig] = bin
    
    return dicc_metabat


def parse_whokaryote_results(results_path):
    dicc_whokaryote = {}
    # open results_path and read it
    with open(results_path, 'r') as f:
        print(f'Proccessing results from "{results_path}"')
        
        # skip header (first line)
        header = f.readline()
        
        for line in f:
            # Skips empty or uncomplete lines
            if not line.strip():
                continue

            parts = line[1:].strip().split('\t')
            contig = parts[0]
            predicted = parts[1]

            dicc_whokaryote[contig] = predicted
    
    return dicc_whokaryote


def parse_gtdbtk_results(results_path):
    dicc_gtdbtk = {}
    # open results_path and read it
    with open(results_path, 'r') as f:
        print(f'Proccessing results from "{results_path}"')
        
        # skip header (first line)
        header = f.readline()
        
        for line in f:
            # Skips empty or uncomplete lines
            if not line.strip():
                continue
            
            # Divides each line into a list of its components
            parts = line.strip().split('\t')
            
            # Each line should contain 20 columns
            if len(parts) < 20:
                continue
            
            # Extract info and fill dictionary
            seqid, classification, warnings = parts[0].strip(), parts[1].strip(), parts[19].strip()
            entry = [classification, warnings]
            dicc_gtdbtk[seqid] = entry
    
    return dicc_gtdbtk


def get_parser(program_name):
    # Ensure the correct use of the program
    program_name = program_name.lower()
    if program_name in ['kraken2', 'k2']:
        return parse_kraken2_results, 'Kraken2'
    elif program_name in ['barrnap', 'b']:
        return parse_barrnap_results, 'Barrnap'
    elif program_name in ['samtools_cov', 'cov']:
        return parse_samtoolsCov_results, 'Samtools_coverage'
    elif program_name in ['metabat', 'metabat2', 'm']:
        return parse_metabat_results, 'Metabat2'
    elif program_name in ['gtdbtk', 'g']:
        return parse_gtdbtk_results, 'GTDBtk'
    elif program_name in ['whokaryote', 'w']:
        return parse_whokaryote_results, 'Whokaryote'
    else:
        # Will raise an error if the program name provided is not supported
        raise ValueError(f"Unsupported program: {program_name}")


def merge_results(contigs, results_dict, program_column, output_path):
    # Dictionary to hold all contig data
    data = {contig: {'Seqid': contig} for contig in contigs}
    existing_fieldnames = set(['Seqid'])  # Start with the default column
    
    # If file exists, load it and update the data dictionary
    if os.path.exists(output_path):
        with open(output_path, 'r', newline='', encoding='utf-8') as infile:
            reader = csv.DictReader(infile, delimiter='\t')
            existing_fieldnames.update(reader.fieldnames or [])
            for row in reader:
                contig = row['Seqid']
                if contig not in data:
                    data[contig] = {'Seqid': contig}
                data[contig].update(row)

    # Update/add new results
    for contig in contigs:
        data[contig][program_column] = results_dict.get(contig, 'NA')

    # Prepare final list of columns (replace old column if it existed)
    final_columns = ['Seqid'] + sorted(c for c in existing_fieldnames.union([program_column]) if c != 'Seqid')
    
    # Write updated data back to file in tab-delimited format
    with open(output_path, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=final_columns, delimiter='\t')
        writer.writeheader()
        for contig in contigs:
            writer.writerow(data[contig])


def main():
    args = parse_arguments()
    contigs = read_contigs_file(args.input)

    parser_func, program_column = get_parser(args.program)
    results_dict = parser_func(args.results)

    output_path = os.path.join(args.directory, args.output)
    merge_results(contigs, results_dict, program_column, output_path)
    
    # Indicate what the program has just done
    print(f"{str(args.program).capitalize()} results written to: {output_path}\t")


if __name__ == "__main__":
    main()

