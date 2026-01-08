import re
from pathlib import Path

def main():
    #Path
    gff_files = "../Barrnap_dir2/*.gff"
    
    # Patrón a buscar en cada línea
    patron = "Name=([^;]+);product=([^;\n]+)(?:;note=([^;\n]+))?"
    
    for f in gff_files:
        # Abrir fichero gff de Barrnap
        input = open(f, "r")
        
        # Diccionario para almacenar la información relevante optenida
        dicc_barrnap = {}
        
        # Output file open
        filename = Path(f).stem
        output = filename + ".tsv"
        output = open(output, "w")
        
        # Leer fichero gff y extraer los resultados
        for line in input:
            if line.startswith("#"):
                continue
            
            # Obtener el match object de regex
            match = re.search(patron, line)
            
            # Del match object obtener info del rrna encontrado
            name = match.group(1)
            product = match.group(2)
            note = match.group(3) if match.group(3) else "NA" # Note puede estar vacío
            
            # Obtener el identificador de la secuencia (seqid)
            seqid = re.search("(.+)(?=\tbarrnap)", line).group(1)
            
            # Escribir línea en output
            
        input.close()
    
    # Escribir resultados Barrnap en un tsv
    with open("Barrnap.tsv", "w") as out:
        out.write("seqid\trna_type\n")
        for key in dicc_barrnap.keys():
            for i in dicc_barrnap[key]:
                out.write(f"{key}\t{i}\n")
    
    return 0


if __name__ == "__main__":
    main()
