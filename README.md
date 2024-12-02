# Práctica 3: Alineamiento con BLAST (*Basic Local Alignment Search Tool*) y Biopython

**Participantes:**
- Ricardo Juan Cárdenes Pérez
- Susana Suárez Mendoza

<div align="justify">

Esta práctica consta de tres ejercicios que ponen en práctica la búsqueda mediente la busqueda BLAST de una secuencia mediante Biopython de forma *online* y *offline*. 

## Instalación de BLAST de forma local

1. Se descarga la versión para el equipo de BLAST desde el siguiente enlace: [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
2. Se ejecuta, en nuestro caso, el .exe descargado y se instala en el equipo.
3. Se configura el path de BLAST como variable de entorno del sistema.
4. Se verifica la instalación con el comando blastn --version.
5. Los diversos ejercicios necesitarán las rutas correspondientes a:
   - La ruta del ejecutable de BLAST según la naturaleza, por ejemplo, blastn o blastp.
   - La ruta de la base de datos que se quiera utilizar.

## Ejercicio 1:

Escribe un programa en Python que use Biopython para hacer una búsqueda BLAST de una secuencia de ADN que introduzcas por teclado. El programa debe mostrar por pantalla el número de resultados obtenidos, el E-value del mejor resultado y la descripción de la secuencia más similar. Hágalo de forma online y local.

**1. Obtención de la secuencia introducida por teclado**

En este apartado se depura la secuencia de ADN introducida por el usuario (`get_sequence()`):
- Se suprimen los espacios: `check_no_spaces(sequence)`
- Se traduce a mayúsculas: `sequence_upper(sequence)`
- En el caso de que se contenga cualquier otra *string* que no corresponda a un nucleótido o la longitud de la secuencia sea 0, se lanzará una excepción.

```python
def get_sequence():
    sequence = input("Introduce la secuencia de ADN: ").strip()
    sequence = check_sequence(sequence)
    return sequence

def check_sequence(sequence):
    sequence = check_no_spaces(sequence)
    sequence = sequence_upper(sequence)
    if not check_nucleotides_DNA(sequence) and len(sequence) == 0:
        raise ValueError("Invalid DNA sequence")
    return sequence
```

**2. Implementación *online***

Se crea la función `search_online_sequence`, la cual es la encargada de realizar lo siguiente:
- Obtención de la secuencia: `get_sequence()`
- Realizar la búsqueda BLAST *online*: `NCBIWWW.qblast("blastn", "nt", sequence)` donde "blastn" es el tipo de búsqueda BLAST que se realiza, y se utiliza para comparar secuencias de nucleótidos (ADN o ARN), "nt" es la base de datos de nucleótidos en la que se va a buscar la secuencia y sequence es la secuencia genética que queremos comparar.
- Obtención de los valores y mostrarlos por pantalla.

```python
def search_online_sequence():
    sequence = get_sequence()

    result = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_records = NCBIXML.read(result)
    
    if not blast_records.alignments:
        print("No se encontraron resultados.")

    # La lista blast_records.alignments contiene todos los alineamientos encontrados
    num_results = len(blast_records.alignments)

    best_alignment = blast_records.alignments[0] # primer alineamiento de la lista de alineamientos (generalmente el mejor)
    best_hsp = best_alignment.hsps[0] # High-scoring Segment Pair que es una comparación entre sub-secuencias de la consulta y la base de datos
    e_value = best_hsp.expect # obtención del e-value del mejor alineamiento
    description = best_alignment.title # descripción del mejor alineamiento
    
    print("\nResultados BLAST:")
    print(f"Número de resultados: {num_results}")
    print(f"E-value del mejor resultado: {e_value}")
    print(f"Descripción de la secuencia más similar: {description}")
```

**3. Implementación *offline***

Se realiza la búsqueda local mediante dos bases de datos de BLAST, disponibles en el NCBI (*National Center for Biotechnology Information*):
- **Base de datos `human_genome`**: Base de datos de ensamblaje actual del genoma humano refseq (GRCh). [Información disponible](https://www.ncbi.nlm.nih.gov/books/NBK62345/#blast_ftp_site.The_blastdb_subdirectory)
- **Base de datos `core_nt`**:  La base de datos core_nt es nt sin la mayoría de las secuencias de cromosomas eucariotas. La mayoría de las búsquedas BLAST de nucleótidos con core_nt serán similares a las de la base de datos nt. Sin embargo, core_nt es mejor que nt para alcanzar los objetivos de búsqueda BLAST más comunes, como la identificación de secuencias relacionadas con genes, como secuencias de transcripción y cromosomas bacterianos completos. Esto se debe a que, en los últimos años, nt ha adquirido más contenido de baja relevancia, no anotado y no génico. Para esta base de datos, debido a las limitaciones de recursos, se han descargado los archivos desde el 00 al 20. [Información disponible](https://bioinformaticsonline.com/news/view/44640/new-blast-core-nucleotide-database-core-nT)

```python
def search_local_sequence():
    # Rutas necesarias para la búsqueda
    os.environ["CMD"] = input("Introduce la ruta del ejecutable de BLAST: ")
    os.environ["BLASTDB"] = input("Introduce la ruta de la base de datos de BLAST: ")

    # Obtención de la secuencia
    sequence = get_sequence()

    # Escribir la frecuencia para realizar la búsqueda
    fasta_file = "query.fasta"
    with open(fasta_file, "w") as f:
        f.write(">query\n")
        f.write(sequence)

    # Configuración y ejecución de la línea de comando BLAST
    blastn_cline = NcbiblastnCommandline(
        cmd=os.environ["CMD"],
        query=fasta_file,
        db=os.environ["BLASTDB"],
        evalue=0.001,
        outfmt=5,
        out="results.xml"
    )

    blastn_cline()

    # Leer y procesar los resultados de BLAST de forma similar a la ejecución online
    with open("results.xml") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            if blast_record.alignments:
                print(f"Número de resultados: {len(blast_record.alignments)}")

                best_alignment = blast_record.alignments[0]
                best_hsp = best_alignment.hsps[0]
                print(f"E-value del mejor resultado: {best_hsp.expect}")
                print(f"Descripción del mejor resultado: {best_alignment.title}")
            else:
                print("No se encontraron coincidencias.")

    # Eliminar los archivos temporales
    os.remove(fasta_file)
    os.remove("results.xml")
```

**4. Comparación de la implementación local vs *online***

Debido a las limitaciones de recursos, se decide comparar los resultados entre las diversas bases de datos, puesto que la implementación online utiliza la base de datos completa de `nt` y la en implementación local, únicamente se dispone de un subconjunto de la base de datos `core_nt` y `human_genome`. Además, se comparará el tiempo de ejecución entre las implementaciones.

Para ello, se decide alinear mediante BLAST parte de una secuencia conocida que corresponde al gen de la hemoglobina:
...`TACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTGCTGGTGGTCTACCCTT`...

En definitiva, la implementación local con la base de datos core_nt es la que más se aproxima a la versión en línea, ya que realiza búsquedas en parte de la base de datos nt. En cambio, la base de datos human_genome genera menos resultados, debido a su enfoque específico en el genoma humano.
Sin embargo, la búsqueda local, debido a la falta de recursos, es la que más tarda. Por lo que las búsquedas locales se deben realizar en el caso de que se dispongan de los recursos necesarios.

<div align="center">
<table border="1" style="border-collapse: collapse; text-align: center;">
  <thead>
    <tr>
      <th>Implementación</th>
      <th>Tiempo de ejecución (s)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>online</td>
      <td>67.487276</td>
    </tr>
    <tr>
      <td>local human genome</td>
      <td>19.800459</td>
    </tr>
    <tr>
      <td>local subconjunto core_nt</td>
      <td>524.341397</td>
    </tr>
  </tbody>
</table>
</div>

<hr>

## Ejercicio 2:

Escribe un programa en Python que use Biopython para hacer una búsqueda BLAST de una secuencia de proteína que introduzcas por teclado. El programa debe guardar en un fichero los resultados que tengan un E-value menor que 0.001. El fichero debe contener el identificador, la longitud, el E-value y el porcentaje de identidad de cada resultado. Hágalo de forma online y local.

**1. Obtención de la secuencia introducida por teclado**

De forma similar al ejercicio 1, se depura la secuencia de aminoácidos introducida por el usuario con la salvedad de que en lugar de comprobar nucleótidos se comprueban aminoácidos:

```python
def check_protein(sequence):
    sequence = check_no_spaces(sequence)
    sequence = sequence_upper(sequence)
    if not check_aminoacids_protein(sequence) and len(sequence) == 0:
        raise ValueError("Invalid protein sequence")
    return sequence

def get_protein():
    protein = input("Enter the protein sequence: ")
    protein = check_protein(protein)
    return protein
```

**2. Implementación *online***

Se crea la función `search_protein_sequence`, la cual es la encargada de realizar lo siguiente:
- Obtención de la secuencia de aminoácidos: `get_protein()`
- Realizar la búsqueda BLAST *online*: `NCBIWWW.qblast("blastp", "nr", protein_sequence)` donde "blastp" indica que se está realizando una búsqueda de proteínas (BLAST para proteínas), "nr" se refiere a la base de datos de proteínas non-redundantes (NR) proporcionada por NCBI, que contiene un gran número de secuencias proteicas de diversas especies y "protein_sequence" es la secuencia proteica que se está buscando en la base de datos.
- Filtrado de los resultados y guardarlos en un archivo json.

</div>
