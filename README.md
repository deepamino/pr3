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
- Filtrado de los resultados y guardar los atributos necesarios en formato JSON. Enlace a la documentación: [BLAST Record](https://biopython.org/docs/1.75/api/Bio.Blast.NCBIWWW.html)

```python
def search_protein_sequence(output_file="./results/blast_results.json"):

    protein_sequence = get_protein()

    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)
    blast_records = NCBIXML.parse(result_handle)

    filtered_results = []

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                # Filtrar los resultados según el E-value
                if hsp.expect < 0.001:
                    # Guardar los atributos requeridos en formato JSON
                    result = {
                        "id": alignment.hit_id,
                        "length": alignment.length,
                        "evalue": hsp.expect,
                        "identity_percentage": (hsp.identities / hsp.align_length) * 100
                    }
                    filtered_results.append(result)

    # Guardar en el archivo
    with open(output_file, "w") as f:
        json.dump(filtered_results, f, indent=4)

    print(f"Resultados guardados en formato JSON en: {output_file}")
```

**3. Implementación *offline***

Al igual que en e ejercicio 1, se realiza la búsqueda local mediante dos bases de datos de BLAST, disponibles en el NCBI (*National Center for Biotechnology Information*):
- **Base de datos `swisprot`**: repositorio de secuencias de proteínas que forma parte de la base de datos UniProt.
- **Base de datos `nr`**:  base de datos no redundante de secuencias de proteínas. Debido a las limitaciones de recursos, se han descargado únicamente los archivos desde el 000 al 020.

```python
def search_local_protein(output_file = "./results/blast_results_locally.json"):
    # Rutas necesarias para la búsqueda
    os.environ["CMD"] = input("Introduce la ruta del ejecutable de BLAST: ")
    os.environ["BLASTDB"] = input("Introduce la ruta de la base de datos de BLAST: ")

    # Obtención de la secuencia
    sequence = input("Introduce la secuencia de proteína en formato FASTA: ")

    # Escribir la frecuencia para realizar la búsqueda
    fasta_file = "query.fasta"
    with open(fasta_file, "w") as f:
        f.write(">query\n")
        f.write(sequence)


    # Configuración y ejecución de la línea de comando BLAST
    blastp_cline = NcbiblastpCommandline(
        cmd=os.environ["CMD"],
        query=fasta_file,
        db=os.environ["BLASTDB"],
        evalue=0.001,
        outfmt=5,
        out="results.xml"
    )

    blastp_cline()


    with open("results.xml") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        filtered_results = []

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Filtrar los resultados según el E-value
                    if hsp.expect < 0.001:
                        # Guardar los atributos necesarios
                        result = {
                            "id": alignment.hit_id,
                            "length": alignment.length,
                            "evalue": hsp.expect,
                            "identity_percentage": (hsp.identities / hsp.align_length) * 100
                        }
                        filtered_results.append(result)

    # Guardar en el archivo
    with open(output_file, "w") as f:
        json.dump(filtered_results, f, indent=4)

    print(f"Resultados guardados en formato JSON en: {output_file}")

    os.remove(fasta_file)
    os.remove("results.xml")
```

**4. Comparación de la implementación local vs *online***

Al igual que el ejercicio 1, se ha decidido realizar una comparativa entre la implementación local y online. Para ello, se utilizará parte de la secuencia de aminoácidos correspondiente a la proteína de la insulina humana.

...`LVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGP`...

En este caso, lo único que se puede comparar es el tiempo de ejecución, ya que los resultados más fiables se consideran los correspondientes a la base de datos online puesto que es la que más información contiene. Por tanto, la ejecución local en la base de datos de Swissprot es la más rápida puesto que es aquella que contiene menos información. Sin embargo, a pesar de que la online es la que más tarda, a su vez es la más completa puesto que realiza la bñusqueda en toda la base de datos de proteínas `nr`. En definitiva, si se busca una mayor precisión, se recomienda la búsqueda online, mientras que si se busca rapidez, se recomienda la búsqueda local en la base de datos de Swissprot. Además, si se disponen de los recursos necesarios, la base de datos `nr` de forma local y completa también es una buena opción.

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
      <td>363.872967</td>
    </tr>
    <tr>
      <td>local swissprot</td>
      <td>14.235358</td>
    </tr>
    <tr>
      <td>local nr subset</td>
      <td>332.662034</td>
    </tr>
  </tbody>
</table>
</div>

<hr>

</div>
