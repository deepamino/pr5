# Estructuras en tres dimensiones con Biopython

<div align="justify">
    
Este documento presenta un conjunto de ejercicios prácticos centrados en la bioinformática estructural y en el uso de herramientas computacionales como Biopython para analizar y manipular datos biológicos, específicamente relacionados con estructuras de proteínas y secuencias.

---

<br>

## **Ejercicio 1: Proteína de la hemoglobina humana**

En primer lugar, el archivo en formato PDB se descarga utilizando la función personalizada `download_pdb`, la cual obtiene el archivo a partir del identificador único de la proteína proporcionado por la base de datos PDB (Protein Data Bank). El código de la función se presenta a continuación:

```python
def download_pdb(pdb_id, route=None):
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    respuesta = requests.get(url)

    if respuesta.status_code == 200:
        if route is None:
            route = f"{pdb_id.upper()}.pdb"
        else:
            if not route.endswith(".pdb"):
                route += ".pdb"
        with open(route, 'w') as archivo:
            archivo.write(respuesta.text)
        print(f"Archivo PDB '{pdb_id.upper()}.pdb' descargado exitosamente en '{route}'.")
        return route
    else:
        raise Exception(f"No se pudo descargar el PDB. Código de estado HTTP: {respuesta.status_code}")
```

La función utiliza el módulo `requests` para enviar una solicitud HTTP al servidor de PDB, valida la respuesta del servidor y, en caso de éxito, guarda el archivo en el directorio especificado por el usuario o en una ubicación predeterminada.

Posteriormente, se emplean dos herramientas para la visualización interactiva de estructuras tridimensionales de proteínas: la librería `nglview` y la librería `py3Dmol`. Estas herramientas permiten explorar las estructuras en detalle y generar representaciones gráficas. En la Figura 1 se muestra la estructura tridimensional de la hemoglobina humana, generada a partir del archivo descargado y visualizada con `py3Dmol`.

<div align="center">
    <img src="results/estructura_1.png" width="60%">
    <p><b>Figura 1.</b> Representación tridimensional de la estructura de la hemoglobina humana.</p>
</div>

Finalmente, el archivo descargado se carga en Biopython para su análisis estructural. Esto se logra creando un objeto de tipo estructura mediante el módulo `PDBParser` de Biopython, como se ilustra en el siguiente código:

```python
parser = PDBParser(QUIET=True)
pdb_code = "1A3N"
structure = parser.get_structure(pdb_code, f"./files/{pdb_code}.pdb")
```

El objeto `structure` contiene toda la información sobre la organización atómica y molecular de la proteína, permitiendo realizar análisis y cálculos adicionales, como los descritos posteriormente.

### (a) Calcula la distancia entre los átomos O del primer y último residuo de la cadena A de la hemoglobina

Para calcular la distancia entre los átomos de oxígeno (O) correspondientes al primer y último residuo de la cadena A, se siguen los pasos descritos a continuación:

1. **Obtención de los residuos**. Se extraen el primer y último residuo de la cadena A utilizando el modelo de la estructura cargada:

```python
model = next(structure.get_models())  # Obtiene el primer modelo
chain_A = model['A']  # Selecciona la cadena A
residues = list(chain_A.get_residues())  # Lista de residuos en la cadena
first_residue = residues[0]  # Primer residuo
last_residue = residues[-1]  # Último residuo
```

2. **Obtención de los átomos de interés**. Se emplea la función `get_atom` para extraer los átomos de oxígeno (O) de los residuos seleccionados:

```python
def get_atom(residue, atom_id):
    try:
        return residue[atom_id]  # Recupera el átomo especificado
    except KeyError:
        print(f"El residuo {residue.get_resname()} {residue.id[1]} no tiene un átomo '{atom_id}'.")
        return None
 ```

3. **Cálculo de la distancia**. Una vez obtenidos los átomos, se calcula la distancia entre ellos en angstroms. Para ello, se define la función `calculate_distance`, que utiliza la operación de resta disponible para los objetos de tipo átomo en Biopython:

```python
def calculate_distance(atom1, atom2):
    distance = atom1 - atom2  # Calcula la distancia en línea recta
    return distance
```

4. **Visualización de la distancia calculada**. La distancia se representa visualmente mediante una imagen tridimensional, donde se destacan los átomos de interés. La **Figura 2** ilustra esta distancia en línea recta entre los átomos O del primer y último residuo:

<div align="center">
    <img src="results/ej_1_a.png" width="60%">
    <p><b>Figura 2.</b> Distancia entre los átomos O del primer y último residuo de la cadena A.</p>
</div>

### (b) Calcula el ángulo diedro entre los átomos N,CA,C, y O del primer residuo de la cadena A

El ángulo diedro es un ángulo de torsión que describe la disposición espacial de cuatro átomos conectados de forma secuencial. Este tipo de ángulo está estrechamente relacionado con las conformaciones locales y las estructuras secundarias de las proteínas, como las hélices alfa y las hojas beta. 

1. **Obtención de los átomos**. Se extraen los átomos especificados (N, CA, C y O) del primer residuo de la cadena A utilizando la función `get_atom`:

```python
atom_N = get_atom(first_residue, 'N')
atom_CA = get_atom(first_residue, 'CA')
atom_C = get_atom(first_residue, 'C')
atom_O = get_atom(first_residue, 'O')
```

2. **Extracción de las coordenadas**. Se obtienen los vectores posicionales de cada uno de los átomos seleccionados:

```python
coord_N = atom_N.get_vector()
coord_CA = atom_CA.get_vector()
coord_C = atom_C.get_vector()
coord_O = atom_O.get_vector()
```

3. **Cálculo del ángulo diedro**. Se utiliza la función `calc_dihedral` de Biopython, que calcula el ángulo de torsión basado en las posiciones de los cuatro átomos. El resultado, inicialmente en radianes, se convierte a grados mediante la función `np.degrees`:

```python
angle = calc_dihedral(coord_N, coord_CA, coord_C, coord_O)
angle_degrees = np.degrees(angle)
```
   
<div align="center">
    <img src="results/ej_1_b.png" width="60%">
    <p><b>Figura 3.</b> Ángulo diedro de los cuatro átomos.</p>
</div>

En este caso, el ángulo diedro calculado es de **-25.67º**. Este valor negativo indica que el ángulo de torsión está en el sentido contrario a las agujas del reloj cuando se observa desde el vector definido por los átomos N y CA hacia el vector definido por los átomos C y O.
Un ángulo diedro de **-25.67º** sugiere una ligera torsión que podría estar asociada con la formación de estructuras secundarias específicas, aunque no es un valor característico de una hélice alfa ni una hoja beta, que suelen tener ángulos más definidos. Este resultado destaca la flexibilidad conformacional del residuo inicial en la cadena A de la hemoglobina, un aspecto importante para su función biológica y estabilidad estructural.

### (c) Calcula el centro de masas de la estructura de la hemoglobina

El **centro de masas (center of mass)** de una estructura molecular se calcula como el promedio ponderado de las coordenadas de todos los átomos en la estructura, donde el peso de cada átomo es su masa. La fórmula es la siguiente:

$$
\text{Centro de masas (COM)} = \frac{\sum (m_i \cdot r_i)}{\sum m_i}
$$

Donde:
- $m_i$ es la masa del átomo $i$ (generalmente el peso atómico de su elemento).
- $r_i$ es el vector de coordenadas (x, y, z) del átomo $i$.
- La suma se realiza sobre todos los átomos de la estructura.

Información de las masas atómicas: [Disponible aquí](https://www.lenntech.es/periodica/masa/masa-atomica.htm)

1. **Identificación de los elementos presentes en la estructura**. Se recorren todos los átomos de la estructura para identificar los elementos únicos presentes. Esto permite asociar las masas atómicas adecuadas a cada elemento:

```python
unique_elements = set()
for atom in structure.get_atoms():
    element = atom.element.strip().upper()
    unique_elements.add(element)
```
En el caso de la hemoglobina, los elementos identificados son: `'FE'`(Hierro), `'O'`(Oxígeno), `'S'`(Sulfuro), `'N'`(Nitrógeno) y `'C'`(Carbono).

2. **Cálculo de las coordenadas ponderadas por la masa**. Para cada átomo, se extraen sus coordenadas ($x, y, z$)) y se ponderan por su masa atómica. Al mismo tiempo, se acumula la masa total de todos los átomos:

```python
atomic_mass = 0.0
sum_x = 0.0
sum_y = 0.0
sum_z = 0.0
for atom in structure.get_atoms():
    element = atom.element.strip().upper()
    masa = get_mass(element)  # Función que retorna la masa atómica
    x, y, z = atom.get_coord()
    atomic_mass += masa
    sum_x += masa * x
    sum_y += masa * y
    sum_z += masa * z
```

3. **Obtención del centro de masas**. Dividiendo las coordenadas ponderadas ($\sum m_i \cdot x, \sum m_i \cdot y, \sum m_i \cdot z$) por la masa total ($\sum m_i$), se obtiene el centro de masas como un vector tridimensional:

```python
centro_masas = np.array([sum_x, sum_y, sum_z]) / atomic_mass
```

4. **Visualización del resultado**. El centro de masas calculado ($14.45, 2.01, 13.18$) se representa gráficamente para resaltar su posición en la estructura de la hemoglobina. En la **Figura 4**, se observa el centro de masas como un punto destacado en la estructura tridimensional:

<div align="center">
    <img src="results/ej_1_c.png" width="60%">
    <p><b>Figura 4.</b> Centro de masas de la hemoglobina representado como un punto destacado.</p>
</div>

En el caso de la hemoglobina, se encuentra cerca del núcleo de la estructura, lo cual es consistente con su configuración compacta y su función biológica de transporte de oxígeno en el organismo. Este análisis puede ser extendido a otras proteínas para evaluar la distribución de masa y su relación con propiedades estructurales y funcionales.

---

<br>

## **Ejercicio 2: Proteína lisozima del huevo de la gallina**

En primer lugar, al igual que en el ejercicio 1, el archivo en formato PDB se descarga utilizando la función personalizada `download_pdb`, la cual obtiene el archivo a partir del identificador único de la proteína proporcionado por la base de datos PDB (Protein Data Bank).

Posteriormente, se emplean dos herramientas para la visualización interactiva de estructuras tridimensionales de proteínas: la librería `nglview` y la librería `py3Dmol`. 

<div align="center">
    <img src="results/estructura_2.png" width="60%">
    <p><b>Figura 5.</b> Representación tridimensional de la estructura de la lisozima.</p>
</div>

Finalmente, el archivo descargado se carga en Biopython para su análisis estructural. Esto se logra creando un objeto de tipo estructura mediante el módulo `PDBParser`. El objeto `structure` contiene toda la información sobre la organización atómica y molecular de la proteína, permitiendo realizar análisis y cálculos adicionales, como los descritos posteriormente.

### (a) Número de átomos y nombres del primer y último átomo en la lista

Este apartado tiene como objetivo determinar el número total de átomos presentes en la estructura de la lisozima y especificar el nombre del primer y último átomo en la lista de átomos de la estructura.

1. **Obtención de la lista de átomos**. Se genera una lista que contiene todos los átomos de la estructura utilizando el método `get_atoms`:

```python
atoms = list(structure.get_atoms())
```

2. **Cálculo del número total de átomos**. La longitud de la lista de átomos se determina mediante la función `len`, lo que proporciona el número total de átomos en la estructura:

```python
num_atoms = len(atoms)
```
En este caso, la proteína de la lisozima contiene **$1.102$** átomos.

3. **Obtención de los nombres del primer y último átomo**. Se acceden al primer y último átomo en la lista mediante índices y se extraen sus nombres utilizando el método `get_name`:

```python
first_atom = atoms[0].get_name()
last_atom = atoms[-1].get_name()
```

En este caso:
- El **primer átomo** es un átomo de nitrógeno (`N`).
- El **último átomo** es un átomo de oxígeno (`O`).

4. **Visualización de los resultados**. La **Figura 6** muestra la estructura tridimensional de la lisozima, donde se destacan claramente el primer y el último átomo:

<div align="center">
    <img src="results/ej_2_a.png" width="60%">
    <p><b>Figura 6.</b> Primer y último átomo de la lisozima.</p>
</div>

El análisis realizado confirma que la estructura de la lisozima contiene un total de **1,102 átomos**. El **primer átomo** identificado es un átomo de nitrógeno (`N`), y el **último átomo** es un átomo de oxígeno (`O`). 

### (b) Cálculo del ángulo entre los tres primeros átomos de la lista

El cálculo del ángulo entre tres átomos se basa en la geometría molecular y utiliza las posiciones de los átomos en el espacio tridimensional. El ángulo se obtiene a partir de los vectores de posición de los tres átomos, calculando el ángulo entre ellos utilizando funciones específicas de Biopython.

1. **Obtención de los vectores de los tres primeros átomos**. Se seleccionan los tres primeros átomos de la lista y se extraen sus vectores de posición tridimensionales ($x, y, z$):

```python
atom1, atom2, atom3 = atoms[0], atoms[1], atoms[2]
vector1 = atom1.get_vector()
vector2 = atom2.get_vector()
vector3 = atom3.get_vector()
```

2. **Cálculo del ángulo en radianes**. Se utiliza la función `calc_angle` de Biopython, que calcula el ángulo formado por tres vectores en radianes:

```python
angle = calc_angle(vector1, vector2, vector3)
```

3. **Conversión del ángulo a grados**. El ángulo calculado en radianes se convierte a grados utilizando la función `math.degrees`:

```python
angle_deg = math.degrees(angle)
```

En este caso, el ángulo calculado es **118.36 grados**.

4. **Visualización del ángulo**. La **Figura 7** ilustra gráficamente la disposición de los tres primeros átomos y el ángulo calculado entre ellos:

<div align="center">
    <img src="results/ej_2_b.png" width="60%">
    <p><b>Figura 7.</b> Representación del ángulo entre los tres primeros átomos de la lisozima.</p>
</div>

### (c) Identificación de la cadena y el residuo del átomo central de la lista

Este apartado tiene como objetivo identificar el átomo que se encuentra en la posición central de la lista de átomos, así como su cadena y residuo correspondiente. Este análisis se basa en las estructuras de datos y clases proporcionadas por Biopython.

1. **Obtención del índice y del átomo central**. El índice del átomo central se calcula dividiendo el número total de átomos (`num_atoms`) entre dos, utilizando la división entera para garantizar un índice válido. A partir de este índice, se obtiene el átomo central:

```python
middle_index = num_atoms // 2
middle_atom = atoms[middle_index]
```

2. **Obtención de la cadena y el residuo del átomo**
Mediante la jerarquía de clases de Biopython, se obtienen los elementos relacionados con el átomo:
- **Residuo**: Se obtiene utilizando el método `get_parent` aplicado al átomo.
- **Cadena**: Se obtiene aplicando `get_parent` al residuo.
- **Información del residuo**: A partir del objeto residuo, se extraen el nombre del residuo (`get_resname`) y su número (`get_id}[1]`).
- **ID de la cadena**: Se obtiene con el método `get_id` del objeto cadena.

```python
residue = middle_atom.get_parent()
chain = residue.get_parent()
chain_id = chain.get_id()
residue_id = residue.get_id()
residue_name = residue.get_resname()
residue_number = residue_id[1]
```

4. **Resultados**  
En este caso, el átomo central se encuentra en la posición **551** de la lista de átomos. Este átomo es un **nitrógeno (`N`)** que pertenece:
- A la **cadena `A`** (única en la estructura).
- Al **residuo `PRO` (prolina)** con número **70**.

5. **Visualización del átomo central**. La **Figura 8** muestra gráficamente la ubicación del átomo central y su residuo correspondiente dentro de la estructura tridimensional:

<div align="center">
    <img src="results/ej_2_c.png" width="60%">
    <p><b>Figura 8.</b> Representación del átomo central y su residuo.</p>
</div>

El átomo central identificado es un **nitrógeno (`N`)** del residuo **prolina (`PRO`)** en la **cadena `A`** de la proteína. La prolina es un aminoácido que desempeña un papel especial en las estructuras proteicas debido a su geometría rígida, que puede influir en la flexibilidad y estabilidad de la proteína. 

---

<br>

### README: Ejercicio 3 - Proteínas del Sueño

#### Introducción
En este ejercicio, se exploran las características estructurales, funcionales y evolutivas de dos proteínas clave implicadas en la regulación del sueño y otros procesos biológicos esenciales: **Orexina-A/Hipocretina-1 (1WSO)** y **Orexina-B/Hipocretina-2 (1CQ0)**. A través del análisis tridimensional, búsquedas de homólogos y la construcción de árboles filogenéticos, se busca entender las similitudes y diferencias entre estas proteínas, sus relaciones evolutivas y su papel en organismos como el ser humano.

La Orexina-A y la Orexina-B son neuropéptidos producidos en el hipotálamo que desempeñan un papel crítico en la regulación del ciclo sueño-vigilia, el apetito y la homeostasis energética. Aunque comparten un precursor común, presentan diferencias notables en su estructura y funcionalidad. Este análisis aborda tanto las similitudes que permiten su funcionalidad compartida como las diferencias que explican sus afinidades específicas por distintos receptores.

---

### Visualización y Exploración de las Estructuras Tridimensionales
Las proteínas 1WSO y 1CQ0 se cargaron utilizando la librería `Bio.PDB` y se visualizaron mediante herramientas como `nglview` y `py3Dmol`. Estas visualizaciones revelaron detalles importantes:

1. **Estructura de 1WSO (Orexina-A):**
   - La proteína muestra una hélice alfa prominente, con una organización compacta y una estructura secundaria bien definida. Esta hélice alfa es crucial para su interacción con los receptores OX1R y OX2R.
   - Se identificaron residuos no estándar, como PCA y NH₂. Estos residuos, aunque pequeños, desempeñan un papel importante en la estabilización de la estructura o en su funcionalidad experimental.

2. **Estructura de 1CQ0 (Orexina-B):**
   - A diferencia de 1WSO, esta proteína es más compacta y tiene menos flexibilidad en sus extremos terminales, lo que refuerza su interacción específica con el receptor OX2R.
   - La estructura presenta una única hélice alfa continua, que también es fundamental para sus funciones biológicas.

Además, se identificaron y describieron ligandos presentes en las estructuras, con un análisis detallado de su rol en la funcionalidad y estabilidad de las proteínas.

---

### Alineación de Proteínas
El alineamiento estructural se llevó a cabo utilizando el algoritmo de superposición de Kabsch implementado en la clase `PDB.Superimposer` de Biopython. Este proceso minimiza la suma de las distancias cuadradas entre átomos equivalentes en ambas proteínas, generando tres componentes clave:

1. **Matriz de Rotación**:
   $$
   R = 
   \\begin{bmatrix}
   1.00000000 & -1.31074131 \\times 10^{-8} &  9.74082657 \\times 10^{-8} \\\\
   1.31074154 \\times 10^{-8} & 1.00000000 & -2.44023484 \\times 10^{-8} \\\\
   -9.74082652 \\times 10^{-8} & 2.44023494 \\times 10^{-8} & 1.00000000
   \\end{bmatrix}
   $$
   Esta matriz muestra que no fue necesario un cambio significativo en la orientación de las proteínas para lograr la alineación, ya que los valores cercanos a 1 en la diagonal indican una conservación espacial notable.

2. **Vector de Traslación**:
   $$
   T = 
   \\begin{bmatrix}
   -1.10392905 \\times 10^{-6} \\\\
   1.58119236 \\times 10^{-7} \\\\
   -7.35061366 \\times 10^{-7}
   \\end{bmatrix}
   $$
   Este vector describe el desplazamiento mínimo necesario para alinear las estructuras en un mismo sistema de coordenadas. Los valores muy pequeños indican que las proteínas ya estaban alineadas de manera aproximada antes del ajuste.

3. **RMSD (Root Mean Square Deviation)**:
   $$
   \\text{RMSD} = \\sqrt{\\frac{1}{N} \\sum_{i=1}^N \\left\\| P_{1i} - P_{2i} \\right\\|^2} = 7.47 \\, \\text{Å}
   $$
   Este valor refleja diferencias estructurales significativas, particularmente en las regiones terminales y flexibles de las proteínas.

La visualización de las proteínas alineadas mostró una conservación en las hélices alfa principales y diferencias en las regiones periféricas, destacando tanto similitudes funcionales como adaptaciones estructurales específicas.

---

### Búsqueda de Homólogos y Construcción de Árboles Filogenéticos
Se utilizó BLAST para identificar proteínas homólogas en la base de datos NR. Las búsquedas para 1WSO y 1CQ0 identificaron múltiples homólogos, que posteriormente se utilizaron para construir árboles filogenéticos con los métodos Neighbor Joining (NJ) y UPGMA. Los resultados incluyeron:

- **1WSO (Orexina-A):** La proteína más cercana fue la Orexina A y B de Ovis aries, lo que sugiere una fuerte conservación evolutiva.
- **1CQ0 (Orexina-B):** Su homólogo más cercano fue la cadena L de 7L1U del Homo sapiens, lo que confirma su funcionalidad específica en esta especie.

El análisis de los árboles filogenéticos incluyó métricas como el número de terminales, las distancias entre nodos y los valores de bootstrap, proporcionando una visión detallada de las relaciones evolutivas.

---

### Cálculo de Máximas Distancias entre Átomos
Se desarrolló una función para calcular la distancia máxima entre átomos en cada estructura. Los resultados indicaron que 1WSO tiene una distancia mayor, reflejando una mayor flexibilidad y tamaño en comparación con 1CQ0. Esto resalta diferencias en su estabilidad estructural y sus interacciones con receptores.

La fórmula utilizada para calcular la distancia entre átomos fue:
$$
\\text{Distancia} = \\sqrt{(x_1 - x_2)^2 + (y_1 - y_2)^2 + (z_1 - z_2)^2}
$$

---

### Conclusión
Este análisis detallado de las proteínas 1WSO y 1CQ0 demuestra cómo herramientas computacionales pueden desentrañar complejidades estructurales y evolutivas. La conservación de las hélices alfa principales subraya su importancia funcional, mientras que las diferencias en regiones flexibles reflejan adaptaciones específicas. Estas proteínas son un ejemplo fascinante de cómo la evolución puede preservar funciones esenciales mientras ajusta estructuras para satisfacer necesidades biológicas específicas.

</div>
