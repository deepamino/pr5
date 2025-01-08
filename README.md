# Estructuras en tres dimensiones con Biopython

Este documento presenta un conjunto de ejercicios prácticos centrados en la bioinformática estructural y en el uso de herramientas computacionales como Biopython para analizar y manipular datos biológicos, específicamente relacionados con estructuras de proteínas y secuencias.

---

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

### (a) Calcula la distancia entre los átomos O del primer y último residuo de la cadena A de la hemoglobina.

Para calcular la distancia entre los átomos de oxígeno (O) correspondientes al primer y último residuo de la cadena A, se siguen los pasos descritos a continuación:

1. **Obtención de los residuos**  
   Se extraen el primer y último residuo de la cadena A utilizando el modelo de la estructura cargada:

   ```python
   model = next(structure.get_models())  # Obtiene el primer modelo
   chain_A = model['A']  # Selecciona la cadena A
   residues = list(chain_A.get_residues())  # Lista de residuos en la cadena
   first_residue = residues[0]  # Primer residuo
   last_residue = residues[-1]  # Último residuo
   ```

2. **Obtención de los átomos de interés**  
   Se emplea la función `get_atom` para extraer los átomos de oxígeno (O) de los residuos seleccionados:

   ```python
   def get_atom(residue, atom_id):
       try:
           return residue[atom_id]  # Recupera el átomo especificado
       except KeyError:
           print(f"El residuo {residue.get_resname()} {residue.id[1]} no tiene un átomo '{atom_id}'.")
           return None
   ```

3. **Cálculo de la distancia**  
   Una vez obtenidos los átomos, se calcula la distancia entre ellos en angstroms. Para ello, se define la función `calculate_distance`, que utiliza la operación de resta disponible para los objetos de tipo átomo en Biopython:

   ```python
   def calculate_distance(atom1, atom2):
       distance = atom1 - atom2  # Calcula la distancia en línea recta
       return distance
   ```

4. **Visualización de la distancia calculada**  
   La distancia se representa visualmente mediante una imagen tridimensional, donde se destacan los átomos de interés. La **Figura 2** ilustra esta distancia en línea recta entre los átomos O del primer y último residuo:

   <div align="center">
       <img src="results/ej_1_a.png" width="60%">
       <p><b>Figura 2.</b> Distancia entre los átomos O del primer y último residuo de la cadena A.</p>
   </div>

### (b) Calcula el ángulo diedro entre los átomos N,CA,C, y O del primer residuo de la cadena A.

El ángulo diedro es un ángulo de torsión que describe la disposición espacial de cuatro átomos conectados de forma secuencial. Este tipo de ángulo está estrechamente relacionado con las conformaciones locales y las estructuras secundarias de las proteínas, como las hélices alfa y las hojas beta. 

1. **Obtención de los átomos**  
   Se extraen los átomos especificados (N, CA, C y O) del primer residuo de la cadena A utilizando la función `get_atom`:

   ```python
   atom_N = get_atom(first_residue, 'N')
   atom_CA = get_atom(first_residue, 'CA')
   atom_C = get_atom(first_residue, 'C')
   atom_O = get_atom(first_residue, 'O')
   ```

2. **Extracción de las coordenadas**  
   Se obtienen los vectores posicionales de cada uno de los átomos seleccionados:

   ```python
   coord_N = atom_N.get_vector()
   coord_CA = atom_CA.get_vector()
   coord_C = atom_C.get_vector()
   coord_O = atom_O.get_vector()
   ```

3. **Cálculo del ángulo diedro**  
   Se utiliza la función `calc_dihedral` de Biopython, que calcula el ángulo de torsión basado en las posiciones de los cuatro átomos. El resultado, inicialmente en radianes, se convierte a grados mediante la función `np.degrees`:

   ```python
   angle = calc_dihedral(coord_N, coord_CA, coord_C, coord_O)
   angle_degrees = np.degrees(angle)
   ```

En este caso, el ángulo diedro calculado es de **-25.67º**. Este valor negativo indica que el ángulo de torsión está en el sentido contrario a las agujas del reloj cuando se observa desde el vector definido por los átomos N y CA hacia el vector definido por los átomos C y O.
Un ángulo diedro de **-25.67º** sugiere una ligera torsión que podría estar asociada con la formación de estructuras secundarias específicas, aunque no es un valor característico de una hélice alfa ni una hoja beta, que suelen tener ángulos más definidos. Este resultado destaca la flexibilidad conformacional del residuo inicial en la cadena A de la hemoglobina, un aspecto importante para su función biológica y estabilidad estructural.

### (c) Calcula el centro de masas de la estructura de la hemoglobina. 

El **centro de masas (center of mass)** de una estructura molecular se calcula como el promedio ponderado de las coordenadas de todos los átomos en la estructura, donde el peso de cada átomo es su masa. La fórmula es la siguiente:

$$
\text{Centro de masas (COM)} = \frac{\sum (m_i \cdot r_i)}{\sum m_i}
$$

Donde:
- $m_i$ es la masa del átomo $i$ (generalmente el peso atómico de su elemento).
- $r_i$ es el vector de coordenadas (x, y, z) del átomo $i$.
- La suma se realiza sobre todos los átomos de la estructura.

Información de las masas atómicas: [Disponible aquí](https://www.lenntech.es/periodica/masa/masa-atomica.htm)

1. **Identificación de los elementos presentes en la estructura**  
   Se recorren todos los átomos de la estructura para identificar los elementos únicos presentes. Esto permite asociar las masas atómicas adecuadas a cada elemento:

   ```python
   unique_elements = set()
   for atom in structure.get_atoms():
       element = atom.element.strip().upper()
       unique_elements.add(element)
   ```

   En el caso de la hemoglobina, los elementos identificados son: `'FE'`(Hierro), `'O'`(Oxígeno), `'S'`(Sulfuro), `'N'`(Nitrógeno) y `'C'`(Carbono).

2. **Cálculo de las coordenadas ponderadas por la masa**  
   Para cada átomo, se extraen sus coordenadas ($x, y, z$)) y se ponderan por su masa atómica. Al mismo tiempo, se acumula la masa total de todos los átomos:

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

3. **Obtención del centro de masas**  
   Dividiendo las coordenadas ponderadas ($\sum m_i \cdot x, \sum m_i \cdot y, \sum m_i \cdot z$) por la masa total ($\sum m_i$), se obtiene el centro de masas como un vector tridimensional:

   ```python
   centro_masas = np.array([sum_x, sum_y, sum_z]) / atomic_mass
   ```

4. **Visualización del resultado**  
   El centro de masas calculado se representa gráficamente para resaltar su posición en la estructura de la hemoglobina. En la **Figura 3**, se observa el centro de masas como un punto destacado en la estructura tridimensional:

   <div align="center">
       <img src="results/ej_1_c.png" width="60%">
       <p><b>Figura 3.</b> Centro de masas de la hemoglobina representado como un punto destacado.</p>
   </div>

En el caso de la hemoglobina, se encuentra cerca del núcleo de la estructura, lo cual es consistente con su configuración compacta y su función biológica de transporte de oxígeno en el organismo. Este análisis puede ser extendido a otras proteínas para evaluar la distribución de masa y su relación con propiedades estructurales y funcionales.


