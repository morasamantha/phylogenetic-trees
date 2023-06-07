# Proyecto Final - Genómica Computacional 7118
## Integrantes:
- Mora Abonce Samantha
- Sosa Huarte Victoria Carolina

## Algoritmos y modelos utilizados:

### Alineación de secuencias
- [X] Algoritmo Needleman-Wunsh
- [X] Multiple Sequence Alignment con Star-Alignment* [Link del proyecto](https://github.com/zahrasalarian/Bioinformatics-Mini-Projects)

### Estimación de la distancia genética
- [X] 'Naive' (Sólo contar las diferencias observadas)
- [X] Modelo Jukes-Cantor (JK69)
- [X] Modelo Motoo Kimura (K71)
- [X] Tamura (T92)

### Construcción del árbol
- [X] Neighbor Joining
- [X] UPGMA

### Comparación entre árboles
- [X] Análisis de parsimonía* <- De R

Los métodos marcados con asterísco vienen de una biblioteca/son externos. 

Para usarlo :- )
- Se necesita tener instalado python3, tabulate, networkx, pylocluster, python-newick.

  `$ pip install tabulate`,

  `$ pip install networkx[default]`

  `$ pip install pylocluster`

  `$ pip install newick`

  Pyvolve debería de instalar las dependencias faltantes, pero igual:

  `$ pip install biopython` 

  `$ pip install scipy`
- Desde la carpeta /src ejecutar 

  `$ python3 main.py`

Si se quiere correr la parte de MSA se debe de descomentar unos pedazos de código, además de traer al directorio el main-2 del repositorio ligado arriba y renombrarlo como msa.py. 