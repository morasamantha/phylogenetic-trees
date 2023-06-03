# Proyecto Final - Genómica Computacional 7118
## Integrantes:
- Mora Abonce Samantha
- Sosa Huarte Victoria Carolina

## Algoritmos y modelos utilizados:

### Alineación de secuencias
- [X] Algoritmo Needleman-Wunsh
- [X] Multiple Sequence Alignment con Star-Alignment* [Link del proyecto]{https://github.com/zahrasalarian/Bioinformatics-Mini-Projects}

### Estimación de la distancia genética
- [X] 'Naive' (Sólo contar las diferencias observadas)
- [X] Modelo Jukes-Cantor (JK69)
- [X] Modelo Motoo Kimura (K71)
- [X] Tamura (T92)
- [X] Modelo GTR*

### Construcción del árbol
- [X] Neighbor Joining
- [ ] UPGMA
- [ ] Maximum likehood
- [x] Método Bayesiano (Cadenas de Markov) <- Este es el de Mr Bayes

### Comparación entre árboles
- [ ] Análisis de parsimonía*

Los métodos marcados con asterísco vienen de una biblioteca/son externos. 

Para usarlo :- )
- Se necesita tener instalado python3, tabulate, networkx, pyvolve (con biopython y scipy).

  `$ pip install tabulate`,

  `$ pip install networkx[default]`

  `$ pip install pyvolve`

  Pyvolve debería de instalar las dependencias faltantes, pero igual:

  `$ pip install biopython` 

  `$ pip install scipy`
- Desde la carpeta /src ejecutar 

  `$ python3 main.py`

