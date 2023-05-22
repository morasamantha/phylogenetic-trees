# Proyecto Final - Genómica Computacional 7118
## Integrantes:
- Mora Abonce Samantha
- Sosa Huarte Victoria Carolina

## Algoritmos y modelos utilizados:

### Alineación de secuencias
- [X] Algoritmo Needleman-Wunsh

### Estimación de la distancia genética
- [X] Modelo Jukes-Cantor (JK69)
- [X] Modelo Motoo Kimura (K71)
- [ ] Modelo Hasegawa, Kishino y Yano (HKY85)
- [ ] Modelo GTR 

### Construcción del árbol
- [ ] Neighbor Joining
- [ ] UPGMA
- [ ] WPGMA
- [ ] Maximum likehood
- [x] Método Bayesiano (Cadenas de Markov) <- Este es el de Mr Bayes

(me gustaría añadir un par más aquí para hacer uno con machine learning y remplazar el 2 y 3, pero no he decidido!)

### Comparación entre árboles
- [ ] Análisis de parsimonía 

## Pendientes extras
- [ ] Subir los FASTA.
- [ ] Desde el main empezar a leer los archivos y usar los métodos que ya están. Lo que deberían hacer hasta ahora es:
  - [ ] Alinear todas las secuencias
  - [ ] Las matrices de distancia con cada módelo, cada una creará un .txt.


Para usarlo
- Se necesita tener instalado python3, tabulate (instalar con pip).
- Desde la carpeta principal ejecutar python3 src/main.py
